#!/usr/bin/env python

"""Detect epigenetic and epitranscriptomic modifications."""

import collections
import contextlib
import math
import os
from collections import deque
from collections.abc import Iterator
from dataclasses import dataclass
from timeit import default_timer
from typing import Any, TextIO

import astropy.stats
import click
import kmeans1d
import numpy as np
import pyslow5
from cachetools import LRUCache, cached
from loguru import logger
from scipy.stats import median_abs_deviation


def distsum(
        stats: Iterator[tuple[int, float]],
        window_size: int = 5) \
        -> Iterator[tuple[int, float]]:
    """Assign each position a distance that is the sum of its window."""
    window = deque(maxlen=window_size)

    # Starting out, populate the entire window.
    for _ in range(window_size):
        stat = next(stats, None)
        if stat is not None:
            window.append(stat)

    # No window? We're done. This happens if filters are too strict.
    if not window:
        return

    # Handle the prefix first as their left neighbours are incomplete.
    for i in range(window_size // 2):
        yield from _distsum_at(window, i)

    # As the rest are fetched, strip out the prefix.
    for stat in stats:
        yield from _distsum_at(window, window_size // 2)
        window.popleft()
        window.append(stat)

    # Handle the suffix, with incomplete right neighbours.
    for i in range(window_size // 2, window_size):
        yield from _distsum_at(window, i)


def _distsum_at(
        window: deque[tuple[int, float]],
        i: int) \
        -> Iterator[tuple[int, float]]:
    assert 0 <= i < len(window)
    span = window.maxlen // 2
    sum_dists = 0
    for j in range(i - span, i + span + 1):
        if 0 <= j < len(window) and abs(window[i][0] - window[j][0]) <= span:
            sum_dists += window[j][1]
    yield (window[i][0], sum_dists, *window[i][2:])


# def position_within_bounds(
#         pos: int,
#         config: Config,
#         with_window: bool=False) \
#         -> bool:
#     """Return whether a position satisfies various filters."""
#     pos += 1 # Perform from/to bounds checks in 1-based indexing.

#     # Check left bound.
#     if config.from_position is not None:
#         left = config.from_position
#         if with_window:
#             left -= config.WINDOW_SIZE // 2
#         if pos < left:
#             return False

#     # Check right bound.
#     if config.to_position is not None:
#         right = config.to_position
#         if with_window:
#             right += config.WINDOW_SIZE // 2
#         if pos > right:
#             return False

#     return True


@contextlib.contextmanager
def timer(message: str) -> None:
    """Output how long a block of code took to execute."""
    time_from = default_timer()
    yield
    seconds = default_timer() - time_from
    logger.info(f"{message} took {seconds:.3f} seconds")


def get_read_ids(slow_path: str) -> set[str]:
    """Get all read IDs from the given SLOW5 file."""
    logger.info("opening and indexing", slow_path)
    f = pyslow5.Open(slow_path, "r")
    ids, num_reads = f.get_read_ids()
    logger.info("found", num_reads, "reads")
    f.close()
    return set(ids)


def index_reads_by_position(idx_path: str, ids: set[str]) -> dict[int, set[tuple[str, int, int]]]:
    """Create an index of reads by position."""
    idx = collections.defaultdict(set)
    uniq_reads = set()
    num_squiggles = 0
    num_missing = 0
    num_skipped = 0
    min_position = math.inf
    max_position = 0

    logger.info(f"indexing squiggles {idx_path} by position")
    with open(idx_path) as tsv:
        tsv.readline()  # discard header
        for line_raw in tsv:
            line = line_raw.strip().split("\t")
            read_id, pos, start, end = line

            # Skip loading ids that are not in the slow5 file to begin with.
            if read_id not in ids:
                num_skipped += 1
                continue

            if start == "." or end == ".":
                num_missing += 1
                continue

            # TODO: Handle DNA vs RNA differences.
            pos, start, end = int(pos), int(start), int(end)
            assert start < end
            idx[pos].add((read_id, start, end))
            uniq_reads.add(read_id)
            num_squiggles += 1

            min_position = min(min_position, pos)
            max_position = max(max_position, pos)

    logger.info("min position is", min_position)
    logger.info("max position is", max_position)
    logger.info("indexed", len(idx), "unique positions")
    logger.info("indexed", len(uniq_reads), "unique reads")
    logger.info("indexed", num_squiggles, "total squiggles")
    logger.info("skipped", num_missing, "squiggles due to missing signal markers")
    logger.info("skipped", num_skipped, "reads due to not having signal data")
    return idx, min_position, max_position


@cached(cache=LRUCache(maxsize=1e9//4, getsizeof=len), key=lambda _, read_id: read_id)
def get_signal_by_read_id(f, read_id: str) -> np.ndarray:
    """Get a read's full signal."""
    read = f.get_read(read_id, pA=True, aux=None)
    signal = read["signal"]
    scale = median_abs_deviation(signal)
    signal -= np.median(signal)
    signal /= scale
    #signal = np.flip(signal)  # TODO: Don't do this if DNA.
    # TODO: Check if works for DNA without flipping, it should.
    return signal


def get_squiggle(f, read_id: str, a: int, b: int) -> np.ndarray:
    """Get a single squiggle between points a and b."""
    assert a <= b
    signal = get_signal_by_read_id(f, read_id)
    # TODO: Handle DNA vs RNA differences.
    signal = signal[a:b]
    signal = np.random.choice(signal, size=30, replace=True)  # TODO: Configurable size!
    return signal


def concat_signals_by_position(f, idx, pos: int) -> tuple[np.ndarray, int]:
    """Get a concatenated signal for all reads of a given position."""
    signals = []

    for squig in idx[pos]:
        signals.extend(get_squiggle(f, *squig))

    return np.array(signals), len(idx[pos])


def kuiper(x: np.ndarray, y: np.ndarray) -> float:
    """Return the Kuiper statistic for two given samples."""
    np.ndarray.sort(x)
    np.ndarray.sort(y)
    dist = astropy.stats._stats.ks_2samp(x.data, y.data)
    return dist


def print_bed_line(
        # chromosome: str,
        # strand: str,
        pos: int,
        dist: float,
        coverage: int,
        out_file: TextIO) \
        -> None:
    """Emit a single line in the bedMethly format."""
    out = out_file.write

    # The first 11 columns are defined by the bedMethyl format.
    out("_ ") # TODO: f"{chromosome} ") # col 1, reference chromosome
    out(f"{pos + 1} ") # col 2, position from
    out(f"{pos + 2} ") # col 3, position to
    out("_ ") # col 4, name of item
    out("_ ") # col 5, score from 1 to 1000, capped number of reads TODO
    out("_ ") # TODO: f"{strand} ") # col 6, strand
    out("_ ") # col 7, start of where display should be thick
    out("_ ") # col 8, end of where display should be thick
    out("_ ") # col 9, color value
    out(f"{coverage} ") # col 10, coverage, or number of reads
    out("_ ") # col 11, percentage of reads that show methylation at this position TODO

    # All further columns are custom extensions of the format.
    out(f"{dist:.5f} ") # col 12, kuiper distance
    out("\n")


def trim_index(idx, min_pos: int, max_pos: int) -> None:
    """Trim the index to contain only the given span of positions."""
    size = len(idx)
    for pos in list(idx):
        if pos < min_pos or pos > max_pos:
            del idx[pos]

    logger.info("trimmed out", size - len(idx), "positions")


def label_output(out_path: str) -> None:
    """Assign positive and negative labels to an annotated BED file.

    This will read the output BED file and overwrite it with a labeled version.
    """
    # Read entire file.
    with open(out_path) as f:
        lines = [line.split() for line in f.read().splitlines()]
        if not lines:
            logger.warning("output file is empty, no clustering performed")
            return

    # Assign labels.
    scores = [float(line[11]) for line in lines]
    labels, _ = kmeans1d.cluster(scores, 2)
    assert len(scores) == len(labels)

    # Kmeans1d gives us labels of 0 and 1. Treat those assigned to larger scores as positive.
    pos_label = labels[scores.index(max(scores))]

    # Output same file, but with assigned labels, to a different file.
    output_annotated = str(out_path) + ".modena-in-progress"
    with open(output_annotated, "w") as f:
        for line, label in zip(lines, labels, strict=False):
            line.append("pos" if label == pos_label else "neg")
            f.write(" ".join(line))
            f.write("\n")

    # Overwrite the original output with the assigned-labels one.
    os.replace(output_annotated, out_path)


@dataclass
class Config:
    """A bundle of all application configuration values."""
    xs: list[str]
    ys: list[str]
    position: str
    coverage: int
    out: str
    resample: int

    @staticmethod
    def build(kwargs: dict[str, Any]) -> "Config":
        """Create a `Config` given `click` commandline arguments."""
        # TODO: Check config for validity, transform values.
        config = Config(**kwargs)
        return config


def run(conf: Config)  -> None:
    """Run a Modena comparison between two datasets."""
    # Load first dataset.
    raw_x, idx_x = conf.xs
    ids_x = get_read_ids(raw_x)
    idx_x, min_pos_x, max_pos_x = index_reads_by_position(idx_x, ids_x)
    del ids_x

    # Load second dataset.
    raw_y, idx_y = conf.ys
    ids_y = get_read_ids(raw_y)
    idx_y, min_pos_y, max_pos_y = index_reads_by_position(idx_y, ids_y)
    del ids_y

    # Find overlapping positions.
    min_pos = max(min_pos_x, min_pos_y)
    max_pos = min(max_pos_x, max_pos_y)
    logger.info("comparing positions between", min_pos, "and", max_pos)

    # Remove non-overlapping positions from indexes.
    logger.info("trimming index for", raw_x)
    trim_index(idx_x, min_pos, max_pos)
    logger.info("trimming index for", raw_y)
    trim_index(idx_y, min_pos, max_pos)

    # Run the comparison.
    f_x = pyslow5.Open(raw_x, "r")
    f_y = pyslow5.Open(raw_y, "r")
    f_out = open(conf.out, "w")

    def results_at_positions():
        for pos in range(min_pos, max_pos + 1):
            if (pos - min_pos) % 1000 == 0:
                logger.info(f"comparing position {pos}...")

            if pos not in idx_x or pos not in idx_y:
                logger.info("position", pos, "skipped, no overlap")
                if pos in idx_x:
                    del idx_x[pos]
                if pos in idx_y:
                    del idx_y[pos]
                continue

            signals_x, coverage_x = concat_signals_by_position(f_x, idx_x, pos)
            signals_y, coverage_y = concat_signals_by_position(f_y, idx_y, pos)

            min_coverage = 75  # TODO: Configurable!
            if coverage_x >= min_coverage and coverage_y >= min_coverage:
                dist = kuiper(signals_x, signals_y)
                yield pos, dist, coverage_x + coverage_y

            del idx_x[pos]
            del idx_y[pos]

    results = results_at_positions()
    results = (r for r in results if r is not None)
    results = distsum(results)

    for r in results:
        print_bed_line(*r, f_out)

    assert len(idx_x) == 0
    assert len(idx_y) == 0

    f_y.close()
    f_x.close()
    f_out.close()

    logger.info("clustering output, labeling as positive or negative")
    label_output(conf.out)
    logger.info("success, results are in", conf.out)


@click.command()
@click.option("-1", "xs", type=str, nargs=2, required=True)
@click.option("-2", "ys", type=str, nargs=2, required=True)
@click.option("-p", "--position", type=str, default="-")
@click.option("-r", "--resample", type=int, default=15)
@click.option("-c", "--coverage", type=int, default=20)
@click.option("-o", "--out", type=str, required=True)
def cli(**kwargs: dict[str, str | int]) -> None:
    """Collect command-line arguments into a `Config` and run Modena."""
    run(Config.build(kwargs))


if __name__ == "__main__":
    cli()
