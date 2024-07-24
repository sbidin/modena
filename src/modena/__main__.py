#!/usr/bin/env python

"""Detect epigenetic and epitranscriptomic modifications."""

import collections
import math
import os
import sys
from collections import deque
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO, TextIO

import astropy.stats
import click
import kmeans1d
import numpy as np
import pyslow5
from cachetools import LRUCache, cached
from loguru import logger
from scipy.stats import median_abs_deviation

# Configure logger to skip outputting code location.
logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> <level>{level}</level> {message}",
    level=os.getenv("LOGLEVEL", "INFO"),
)

@dataclass
class Config:
    """A bundle of all application configuration values."""

    blow5_x: Path  # 1st dataset's blow5/slow5 file.
    blow5_y: Path  # 2nd dataset's blow5/slow5 file.
    squig_x: Path  # 1st dataset's f5c resquiggle tsv file.
    squig_y: Path  # 2nd dataset's f5c resquiggle tsv file.
    coverage: int  # Minimum coverage a position needs to have, per-dataset.
    resample: int  # How many times to resample each position's signal.
    output: Path   # Path of output labeled tsv file.

    @staticmethod
    def build(args: dict[str, Any]) -> "Config":
        """Create a `Config` given `click` commandline arguments."""
        # Check file inputs.
        blow5_x, squig_x = Path(args["blow5_x"]), Path(args["squig_x"])
        blow5_y, squig_y = Path(args["blow5_y"]), Path(args["squig_y"])
        for path in (blow5_x, blow5_y, squig_x, squig_y):
            if not path.exists():
                logger.error(f"no such file: {path}")
                sys.exit(1)

        # Check minimum coverage.
        coverage = args["coverage"]
        if coverage <= 0:
            logger.error("coverage needs to be at least 1")
            sys.exit(1)

        # Check resample size.
        resample = args["resample"]
        if resample < 0:
            logger.error("resample needs to be at least 0 (disabled) or higher")
            sys.exit(1)
        if resample == 0:
            resample = None

        # Check output path.
        output = Path(args["output"])
        if output.is_dir():
            logger.error("output path must point to file, not directory")
            sys.exit(1)
        if output.exists():
            logger.warning("output path already exists, will overwrite")

        return Config(
            blow5_x=blow5_x,
            blow5_y=blow5_y,
            squig_x=squig_x,
            squig_y=squig_y,
            coverage=coverage,
            resample=resample,
            output=output,
        )


@click.command()
@click.argument("blow5_x")
@click.argument("squig_x")
@click.argument("blow5_y")
@click.argument("squig_y")
@click.option("-r", "--resample", type=int, default=15, help="resample signals to this size (default 15)")
@click.option("-c", "--coverage", type=int, default=20, help="ignore positions below this coverage (default 20)")
@click.option("-o", "--output", type=str, required=True, help="output file path")
def cli(**kwargs: dict[str, str | int]) -> None:
    """Detect modifications and output per-position results."""
    run(Config.build(kwargs))


def run(conf: Config)  -> None:
    """Detect modifications and output per-position results."""
    # Load and index both datasets.
    ids_x = get_all_read_ids(conf.blow5_x)
    idx_x, min_pos_x, max_pos_x = index_squiggles_by_position(conf.squig_x, ids_x)
    del ids_x
    ids_y = get_all_read_ids(conf.blow5_y)
    idx_y, min_pos_y, max_pos_y = index_squiggles_by_position(conf.squig_y, ids_y)
    del ids_y

    # Find overlapping positions.
    min_pos = max(min_pos_x, min_pos_y)
    max_pos = min(max_pos_x, max_pos_y)
    logger.info(f"processing positions between {min_pos} and {max_pos}")

    # Remove non-overlapping positions from both indexes.
    trim_index(conf.squig_x, idx_x, min_pos, max_pos)
    trim_index(conf.squig_y, idx_y, min_pos, max_pos)

    # Run the comparison.
    f_x = pyslow5.Open(str(conf.blow5_x), "r")
    f_y = pyslow5.Open(str(conf.blow5_y), "r")
    with open(conf.output, "w") as f_out:
        print_header(f_out)
        results = compare_position_pairs(conf, f_x, f_y, idx_x, idx_y, min_pos, max_pos)
        results = distance_summing(results)
        count = 0
        for result in results:
            print_result(*result, f_out)
            count += 1

    logger.info(f"compared {count} total positions")

    # Sanity checking.
    assert len(idx_x) == 0
    assert len(idx_y) == 0

    # Cleanup.
    f_y.close()
    f_x.close()
    f_out.close()

    # Apply kmeans1d clustering.
    cluster_output(conf.output)


def get_all_read_ids(blow5_path: Path) -> set[str]:
    """Get all read IDs from the given BLOW5 file."""
    logger.info(f"indexing dataset {blow5_path}")
    f = pyslow5.Open(str(blow5_path), "r")
    ids, num_reads = f.get_read_ids()
    logger.info(f"found {num_reads} unique reads")
    f.close()
    if not num_reads:
        logger.error(f"no reads found in {blow5_path}")
        sys.exit(1)
    return set(ids)


# An index of positions to (read id, squiggle from, squiggle to) pairs.
PositionIndex = dict[int, set[tuple[str, int, int]]]


def index_squiggles_by_position(
        idx_path: Path,
        read_ids: set[str]) \
        -> PositionIndex:
    """Create an index of squiggles across all reads, by position."""
    logger.info(f"indexing squiggles {idx_path}")
    idx = collections.defaultdict(set)  # The index itself.

    # Various counters for logging purposes.
    uniq_reads = set()
    num_squiggles = 0
    num_missing = 0
    num_skipped = 0
    min_position = math.inf
    max_position = 0

    # Build the index by iterating through all read/position lines.
    with open(idx_path) as tsv:
        tsv.readline()  # Discard the header.

        for line in tsv:
            # Unpack columns.
            read_id, pos, start, end = line.strip().split("\t")

            # Skip loading ids that are not in the blow5 file to begin with.
            if read_id not in read_ids:
                num_skipped += 1
                continue

            # Skip unknown squiggle positions.
            if start == "." or end == ".":
                num_missing += 1
                continue

            pos, start, end = int(pos), int(start), int(end)
            assert start < end, "squiggle start position is after its end"

            # Update index and counters.
            idx[pos].add((read_id, start, end))
            uniq_reads.add(read_id)
            num_squiggles += 1
            min_position = min(min_position, pos)
            max_position = max(max_position, pos)

    # Check and log a bunch of stuff to make sure we have somewhat correct data.
    logger.info(f"minimum found position is {min_position}")
    logger.info(f"maximum found position is {max_position}")
    logger.info(f"indexed {len(idx)} unique positions")
    logger.info(f"indexed {len(uniq_reads)} unique reads")
    if not uniq_reads:
        logger.error(f"no reads found in {idx_path}")
        sys.exit(1)
    logger.info(f"indexed {num_squiggles} total squiggles")
    if num_missing:
        logger.warning(f"skipped {num_missing} squiggles due to missing markers")
    if num_skipped:
        logger.warning(f"skipped {num_skipped} reads due to missing in blow5 file")

    return idx, min_position, max_position


def trim_index(
        idx_path: Path,
        idx: PositionIndex,
        min_pos: int,
        max_pos: int) \
        -> None:
    """Trim the index to contain only the given span of positions."""
    logger.info(f"trimming index for {idx_path}")
    size_old = len(idx)
    for pos in list(idx):
        if pos < min_pos or pos > max_pos:
            del idx[pos]
    logger.info(f"removed {size_old - len(idx)} positions")


# A comparison result: position, Kuiper distance and total coverage.
Result = tuple[int, float, int]


def compare_position_pairs(
        conf: Config,
        f_x: BinaryIO,
        f_y: BinaryIO,
        idx_x: PositionIndex,
        idx_y: PositionIndex,
        min_pos: int,
        max_pos: int) \
        -> Iterator[Result]:
    """Compare signals at same positions from both datasets."""
    for pos in range(min_pos, max_pos + 1):
        # A bit of progress tracking.
        if (pos - min_pos) % 1000 == 0:
            logger.info(f"comparing position {pos + 1}...")

        # Compare signals if filters satisfied.
        if pos not in idx_x or pos not in idx_y:
            logger.debug(f"skipped position {pos + 1} due to no overlap")
        elif len(idx_x[pos]) < conf.coverage or len(idx_y[pos]) < conf.coverage:
            logger.debug(f"skipped position {pos + 1} due to insufficient coverage")
        else:
            signals_x, signals_y = [], []
            for squig in (get_squiggle(conf, f_x, *s) for s in idx_x[pos]):
                signals_x.extend(squig)
            for squig in (get_squiggle(conf, f_y, *s) for s in idx_y[pos]):
                signals_y.extend(squig)
            signals_x, signals_y = np.array(signals_x), np.array(signals_y)
            dist = kuiper(signals_x, signals_y)
            yield pos, dist, len(idx_x[pos]) + len(idx_y[pos])

        # Clear index for the processed position.
        if pos in idx_x:
            del idx_x[pos]
        if pos in idx_y:
            del idx_y[pos]


def distance_summing(
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


@cached(cache=LRUCache(maxsize=1e9, getsizeof=len), key=lambda _, read_id: read_id)
def get_signal_by_read_id(f: BinaryIO, read_id: str) -> np.ndarray:
    """Get a read's full signal, caching the output."""
    read = f.get_read(read_id, pA=True, aux=None)
    signal = read["signal"]
    scale = median_abs_deviation(signal)
    signal -= np.median(signal)
    signal /= scale
    return signal


def get_squiggle(
        conf: Config,
        f: BinaryIO,
        read_id: str,
        a: int,
        b: int) \
        -> np.ndarray:
    """Get a single squiggle between points a and b."""
    assert a <= b
    signal = get_signal_by_read_id(f, read_id)
    signal = signal[a:b]
    signal = np.random.choice(signal, size=conf.resample, replace=True)  # Resampling.
    return signal


def kuiper(x: np.ndarray, y: np.ndarray) -> float:
    """Return the Kuiper statistic for two given samples."""
    np.ndarray.sort(x)
    np.ndarray.sort(y)
    dist = astropy.stats._stats.ks_2samp(x.data, y.data)
    return dist


def print_header(out_file: TextIO) -> None:
    """Print the header for the output file."""
    print("pos\tcov\tdist\tlabel", file=out_file)


def print_result(
        pos: int,
        dist: float,
        coverage: int,
        out_file: TextIO) \
        -> None:
    """Emit a single line for a given result."""
    out = out_file.write
    out(f"{pos + 1}\t") # col 1, position, 1-based
    out(f"{coverage}\t") # col 2, total position coverage, summed across both datasets
    out(f"{dist:.5f}") # col 3, kuiper distance
    # col 4 is added later by clustering and is a pos/neg label
    out("\n")


def cluster_output(out_path: Path) -> None:
    """Assign positive and negative labels to the unlabeled tsv file."""
    logger.info("clustering output into positive/negative via kmeans1d")

    # Read entire file.
    with open(out_path) as f:
        lines = [line.split("\t") for line in f.read().splitlines()]
        lines = lines[1:]  # Skip the header.
        if not lines:
            logger.warning("output file is empty, no clustering performed")
            return

    # Assign labels.
    scores = [float(line[2]) for line in lines]
    labels, _ = kmeans1d.cluster(scores, 2)
    assert len(scores) == len(labels)

    # Kmeans1d gives us labels of 0 and 1. Treat those assigned to larger scores as positive.
    pos_label = labels[scores.index(max(scores))]

    # Output same file, but with assigned labels, to a different file.
    output_annotated = str(out_path) + ".modena-in-progress"
    with open(output_annotated, "w") as f:
        print_header(f)
        for line, label in zip(lines, labels, strict=False):
            line.append("pos" if label == pos_label else "neg")
            f.write("\t".join(line))
            f.write("\n")

    # Overwrite the original output with the assigned-labels one.
    os.replace(output_annotated, out_path)
    logger.info(f"all done, results are stored into {out_path}")


if __name__ == "__main__":
    cli()
