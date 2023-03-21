"""Exposes the `run_on_datasets` function."""

import multiprocessing
from pathlib import Path
from typing import TextIO

import astropy.stats

from signals.concat import Signal, concat_pairs
from signals.fast5 import Fast5
from signals.helpers import order_paths_by_size


def emit_line_csv(x: Signal, dist: float, out: TextIO) -> None:
    """Emit a single value in the CSV format."""
    out.write(f"{x.position},{dist:.5f}\n")


def emit_line_bedgraph(x: Signal, dist: float, out: TextIO) -> None:
    """Emit a single value in the bedGraph format."""
    out.write(f"c1 {x.position} {x.position + 1} {dist:.5f}\n")


def run_on_datasets(
        xs_path: Path,
        ys_path: Path,
        min_coverage: int,
        out: TextIO,
        out_format: str) \
        -> None:
    """Run the application on provided datasets, outputing to a file."""
    xs_path, ys_path, flipped = order_paths_by_size(xs_path, ys_path)
    xs, ys = index_datasets(xs_path, ys_path)
    if flipped: # Undo the order flip so the output doesn't get mirrored.
        xs, ys = ys, xs

    # Format-specific configuration.
    if out_format == "bedgraph":
        out.write("track type=bedGraph\n")
        emit_func = emit_line_bedgraph
    else:
        out.write("position,distance\n")
        emit_func = emit_line_csv

    with multiprocessing.Pool() as pool:
        pairs = concat_pairs(xs, ys, min_coverage)
        for x, _, dist, _ in pool.imap(kuiper, pairs, chunksize=10_000):
            emit_func(x, dist, out)


def index_datasets(xs_path: Path, ys_path: Path) -> tuple[list[Fast5], list[Fast5]]:
    """Preload, filter and sort relevant dataset subsets."""
    # The first (and smaller) dataset's basic data is read fully and sorted by
    # position. The second dataset is read in a streaming fashion, unordered.
    # We skip any file from the second dataset whose positions do not overlap
    # with those of the first dataset.
    xs = sorted(Fast5.from_path(xs_path), key=lambda x: (x.start, -x.end))
    ys = []
    ys_orig_len = 0
    for y in Fast5.from_path(ys_path):
        ys_orig_len += 1
        if y.overlap(xs):
            ys.append(y)

    # We then reverse the process: there are files from the first dataset that
    # never overlap with any of the ones from the second dataset. Get rid of
    # them since their positions will never be matched anyway.
    ys.sort(key=lambda y: (y.start, -y.end))
    xs = [x for x in xs if x.overlap(ys)]
    return xs, ys


def kuiper(t: tuple[Signal, Signal]) -> tuple[Signal, Signal, float, float]:
    """Return the Kuiper statistic for two given samples."""
    x, y = t
    d, p = astropy.stats.kuiper_two(x.data, y.data)
    return x, y, d, p
