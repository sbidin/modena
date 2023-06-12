"""Exposes the `run_on_datasets` function."""

import logging
import sys
from pathlib import Path
from typing import TextIO

import astropy.stats
import numpy as np

from nodclust.concat import Signal, concat_pairs
from nodclust.distsum import distsum
from nodclust.fast5 import Fast5

log = logging.getLogger("nodclust")


def emit_line_bed_methyl(
        chrom: str,
        strand: str,
        pos: int,
        dist: float,
        out: TextIO) \
        -> None:
    """Emit a single line in the bedMethly format."""
    # The first 11 columns are defined by the bedMethyl format.
    out.write(f"{chrom} ")  # col 1, reference chromosome
    out.write(f"{pos + 1} ")  # col 2, position from
    out.write(f"{pos + 2} ")  # col 3, position to
    out.write("_ ")  # col 4, name of item
    out.write("_ ")  # col 5, score from 1 to 1000, capped number of reads TODO
    out.write(f"{strand} ") # col 6, strand
    out.write("_ ")  # col 7, start of where display should be thick
    out.write("_ ")  # col 8, end of where display should be thick
    out.write("_ ")  # col 9, color value
    out.write("_ ")  # col 10, coverage, or number of reads TODO
    out.write("_ ")  # col 11, percentage of reads that show methylation at this position TODO

    # All further columns are custom extensions of the format.
    out.write(f"{dist:.5f}")  # col 12, kuiper distance
    out.write("\n")


def run_on_datasets(
        xs_path: Path,
        ys_path: Path,
        type: str,
        strand: str | None,
        chrom: str | None,
        force_type: bool,
        min_coverage: int,
        resample: int,
        distance_sum: bool,
        out: TextIO,
        random_seed: int | None) \
        -> None:
    """Run the application on provided datasets, outputing to a file."""
    if random_seed is not None:
        np.random.seed(random_seed)

    xs, ys = index_datasets(xs_path, ys_path, type, strand, chrom, force_type)

    # Some metadata gets output once per every line. We store it here since
    # otherwise it'll get deleted during the process below.
    chrom, strand = xs[0].chrom, xs[0].strand

    # Get dist-summed Kuiper statistics for pairs of positions' signals.
    pairs = concat_pairs(xs, ys, min_coverage, resample)
    stats = map(kuiper, pairs)
    if distance_sum:
        stats = distsum(stats, 5)

    # Output the BED file.
    for pos, dist in stats:
        emit_line_bed_methyl(chrom, strand, pos, dist, out)


def index_datasets(
        xs_path: Path,
        ys_path: Path,
        type: str,
        strand: str | None,
        chrom: str | None,
        force_type: bool) \
        -> tuple[list[Fast5], list[Fast5]]:
    """Preload, filter and sort relevant dataset subsets.

    Accepts strand and chromosome filters.
    """
    # The first (and smaller) dataset's basic data is read fully and sorted by
    # position. The second dataset is read in a streaming fashion, unordered.
    # We skip any file from the second dataset whose positions do not overlap
    # with those of the first dataset.
    xs = Fast5.from_path(xs_path, type, strand, chrom, force_type)
    xs = sorted(xs, key=lambda x: (x.start, -x.end))

    if not xs:
        log.error(f"no valid files in {xs_path}")
        sys.exit(1)

    ys = []
    ys_orig_len = 0
    for y in Fast5.from_path(ys_path, xs[0].type, xs[0].strand, fr"^{xs[0].chrom}$", force_type):
        ys_orig_len += 1
        if y.overlap(xs):
            ys.append(y)

    if not ys:
        log.error(f"no valid files with chrom {xs[0].chrom} in {ys_path}")
        sys.exit(1)

    # We then reverse the process: there are files from the first dataset that
    # never overlap with any of the ones from the second dataset. Get rid of
    # them since their positions will never be matched anyway.
    ys.sort(key=lambda y: (y.start, -y.end))
    xs = [x for x in xs if x.overlap(ys)]
    return xs, ys


def kuiper(t: tuple[Signal, Signal]) -> tuple[int, float]:
    """Return the Kuiper statistic for two given samples."""
    x, y = t
    np.ndarray.sort(x.data)
    np.ndarray.sort(y.data)
    d = astropy.stats._stats.ks_2samp(x.data, y.data)
    return x.position, d
