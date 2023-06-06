"""Exposes the `run_on_datasets` function."""

from pathlib import Path
from typing import TextIO

import astropy.stats
import numpy as np

from signals.concat import Signal, concat_pairs
from signals.distsum import distsum
from signals.fast5 import Fast5
from signals.helpers import order_paths_by_size


def emit_line_bed_methyl(pos: int, dist: float, out: TextIO) -> None:
    """Emit a single line in the bedMethly format."""
    # Col 1 Reference chromosome or scaffold TODO, currently placeholder
    # Col 4 Name of item TODO
    # Col 5 Score from 0-1000. Capped number of reads TODO
    # Col 6 Strandedness, plus (+), minus (-), or unknown (.) TODO proslijediti iz dataseta
    # Col 7 Start of where display should be thick (start codon) TODO
    # Col 8 End of where display should be thick (stop codon) TODO
    # Col 9 Color value (RGB) TODO
    # Col 10 Coverage, or number of reads TODO suma prvog i drugog dataseta?
    # Col 11 Percentage of reads that show methylation at this position in the genome TODO
    # Col 12+ su arbitrary i u njih stavljamo svoje ekstra things, trenutno distanca
    out.write(f"c1 {pos + 1} {pos + 2} _ _ _ _ _ _ _ _ {dist:.5f}\n")


def run_on_datasets(
        xs_path: Path,
        ys_path: Path,
        strand: str,
        chrom: str,
        min_coverage: int,
        resample: int,
        distance_sum: bool,
        out: TextIO,
        random_seed: int | None) \
        -> None:
    """Run the application on provided datasets, outputing to a file."""
    if random_seed is not None:
        np.random.seed(random_seed)

    xs_path, ys_path, flipped = order_paths_by_size(xs_path, ys_path)
    xs, ys = index_datasets(xs_path, ys_path, strand, chrom)
    if flipped: # Undo the order flip so the output doesn't get mirrored.
        xs, ys = ys, xs

    pairs = concat_pairs(xs, ys, min_coverage, resample)
    stats = map(kuiper, pairs)
    if distance_sum:
        stats = distsum(stats, 5)

    # out.write("track ...") TODO: Do we need headers?
    for pos, dist in stats:
        emit_line_bed_methyl(pos, dist, out)


def index_datasets(
        xs_path: Path,
        ys_path: Path,
        strand: str | None,
        chrom: str | None) \
        -> tuple[list[Fast5], list[Fast5]]:
    """Preload, filter and sort relevant dataset subsets.

    Accepts strand and chromosome filters.
    """
    # The first (and smaller) dataset's basic data is read fully and sorted by
    # position. The second dataset is read in a streaming fashion, unordered.
    # We skip any file from the second dataset whose positions do not overlap
    # with those of the first dataset.
    xs = Fast5.from_path(xs_path, strand, chrom)
    xs = sorted(xs, key=lambda x: (x.start, -x.end))

    ys = []
    ys_orig_len = 0
    for y in Fast5.from_path(ys_path, strand, chrom):
        ys_orig_len += 1
        if y.overlap(xs):
            ys.append(y)

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
