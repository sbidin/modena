"""Exposes the `run_on_datasets` function."""

import dataclasses
import hashlib
import logging
import re
import subprocess
import sys
from pathlib import Path

import astropy.stats
import numpy as np

from nodclust.concat import Signal, concat_pairs
from nodclust.config import Config
from nodclust.distsum import distsum
from nodclust.fast5 import Fast5
from nodclust.helpers import position_within_bounds

log = logging.getLogger("nodclust")


def compare_datasets(config: Config) -> None:
    """Run the application on provided datasets, outputing to a file."""
    # Index the datasets by overlapping ranges of positions.
    xs_path, ys_path = _order_paths_by_size(config.dataset1, config.dataset2)
    xs, ys = _index_datasets(xs_path, ys_path, config)

    # Some metadata gets output once per every line.
    chromosome, strand = xs[0].chromosome, xs[0].strand

    # Get Kuiper values for each position signal pair.
    pairs = concat_pairs(xs, ys, config)
    stats = map(_kuiper, pairs)
    if not config.no_distance_sum:
        stats = distsum(stats, config)

    # Output the resulting BED file.
    for pos, dist in stats:
        if position_within_bounds(pos, config):
            _emit_line_bed_methyl(chromosome, strand, pos, dist, config)


def _order_paths_by_size(a: Path, b: Path) -> tuple[Path, Path]:
    """Reorder given paths by size in ascending order."""
    datasets = [a, b]
    datasets.sort(key=_size_at_path)
    return datasets


def _size_at_path(path: Path) -> int:
    """Get the size in bytes of a file or directory at `path`."""
    try:
        size = subprocess.check_output(["du", "-sb", str(path)])
        size = int(size.split()[0].decode("utf-8"))
        assert size > 0, "size is zero"
        return size
    except:
        # We don't have access to size information so just make it consistent.
        return hashlib.md5(str(path).encode()).digest()[0] + 1


def _index_datasets(
        xs_path: Path,
        ys_path: Path,
        config: Config) \
        -> tuple[list[Fast5], list[Fast5]]:
    """Preload, filter and sort relevant dataset subsets."""
    # The first (and smaller) dataset's basic data is read fully and sorted by
    # position. The second dataset is read in a streaming fashion, unordered.
    # We skip any file from the second dataset whose positions do not overlap
    # with those of the first dataset.
    xs = Fast5.from_path(xs_path, config)
    xs = sorted(xs, key=lambda x: (x.start, -x.end))

    if not xs:
        log.error(f"no valid files in {xs_path}")
        sys.exit(1)

    ys = []
    ys_orig_len = 0

    # Infer second dataset's config to match the first.
    config = Config(**config.__dict__) # Make a copy so as to not destroy original.
    config.acid = xs[0].acid
    config.chromosome = re.compile(fr"^{xs[0].chromosome}$")
    config.strand = xs[0].strand

    for y in Fast5.from_path(ys_path, config):
        ys_orig_len += 1
        if y.overlap(xs):
            ys.append(y)

    if not ys:
        log.error(f"no selected files from {xs_path} overlap those in {ys_path}")
        sys.exit(1)

    # We then reverse the process: there are files from the first dataset that
    # never overlap with any of the ones from the second dataset. Get rid of
    # them since their positions will never be matched anyway.
    ys.sort(key=lambda y: (y.start, -y.end))
    xs = [x for x in xs if x.overlap(ys)]
    return xs, ys


def _kuiper(t: tuple[Signal, Signal]) -> tuple[int, float]:
    """Return the Kuiper statistic for two given samples."""
    x, y = t
    np.ndarray.sort(x.data)
    np.ndarray.sort(y.data)
    d = astropy.stats._stats.ks_2samp(x.data, y.data)
    return x.position, d


def _emit_line_bed_methyl(
        chromosome: str,
        strand: str,
        pos: int,
        dist: float,
        config: Config) \
        -> None:
    """Emit a single line in the bedMethly format."""
    out = config.out.write

    # The first 11 columns are defined by the bedMethyl format.
    out(f"{chromosome} ")  # col 1, reference chromosome
    out(f"{pos + 1} ")  # col 2, position from
    out(f"{pos + 2} ")  # col 3, position to
    out("_ ")  # col 4, name of item
    out("_ ")  # col 5, score from 1 to 1000, capped number of reads TODO
    out(f"{strand} ") # col 6, strand
    out("_ ")  # col 7, start of where display should be thick
    out("_ ")  # col 8, end of where display should be thick
    out("_ ")  # col 9, color value
    out("_ ")  # col 10, coverage, or number of reads TODO
    out("_ ")  # col 11, percentage of reads that show methylation at this position TODO

    # All further columns are custom extensions of the format.
    out(f"{dist:.5f}")  # col 12, kuiper distance
    out("\n")
