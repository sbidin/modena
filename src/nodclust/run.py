"""Exposes the `run_on_datasets` function."""

import hashlib
import logging
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import List, TextIO, Tuple

import astropy.stats
import kmeans1d
import numpy as np

from nodclust.concat import Signal, concat_pairs
from nodclust.config import Config
from nodclust.distsum import distsum
from nodclust.fast5 import Fast5
from nodclust.helpers import position_within_bounds, timer

log = logging.getLogger("nodclust")


def compare_datasets(config: Config) -> None:
    """Run the application on provided datasets, outputing to a file."""
    # Index the datasets by overlapping ranges of positions.
    with timer("ordering datasets by size"):
        xs_path, ys_path = _order_paths_by_size(config.dataset1, config.dataset2)

    with timer("indexing datasets"):
        xs, ys = _index_datasets(xs_path, ys_path, config)

    # Some metadata gets output once per every line.
    chromosome, strand = xs[0].chromosome, xs[0].strand

    # Get Kuiper values for each position signal pair.
    pairs = concat_pairs(xs, ys, config)
    stats = map(_kuiper, pairs)
    if not config.no_distance_sum:
        stats = distsum(stats, config)

    # Output the resulting BED file.
    with timer("calculating position distances"), open(config.output_bed, "w") as f:
        for pos, dist in stats:
            if position_within_bounds(pos, config):
                _emit_line_bed_methyl(chromosome, strand, pos, dist, f)

    # Overwrite it with an annotated BED file.
    with timer("clustering output"):
        _cluster_output(config)


def _order_paths_by_size(a: Path, b: Path) -> Tuple[Path, Path]:
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
        log.debug(f"size of {path} is {size} bytes")
        return size
    except:
        # We don't have access to size information so just make it consistent.
        log.debug(f"could not determine size of {path}")
        return hashlib.md5(str(path).encode()).digest()[0] + 1


def _index_datasets(
        xs_path: Path,
        ys_path: Path,
        config: Config) \
        -> Tuple[List[Fast5], List[Fast5]]:
    """Preload, filter and sort relevant dataset subsets."""
    # The first (and smaller) dataset's basic data is read fully and sorted by
    # position. The second dataset is read in a streaming fashion, unordered.
    # We skip any file from the second dataset whose positions do not overlap
    # with those of the first dataset.
    xs = Fast5.from_path(xs_path, config)
    xs = sorted(xs, key=lambda x: (x.start, -x.end))

    if not xs:
        log.error(f"no files satisfying filters found in {xs_path}")
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


def _kuiper(t: Tuple[Signal, Signal]) -> Tuple[int, float]:
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
        out_file: TextIO) \
        -> None:
    """Emit a single line in the bedMethly format."""
    out = out_file.write

    # The first 11 columns are defined by the bedMethyl format.
    out(f"{chromosome} ") # col 1, reference chromosome
    out(f"{pos + 1} ") # col 2, position from
    out(f"{pos + 2} ") # col 3, position to
    out("_ ") # col 4, name of item
    out("_ ") # col 5, score from 1 to 1000, capped number of reads TODO
    out(f"{strand} ") # col 6, strand
    out("_ ") # col 7, start of where display should be thick
    out("_ ") # col 8, end of where display should be thick
    out("_ ") # col 9, color value
    out("_ ") # col 10, coverage, or number of reads TODO
    out("_ ") # col 11, percentage of reads that show methylation at this position TODO

    # All further columns are custom extensions of the format.
    out(f"{dist:.5f} ") # col 12, kuiper distance
    out("\n")


def _cluster_output(config: Config) -> None:
    """Assign positive and negative labels to an annotated BED file.

    This will read the output BED file and overwrite it with a labeled version.
    """
    # Read entire file.
    with config.output_bed.open() as f:
        lines = [line.split() for line in f.read().splitlines()]
        if not lines:
            log.warning("output file is empty, no clustering performed")
            return

    # Assign labels.
    scores = [float(line[11]) for line in lines]
    labels, _ = kmeans1d.cluster(scores, 2)
    assert len(scores) == len(labels)

    # Kmeans1d gives us labels of 0 and 1. Treat those assigned to larger scores as positive.
    pos_label = labels[scores.index(max(scores))]

    # Output same file, but with assigned labels, to a different file.
    output_annotated = str(config.output_bed) + ".cluster-in-progress"
    with open(output_annotated, "w") as f:
        for line, label in zip(lines, labels):
            line.append("pos" if label == pos_label else "neg")
            f.write(" ".join(line))
            f.write("\n")

    # Overwrite the original output with the assigned-labels one.
    os.rename(output_annotated, config.output_bed)
