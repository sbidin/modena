#!/usr/bin/env python

import logging
import os
import subprocess
from pathlib import Path

import astropy.stats
import click

from signals import concat
from signals.fast5 import Fast5

try:
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
except ValueError:
    logging.basicConfig(level="INFO")
finally:
    log = logging.getLogger("signals")


def size_at_path(path: Path) -> int | None:
    """Get the size in bytes of a file or directory at path."""
    try:
        size = subprocess.check_output(["du", "-sb", str(path)])
        size = int(size.split()[0].decode("utf-8"))
        assert size > 0, "size is zero"
        log.debug(f"size of {path} is {size}")
        return size
    except Exception as err:
        log.warning(f"could not estimate size at path {path}: {err}")
        return 0


def reorder_datasets_by_size(a: Path, b: Path) -> tuple[Path, Path]:
    """Reorder given datasets by size; we process the smaller one first."""
    datasets = [a, b]
    datasets.sort(key=size_at_path)
    return tuple(datasets)


@click.command()
@click.argument("dataset1")
@click.argument("dataset2")
@click.option("--min-coverage", default=3)
def _main(dataset1: str, dataset2: str, min_coverage: int) -> None:
    """Run the main application."""
    xs_path, ys_path = Path(dataset1), Path(dataset2)
    assert xs_path.exists(), f"no such path exists: {xs_path}"
    assert ys_path.exists(), f"no such path exists: {ys_path}"
    xs_path, ys_path = reorder_datasets_by_size(xs_path, ys_path)

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

    for x, y in concat.concat_pairs(xs, ys, min_coverage):
        print(
            x.position,
            x.coverage,
            y.coverage,
            *astropy.stats.kuiper_two(x.data, y.data))


if __name__ == "__main__":
    try:
        _main()
    except AssertionError as err:
        log.error(f"assertion failed: {err}")
