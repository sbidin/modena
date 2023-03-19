#!/usr/bin/env python

import logging
import os
import subprocess
from pathlib import Path

import click

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
def _main(dataset1: str, dataset2: str) -> None:
    """Run the main application."""
    xs_path, ys_path = Path(dataset1), Path(dataset2)
    assert xs_path.exists(), f"no such path exists: {xs_path}"
    assert ys_path.exists(), f"no such path exists: {ys_path}"
    xs_path, ys_path = reorder_datasets_by_size(xs_path, ys_path)

    # The first (and smaller) dataset's basic data is read fully and sorted by
    # position. The second dataset is read in a streaming fashion, unordered.
    # We skip any file from the second dataset whose positions do not overlap
    # with those of the first dataset.
    log.info(f"dataset1: loading: {xs_path}")
    xs = sorted(Fast5.from_path(xs_path), key=lambda x: (x.start, -x.size))
    xs_orig_len = len(xs)
    log.info(f"dataset1: loaded files: {xs_orig_len}")

    ys = []
    ys_orig_len = 0
    log.info(f"dataset2: loading: {ys_path}")
    for y in Fast5.from_path(ys_path):
        ys_orig_len += 1
        if y.overlap(xs):
            ys.append(y)
    log.info(f"dataset2: loaded files: {ys_orig_len}")
    log.info(f"dataset2: filtered down to: {len(ys)}")

    # We then reverse the process: there are files from the first dataset that
    # never overlap with any of the ones from the second dataset. Get rid of
    # them since their positions will never be matched anyway.
    log.info(f"dataset1: filtering: {xs_path}")
    ys.sort(key=lambda y: (y.start, -y.size))
    xs = [x for x in xs if x.overlap(ys)]
    log.info(f"dataset1: filtered down to: {len(xs)}")


if __name__ == "__main__":
    try:
        _main()
    except AssertionError as err:
        log.error(f"assertion failed: {err}")
