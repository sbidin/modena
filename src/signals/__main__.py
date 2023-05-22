#!/usr/bin/env python

"""Runs the main application."""

import logging
import os
import sys
from pathlib import Path

import click

from signals.run import run_on_datasets

try:
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
except ValueError:
    logging.basicConfig(level="INFO")
finally:
    log = logging.getLogger("signals")


@click.command()
@click.argument("dataset1")
@click.argument("dataset2")
@click.option("-c", "--min-coverage", type=int, default=5)
@click.option("-r", "--resample", type=int, default=None)
@click.option("--no-distance-sum", is_flag=True, type=bool, default=False)
@click.option("-o", "--out", default="-")
def _main(
        dataset1: str,
        dataset2: str,
        min_coverage: int,
        resample: int | None,
        no_distance_sum: bool | None,
        out: str) \
        -> None:
    """Run the main application."""
    xs_path, ys_path = Path(dataset1), Path(dataset2)
    assert xs_path.exists(), f"no such path exists: {xs_path}"
    assert ys_path.exists(), f"no such path exists: {ys_path}"
    out = sys.stdout if out == "-" else Path(out).open("w")
    run_on_datasets(
        xs_path,
        ys_path,
        min_coverage,
        resample,
        not no_distance_sum,
        out)


if __name__ == "__main__":
    _main()
