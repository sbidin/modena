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
@click.option("--out", default="-")
@click.option("--out-format", default="bedgraph")
@click.option("--min-coverage", default=3)
def _main(
        dataset1: str,
        dataset2: str,
        out: str,
        out_format: str,
        min_coverage: int) \
        -> None:
    """Run the main application."""
    xs_path, ys_path = Path(dataset1), Path(dataset2)
    assert xs_path.exists(), f"no such path exists: {xs_path}"
    assert ys_path.exists(), f"no such path exists: {ys_path}"
    assert out_format in ("bedgraph", "csv"), f"bad format: {out_format}"
    out = sys.stdout if out == "-" else Path(out).open("w")
    run_on_datasets(xs_path, ys_path, min_coverage, out, out_format)


if __name__ == "__main__":
    _main()
