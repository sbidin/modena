#!/usr/bin/env python

"""Runs the main application."""

import logging
import os
import sys
from pathlib import Path
from typing import Any

import click
import jenkspy

from nodclust.config import Config
from nodclust.run import compare_datasets

try:
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
except ValueError:
    logging.basicConfig(level="INFO")
finally:
    log = logging.getLogger("nodclust")


@click.command()
@click.argument("dataset1")
@click.argument("dataset2")
@click.output("output_bed")
@click.option("-a", "--acid", type=str, default="autodetect", help="Filter by acid, dna or rna")
@click.option("-c", "--chromosome", type=str, default=None, help="Filter by chromosome regex")
@click.option("--force-acid", is_flag=True, type=bool, default=False, help="Force read files as specified by --acid")
@click.option("-f", "--from-position", type=int, default=None, help="Filter by minimum position (inclusive)")
@click.option("-m", "--min-coverage", type=int, default=5, help="Filter by minimum coverage (default 5)")
@click.option("--no-distance-sum", is_flag=True, type=bool, default=False, help="Don't sum neighbour position distances")
@click.option("--random-seed", type=int, default=None, help="Force a random seed, for reproducibility")
@click.option("-r", "--resample-size", type=int, default=15, help="Signal resample size; 0 to disable (default 15)")
@click.option("-s", "--strand", type=str, default=None, help="Filter by strand, '+' or '-'")
@click.option("-t", "--to-position", type=int, default=None, help="Filter by maximum position (inclusive)")
def cli(**kwargs: dict[str, Any]) -> None:
    """Compare two datasets & output an annotated bedMethyl file."""
    compare_datasets(Config.from_cmd_args(**kwargs))


@click.command()
@click.argument("bed_file")
@click.option("-o", "--out", default="-", help="Output to a given path (default stdout)")
def label(bed_file: str, out: str) -> None:
    """Assign positive and negative labels to an annotated BED file.

    Note that this process can take a long time for very large datasets.
    """
    path = Path(bed_file)
    assert path.exists(), f"path {path} does not exist"
    out = sys.stdout if out == "-" else Path(out).open("w")

    # Read entire file.
    with path.open() as f:
        lines = [line.split() for line in f.read().splitlines()]

    # Calculate breaks.
    scores = sorted(-float(line[11]) for line in lines)
    breaks = jenkspy.jenks_breaks(scores, n_classes=2)

    # Assign labels and output.
    for line in lines:
        line.append("pos" if -float(line[11]) < breaks[1] else "neg")
        out.write(" ".join(line))
        out.write("\n")


if __name__ == "__main__":
    cli()
