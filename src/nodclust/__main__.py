#!/usr/bin/env python

"""Runs the main application."""

import logging
import os
import sys
from pathlib import Path

import click
import jenkspy

from nodclust.run import run_on_datasets

try:
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
except ValueError:
    logging.basicConfig(level="INFO")
finally:
    log = logging.getLogger("nodclust")


@click.group()
def cli() -> None:
    pass


@click.command()
@click.argument("dataset1")
@click.argument("dataset2")
@click.option("-a", "--acid", type=str, default="autodetect", help="Filter by acid, dna or rna")
@click.option("--force-acid", is_flag=True, type=bool, default=False, help="Force read files as specified by --acid")
@click.option("-s", "--strand", type=str, default=None, help="Filter by strand, '+' or '-'")
@click.option("-c", "--chromosome", type=str, default=None, help="Filter by chromosome regex")
@click.option("-f", "--from-position", type=int, default=None, help="Filter by minimum position (inclusive)")
@click.option("-t", "--to-position", type=int, default=None, help="Filter by maximum position (inclusive)")
@click.option("-m", "--min-coverage", type=int, default=5, help="Filter by minimum coverage (default 5)")
@click.option("-r", "--resample-size", type=int, default=10, help="Signal resample size; 0 to disable (default 10)")
@click.option("-o", "--out", default="-", help="Output to a given path")
@click.option("--no-distance-sum", is_flag=True, type=bool, default=False, help="Don't sum neighbour position distances")
@click.option("--random-seed", type=int, default=None, help="Force a random seed, for reproducibility")
def compare(
        dataset1: str,
        dataset2: str,
        acid: str,
        force_acid: bool,
        strand: str | None,
        chromosome: str | None,
        from_position: int | None,
        to_position: int | None,
        min_coverage: int,
        resample_size: int,
        out: str,
        no_distance_sum: bool | None,
        random_seed: int | None) \
        -> None:
    """Compare two datasets & output an annotated bedMethyl file."""
    xs_path, ys_path = Path(dataset1), Path(dataset2)
    assert xs_path.exists(), f"no such path exists: {xs_path}"
    assert ys_path.exists(), f"no such path exists: {ys_path}"

    # Perform various input argument sanity checks.
    assert acid in ("dna", "rna", "autodetect"), "--acid must be one of dna, rna or autodetect"
    assert strand in (None, "+", "-"), "--strand must be one of '+' or '-'"
    if force_acid:
        assert acid != "autodetect", "cannot --force-acid without specifying --acid"
    if from_position is not None:
        assert from_position >= 0, "--from-position must be at least zero"
    if to_position is not None:
        assert to_position >= 0, "--from-position must be at least zero"

    out = sys.stdout if out == "-" else Path(out).open("w")
    run_on_datasets(
        xs_path,
        ys_path,
        acid,
        strand,
        chromosome,
        force_acid,
        from_position,
        to_position,
        min_coverage,
        resample_size if resample_size > 0 else 0,
        not no_distance_sum,
        out,
        random_seed)


@click.command()
@click.argument("bed_file")
@click.option("-o", "--out", default="-", help="Output to a given path (default stdout)")
def jenks(bed_file: str, out: str) -> None:
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
    cli.add_command(compare)
    cli.add_command(jenks)
    cli()
