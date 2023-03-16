import sys
from pathlib import Path

import click

from signals.fast5 import Fast5


@click.command()
@click.argument("dataset1")
@click.argument("dataset2")
def _main(dataset1: str, dataset2: str) -> None:
    dataset1, dataset2 = Path(dataset1), Path(dataset2)
    assert dataset1.exists(), f"no such path: {dataset1}"
    assert dataset2.exists(), f"no such path: {dataset2}"
    reads1 = [(f.position_range(), f) for f in Fast5.from_path(dataset1)]
    reads2 = [(f.position_range(), f) for f in Fast5.from_path(dataset2)]
    reads1.sort()
    reads2.sort()
    print(reads1)
    print(reads2)


if __name__ == "__main__":
    _main()
