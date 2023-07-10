"""All tests."""

import os
from pathlib import Path
from typing import Iterator

import pytest
from nodclust.config import Config
from nodclust.run import compare_datasets


def _param_combos() -> Iterator[tuple]:
    """Get various combinations of test parameters."""
    for path in Path("tests/out").glob("*.bed"):
        name = path.name.split(".")
        yield (
            name[0], # RNA or DNA.
            int(name[1].split("-")[1]), # Minimum coverage.
            int(name[2].split("-")[1])) # Resample size.


@pytest.fixture
def _cleanup() -> None:
    """Remove the resulting "out.test.bed" file, if present."""
    yield
    out = Path("out.test.bed")
    if out.exists():
        os.remove(out)


@pytest.mark.parametrize(
    ("acid", "coverage", "resample"),
    _param_combos())
def test_output(acid: str, coverage: int, resample: int, _cleanup) -> None:
    """Test that application output matches one on disk."""
    config = Config.from_cmd_args(
        dataset1=Path(f"tests/inp/{acid}-1"),
        dataset2=Path(f"tests/inp/{acid}-2"),
        acid=acid,
        chromosome=None,
        force_acid=True,
        from_position=None,
        min_coverage=coverage,
        no_distance_sum=False,
        output_bed="out.test.bed",
        pattern=None,
        random_seed=42, # Needed to make the test deterministic.
        resample_size=resample,
        strand="+", # Hardcoded to + because CI vs local can be non-deterministic.
        to_position=None)

    # This will write results into "out.test.bed".
    compare_datasets(config)

    # Compare expected with actual output.
    with open(f"tests/out/{acid}.coverage-{coverage}.resample-{resample}.bed") as out:
        with open("out.test.bed") as exp:
            for out_line, exp_line in zip(out, exp):
                assert out_line == exp_line, "output vs expected line mismatch"
