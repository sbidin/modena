"""All tests."""

from collections.abc import Iterator
from io import StringIO
from pathlib import Path

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


@pytest.mark.parametrize(
    ("acid", "coverage", "resample"),
    _param_combos())
def test_output(acid: str, coverage: int, resample: int) -> None:
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
        out="-",
        random_seed=42, # Needed to make the test deterministic.
        resample_size=resample,
        strand="+", # Hardcoded to + because CI vs local can be non-deterministic.
        to_position=None)
    config.out = StringIO() # Override, since we're writing to memory.

    compare_datasets(config)

    # Read expected output.
    with open(f"tests/out/{acid}.coverage-{coverage}.resample-{resample}.bed") as f:
        exp = f.read()

    # Compare expected with actual output.
    for a, b in zip(config.out.getvalue().splitlines(), exp.splitlines(), strict=False):
        assert a == b, "output mismatch"
