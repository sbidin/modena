"""All tests."""

from io import StringIO
from pathlib import Path

import pytest
from signals.run import run_on_datasets


@pytest.mark.parametrize(
    ("suffix", "resample", "distsum"),
    (
        ("default", 10, True),
        ("no-distsum", 10, False),
        ("no-resample", 0, True),
    )
)
def test_default(suffix: str, resample: int, distsum: bool) -> None:
    """CSV output equals output of old script."""
    out = StringIO()
    inp1, inp2 = Path("tests/inp/dna/1"), Path("tests/inp/dna/2")
    run_on_datasets(
        inp1,
        inp2,
        "autodetect",
        None,
        None,
        False,
        1,
        resample,
        distsum,
        out,
        random_seed=42)
    out = out.getvalue()
    with open(f"tests/out/dna/out.{suffix}.bed") as f:
        exp = f.read()
    for a, b in zip(out.splitlines(), exp.splitlines(), strict=False):
        assert a == b, "output mismatch"
