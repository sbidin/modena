"""All tests."""

from io import StringIO
from pathlib import Path

import pytest
from signals.run import run_on_datasets


@pytest.mark.parametrize("out_format", ("csv", "bedgraph"))
def test_output_equals_old_script(out_format: str) -> None:
    """CSV output equals output of old script."""
    out = StringIO()
    inp1, inp2 = Path("tests/inp/1"), Path("tests/inp/2")
    run_on_datasets(inp1, inp2, 1, out, out_format)
    out = out.getvalue()
    with open(f"tests/out.{out_format}") as f:
        exp = f.read()
    for a, b in zip(out.splitlines(), exp.splitlines(), strict=False):
        assert a == b, "output mismatch"
