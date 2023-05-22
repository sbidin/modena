"""All tests."""

from io import StringIO
from pathlib import Path

from signals.run import run_on_datasets


def test_output_equals_old_script() -> None:
    """CSV output equals output of old script."""
    out = StringIO()
    inp1, inp2 = Path("tests/inp/1"), Path("tests/inp/2")
    run_on_datasets(inp1, inp2, 1, None, False, out)
    out = out.getvalue()
    with open("tests/out/out.bed") as f:
        exp = f.read()
    for a, b in zip(out.splitlines(), exp.splitlines(), strict=False):
        assert a == b, "output mismatch"
