"""All tests."""

from io import StringIO
from pathlib import Path

from signals.run import run_on_datasets


def test_csv_equals_old_script() -> None:
    """CSV output equals output of old script."""
    out = StringIO()
    inp1, inp2 = Path("tests/inp/1"), Path("tests/inp/2")
    run_on_datasets(inp1, inp2, 1, out)
    out = out.getvalue()
    with open("tests/out.csv") as f:
        exp = f.read()
        assert out == exp, "output mismatch"