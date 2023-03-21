"""Provides various helper functionality."""

import subprocess
from pathlib import Path


def size_at_path(path: Path) -> int:
    """Get the size in bytes of a file or directory at `path`."""
    try:
        size = subprocess.check_output(["du", "-sb", str(path)])
        size = int(size.split()[0].decode("utf-8"))
        assert size > 0, "size is zero"
        return size
    except:
        return 0


def order_paths_by_size(a: Path, b: Path) -> tuple[Path, Path, bool]:
    """Order given paths by size in ascending order.

    Additionally return whether the paths have been reordered.
    """
    datasets = [a, b]
    datasets.sort(key=size_at_path)
    return *datasets, datasets[0] != a

