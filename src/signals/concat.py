"""Provides various signal-related functionality."""

import logging
from collections import deque
from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np

from signals.fast5 import Fast5

log = logging.getLogger("signals.signal")


@dataclass
class Signal:
    """A concatenated signal at a certain position."""
    position: int
    coverage: int
    data: np.ndarray


def concat(fasts: list[Fast5]) -> Iterator[Signal]:
    """Concat signals from a sorted iterable of FAST5 files and yield them."""
    fasts.sort(key=lambda f: (f.start, f.end))
    start = 0 # Up to which position has been processed already?

    for i, f in enumerate(fasts):
        # Find minimum overlap of current file with next files.
        for j in range(i + 1, len(fasts)):
            if fasts[j].start >= f.end:
                break

        # Process if at least one file within overlap.
        overlap = fasts[i:j]
        if overlap:
            start = max(start, f.start)
            if start < f.end:
                yield from _concat_overlap(overlap, start, f.end)
                start = max(start, f.end - 1)

        # Get rid of the processed range.
        fasts[i] = None


def _concat_overlap(fasts: list[Fast5], start: int, end: int) -> Iterator[Signal]:
    """Yield concatenated signals from a FAST5 subset, but only for the given range."""
    for i in range(start, end):
        signal = Signal(position=i, coverage=0, data=[])
        for f in fasts:
            if f.start <= i < f.end:
                data, length = f.signal
                length = length[i - f.start]
                signal.data.extend(data[f.position:f.position + length])
                signal.coverage += 1
                if 1 < i < end - 2:
                    f.position += length
        signal.data = np.array(signal.data)
        yield signal
