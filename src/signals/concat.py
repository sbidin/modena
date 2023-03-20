"""Provides an iterator concatenating signals from multiple FAST5 files."""

from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np

from signals.fast5 import Fast5


@dataclass
class Signal:
    """A concatenated signal at a certain position."""
    position: int
    coverage: int
    data: np.ndarray


def concat_pairs(
        xs: list[Fast5],
        ys: list[Fast5],
        min_coverage: int) \
        -> Iterator[tuple[Signal, Signal]]:
    """Yield tuples of same-position concatenated signals."""
    xs = concat(xs, min_coverage)
    ys = concat(ys, min_coverage)
    x = next(xs, None)
    y = next(ys, None)

    while True:
        if x is None or y is None:
            break

        if x.position > y.position:
            y = next(ys, None)
            continue

        if y.position > x.position:
            x = next(xs, None)
            continue

        assert x.position == y.position
        yield x, y
        x = next(xs, None)
        y = next(ys, None)


def concat(fasts: list[Fast5], min_coverage: int) -> Iterator[Signal]:
    """Yield concatenated signals from FAST5 files."""
    fasts.sort(key=lambda f: (f.start, f.end))
    start = 0 # Up to which position has been processed already?

    for i, f in enumerate(fasts):
        # Find minimum overlap of current file with next files.
        j = i + 1
        for j in range(i + 1, len(fasts)):
            if fasts[j].start >= f.end:
                break

        # Process if at least one file within overlap.
        overlap = fasts[i:j]
        if overlap:
            start = max(start, f.start)
            if start < f.end:
                yield from _concat_range(overlap, start, f.end, min_coverage)
                start = max(start, f.end - 1)

        # Get rid of the processed range.
        fasts[i] = None


def _concat_range(
        fasts: list[Fast5],
        start: int,
        end: int,
        min_coverage: int) \
        -> Iterator[Signal]:
    """Yield concatenated signals for a range of positions."""
    # For each position...
    for i in range(start, end):
        signal = Signal(position=i, coverage=0, data=[])

        # For each file containing this position...
        for f in (f for f in fasts if f.start <= i < f.end):
            signal.data.extend(f.signal_at(i))
            signal.coverage += 1

        if signal.coverage >= min_coverage:
            signal.data = np.array(signal.data)
            yield signal
