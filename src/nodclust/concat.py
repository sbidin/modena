"""Provides an iterator concatenating signals from multiple FAST5 files."""

import logging
from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np

from nodclust.fast5 import Fast5

log = logging.getLogger("nodclust")


@dataclass
class Signal:
    """A concatenated signal at a certain position."""
    position: int
    coverage: int
    data: np.ndarray


def concat_pairs(
        xs: list[Fast5],
        ys: list[Fast5],
        min_coverage: int,
        resample: int,
        from_position: int | None,
        to_position: int | None) \
        -> Iterator[tuple[Signal, Signal]]:
    """Yield tuples of same-position concatenated signals."""
    xs = _concat(xs, min_coverage, resample, from_position, to_position)
    ys = _concat(ys, min_coverage, resample, from_position, to_position)
    x, y = next(xs, None), next(ys, None)
    while True:
        if x is None or y is None:
            break
        elif x.position > y.position:
            y = next(ys, None)
        elif y.position > x.position:
            x = next(xs, None)
        else:
            yield x, y
            x, y = next(xs, None), next(ys, None)


def _concat(
        fasts: list[Fast5],
        min_coverage: int,
        resample: int,
        from_position: int | None,
        to_position: int | None) \
        -> Iterator[Signal]:
    """Yield concatenated signals from FAST5 files."""
    fasts.sort(key=lambda f: (f.start, f.end))
    start = 0 # Up to which position has been processed already?

    for i, f in enumerate(fasts):
        # Find minimum overlap of current file with next files.
        j = i + 1
        while j < len(fasts) and fasts[j].start < f.end:
            j += 1

        # Process if at least one file within overlap.
        overlap = fasts[i:j]
        start = max(start, f.start)
        if start < f.end and j - i >= min_coverage:
            yield from _concat_range(
                overlap,
                start,
                f.end,
                min_coverage,
                resample,
                from_position,
                to_position)
        start = max(start, f.end)

        # Get rid of the processed range.
        fasts[i] = None


def _concat_range(
        fasts: list[Fast5],
        start: int,
        end: int,
        min_coverage: int,
        resample: int,
        from_position: int | None,
        to_position: int | None) \
        -> Iterator[Signal]:
    """Yield concatenated signals for a range of positions."""
    # For each position...
    for i in range(start, end):

        # Skip positions outside the desired range. However, still include at
        # least window size positions around the desired range, since this is
        # important for purposes of distance summing.
        if from_position is not None and from_position - 5 > i + 1 or \
           to_position is not None and i + 1 > to_position + 5:
            continue

        signal = Signal(position=i, coverage=0, data=[])

        # For each file containing this position...
        for f in (f for f in fasts if f.start <= i < f.end):
            signal.data.extend(f.signal_at(i, resample))
            signal.coverage += 1

        if signal.coverage >= min_coverage:
            signal.data = np.array(signal.data)
            yield signal
