"""Provides an iterator concatenating signals from multiple FAST5 files."""

import logging
from dataclasses import dataclass
from typing import Iterator, List, Tuple

import numpy as np

from modena.config import Config
from modena.fast5 import Fast5
from modena.helpers import position_within_bounds

log = logging.getLogger("modena")


@dataclass
class Signal:
    """A concatenated signal at a certain position."""
    position: int
    coverage: int
    data: np.ndarray


def concat_pairs(
        xs: List[Fast5],
        ys: List[Fast5],
        config: Config) \
        -> Iterator[Tuple[Signal, Signal]]:
    """Yield tuples of same-position concatenated signals."""
    xs, ys = _concat(xs, config), _concat(ys, config)
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


def _concat(fasts: List[Fast5], config: Config) -> Iterator[Signal]:
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
        if start < f.end and j - i >= config.min_coverage:
            yield from _concat_range(overlap, start, f.end, config)
        start = max(start, f.end)

        # Get rid of the processed range.
        fasts[i] = None


def _concat_range(
        fasts: List[Fast5],
        start: int,
        end: int,
        config: Config) \
        -> Iterator[Signal]:
    """Yield concatenated signals for a range of positions."""
    # For each position...
    for pos in range(start, end):

        # Skip positions outside of bounds.
        if not position_within_bounds(pos, config, with_window=True):
            continue

        signal = Signal(position=pos, coverage=0, data=[])

        # For each file containing this position...
        for f in (f for f in fasts if f.start <= pos < f.end):
            signal.data.extend(f.signal_at(pos, config.resample_size))
            signal.coverage += 1

        if signal.coverage >= config.min_coverage:
            signal.data = np.array(signal.data)
            yield signal
