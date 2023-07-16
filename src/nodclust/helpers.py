"""Provides various helpers."""

import contextlib
import logging
from timeit import default_timer

from modena.config import Config
from modena.fast5 import PATTERN_MATCH_POSITIONS

log = logging.getLogger("modena")


def position_within_bounds(
        pos: int,
        config: Config,
        with_window: bool=False) \
        -> bool:
    """Return whether a position satisfies various filters."""
    # Check for pattern position match.
    if config.pattern:
        left, right = pos, pos + 1
        if with_window:
            left -= config.WINDOW_SIZE // 2
            right += config.WINDOW_SIZE // 2
        if not any(pos in PATTERN_MATCH_POSITIONS for pos in range(left, right)):
            return False

    pos += 1 # Perform from/to bounds checks in 1-based indexing.

    # Check left bound.
    if config.from_position is not None:
        left = config.from_position
        if with_window:
            left -= config.WINDOW_SIZE // 2
        if pos < left:
            return False

    # Check right bound.
    if config.to_position is not None:
        right = config.to_position
        if with_window:
            right += config.WINDOW_SIZE // 2
        if pos > right:
            return False

    return True


@contextlib.contextmanager
def timer(message: str) -> None:
    """Output how long a block of code took to execute."""
    time_from = default_timer()
    yield
    seconds = default_timer() - time_from
    log.debug(f"{message} took {seconds:.3f} seconds")
