"""Provides various helpers."""

import contextlib
import logging
from timeit import default_timer

from nodclust.config import Config

log = logging.getLogger("nodclust")


def position_within_bounds(
        pos: int,
        config: Config,
        with_window: bool=False) \
        -> bool:
    """Return whether a position satisfies from/to bounds, inclusively."""
    pos += 1 # Comparisons are done in 1-based indexing.

    # Check left bound.
    if config.from_position is not None:
        left = config.from_position
        if with_window:
            left -= config.WINDOW_SIZE
        if pos < left:
            return False

    # Check right bound.
    if config.to_position is not None:
        right = config.to_position
        if with_window:
            right += config.WINDOW_SIZE
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
