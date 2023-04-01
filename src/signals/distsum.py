"""Provides a distance sum generator."""

from collections import deque
from collections.abc import Iterator


def distsum(
        stats: Iterator[tuple[int, float]],
        window_size: int) \
        -> Iterator[tuple[int, float]]:
    """Assign each position a distance that is the sum of its window."""
    assert window_size > 1
    assert window_size % 2 == 1
    window = deque(maxlen=window_size)

    # Starting out, populate the entire window.
    for _ in range(window_size):
        stat = next(stats, None)
        if stat is not None:
            window.append(stat)

    # Handle the prefix first as their left neighbours are incomplete.
    for i in range(window_size // 2):
        yield from _distsum_at(window, i)

    # As the rest are fetched, strip out the prefix.
    for stat in stats:
        yield from _distsum_at(window, window_size // 2)
        window.popleft()
        window.append(stat)

    # Handle the suffix, with incomplete right neighbours.
    for i in range(window_size // 2, window_size):
        yield from _distsum_at(window, i)






def _distsum_at(
        window: deque[tuple[int, float]],
        i: int) \
        -> Iterator[tuple[int, float]]:
    assert 0 <= i < len(window)
    span = window.maxlen // 2
    sum_dists = 0
    for j in range(i - span, i + span + 1):
        if 0 <= j < len(window) and abs(window[i][0] - window[j][0]) <= span:
            sum_dists += window[j][1]
    yield window[i][0], sum_dists
