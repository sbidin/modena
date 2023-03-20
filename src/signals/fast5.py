"""Provides various FAST5 parsing functionality."""

from __future__ import annotations

import logging
from bisect import bisect_left
from collections.abc import Iterator
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import h5py
import numpy as np

log = logging.getLogger("signals.fast5")


@dataclass
class Fast5:
    """A tiny subset of a parsed FAST5 file.

    This is fetched prior to fetching any signal data. Position data is used to
    later on pick which files we actually want to process and in which order.
    """
    path: Path
    start: int
    end: int

    position: int = 0 # Keeps track of our position within the file while reading.
    last_i: int = None # Last position accessed.

    @staticmethod
    def from_path(path: Path) -> Iterator[Fast5]:
        """Yield all FAST5 files found under the given path."""
        if path.is_file():
            yield from Fast5._maybe_from_file_path(path)
            return
        for sub in path.glob("**/*.*"):
            if sub.name.lower().endswith(".fast5") and sub.is_file():
                yield from Fast5._maybe_from_file_path(sub)

    @staticmethod
    def _maybe_from_file_path(path: Path) -> Iterator[Fast5]:
        """Attempt parsing the basics of a FAST5 file."""
        try:
            log.debug(f"parsing attributes from {path}")
            with h5py.File(path) as f:
                yield from Fast5._parse_file(f, path)
        except AssertionError as err:
            log.debug(f"skipped file {path} because {err}")

    @staticmethod
    def _parse_file(f: h5py.File, path: Path) -> Iterator[Fast5]:
        """Get all basic file contents that are useful to us."""
        tpl = f.get("Analyses/RawGenomeCorrected_000/BaseCalled_template")
        assert tpl, "missing BaseCalled_template"
        aln = tpl.get("Alignment")
        assert aln, "missing Alignment"
        strand = aln.attrs.get("mapped_strand")
        assert strand == "+", "strand not positive"
        evt = tpl.get("Events")
        assert evt, "missing Events"
        start = aln.attrs.get("mapped_start")
        assert isinstance(start, np.int64), "missing mapped_start"
        size = evt.shape[0]
        assert size > 0, "empty Events"
        fast5 = Fast5(path=path, start=start, end=start + size, last_i=start)
        log.debug(f"parsed file {path} {fast5}")
        yield fast5

    def signal_at(self, i: int) -> list[float] | None:
        """Return signal data from a certain position in the file."""
        signal, length = self.signal
        assert signal is not None, f"file {self.path} missing signal data"

        # Convert absolute positions to within-signal ones.
        i -= self.start
        size = self.end - self.start

        # We might have skipped a portion of the positions. Time to catch up.
        if i > 2 and self.last_i < i:
            self.position += sum(length[self.last_i + 1:min(i, size - 2)])

        signal = signal[self.position:self.position + length[i]]
        if 1 < i < size - 2: # Beginnings and ends are special.
            self.position += length[i]

        self.last_i = i
        return signal

    @cached_property
    def signal(self) -> tuple[np.ndarray, np.ndarray] | None:
        """Get signal and base data."""
        try:
            with h5py.File(self.path) as f:
                # Get useful groups.
                tpl = f.get("Analyses/RawGenomeCorrected_000/BaseCalled_template")
                assert tpl, "missing BaseCalled_template"
                aln = tpl.get("Alignment")
                assert aln, "missing Alignment"
                evt = tpl.get("Events")
                assert evt, "missing Events"

                # Get signal modifiers.
                read_start_rel_to_raw = evt.attrs.get("read_start_rel_to_raw")
                assert isinstance(read_start_rel_to_raw, np.int64), "missing read_start_rel_to_raw"
                shift = tpl.attrs.get("shift")
                assert isinstance(shift, np.float64),  "missing shift"
                scale = tpl.attrs.get("scale")
                assert isinstance(scale, np.float64), "missing scale"

                # Get the signal and base data.
                rds = f.get("Raw/Reads")
                assert rds, "missing Reads"
                rdn = list(rds.keys())
                assert len(rdn) == 1, "not exactly one Read_XYZ found"
                sgn = rds[rdn[0]]["Signal"]
                assert sgn, "missing Signal"
                signal = (sgn[read_start_rel_to_raw:] - shift) / scale
                bases = evt["base"]
                assert isinstance(bases, np.ndarray), "missing base"
                lengths = evt["length"]
                assert isinstance(lengths, np.ndarray), "missing length"
                return signal, lengths
        except AssertionError as err:
            log.debug(f"skipped file {self.path} because {err}")
            return None, None

    def __repr__(self) -> str:
        """Get a simple Fast5 representation for debugging."""
        return f"F({self.start}-{self.end})"

    def overlap(self, fasts: list[Fast5]) -> Fast5 | None:
        """Return the first file whose positions overlap with self."""
        j = bisect_left(fasts, (self.start, -self.end), key=lambda f: (f.start, -f.end))
        for i in range(max(j - 1, 0), len(fasts)):
            f = fasts[i]
            if f.start >= self.end:
                break
            if self.start < f.end and f.start < self.end:
                return f
