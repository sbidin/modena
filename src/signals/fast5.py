"""Provides various FAST5 parsing functionality."""

import logging
import os
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
log = logging.getLogger("signals.fast5")


@dataclass
class Fast5:
    """A parsed FAST5 file."""
    path: Path
    mapped_start: int
    signal: np.ndarray | None
    bases: np.ndarray

    @staticmethod
    def from_path(path: Path) -> Iterator["Fast5"]:
        """Yield all FAST5 files found under the given path."""
        if path.is_file():
            yield from Fast5._maybe_from_file_path(path)
            return
        for sub in path.glob("**/*.*"):
            if sub.name.lower().endswith(".fast5") and sub.is_file():
                yield from Fast5._maybe_from_file_path(sub)

    def position_range(self) -> tuple[int, int]:
        """Get the range of positions found within the file."""
        return self.mapped_start, self.mapped_start + self.bases.shape[0]

    @staticmethod
    def _maybe_from_file_path(path: Path) -> Iterator["Fast5"]:
        """Attempt parsing a FAST5 file."""
        try:
            log.debug(f"parsing attributes from {path}")
            with h5py.File(path) as f:
                yield from Fast5._parse_file(f, path)
        except AssertionError as err:
            log.warning(f"skipped file {path} because {err}")
    
    @staticmethod
    def _parse_file(f: h5py.File, path: Path) -> Iterator["Fast5"]:
        """Get all file contents that are useful to us."""
        # Get various useful groups.
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

        # Get the signal.
        rds = f.get("Raw/Reads")
        assert rds, "missing Reads"
        rdn = list(rds.keys())
        assert len(rdn) == 1, "not exactly one Read_XYZ found"
        sgn = rds[rdn[0]]["Signal"]
        assert sgn, "missing Signal"
        signal = (sgn[read_start_rel_to_raw:] - shift) / scale

        # Get other useful attributes.
        mapped_start = aln.attrs.get("mapped_start")
        assert isinstance(mapped_start, np.int64), "missing mapped_start"
        bases = evt["base"]
        assert isinstance(bases, np.ndarray), "missing base"
    
        yield Fast5(
            path=path,
            mapped_start=mapped_start,
            signal=signal,
            bases=bases)
