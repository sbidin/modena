"""Provides various FAST5 parsing functionality."""

from __future__ import annotations

import logging
import re
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
    """A subset of a parsed FAST5 file.

    This is fetched prior to fetching any signal data. Position data is used to
    later on pick which files we actually want to process and in which order.
    """
    path: Path
    is_rna: bool
    strand: str
    chrom: str
    start: int
    end: int

    @staticmethod
    def from_path(path: Path, strand: str, chrom: str) -> Iterator[Fast5]:
        """Yield all FAST5 files found under the given path.

        Optionally allows filtering by strand or chromosome.
        """
        for f in Fast5._from_path(path):
            satisfies = True
            if strand is not None:
                satisfies = satisfies and strand == f.strand
            if chrom is not None:
                satisfies = satisfies and re.search(chrom, f.chrom) is not None
            if satisfies:
                yield f

    @staticmethod
    def _from_path(path: Path) -> Iterator[Fast5]:
        """Yield all FAST5 files found under the given path."""
        if path.is_file():
            yield from Fast5._maybe_from_file_path(path)
            return
        for sub in path.glob("**/*.*"):
            if sub.name.lower().endswith(".fast5") and sub.is_file():
                yield from Fast5._maybe_from_file_path(sub)

    def overlap(self, fasts: list[Fast5]) -> Fast5 | None:
        """Return the first file whose positions overlap with self."""
        j = bisect_left(fasts, (self.start, -self.end), key=lambda f: (f.start, -f.end))
        for i in range(max(j - 1, 0), len(fasts)):
            f = fasts[i]
            if f.start >= self.end:
                break
            if self.start < f.end and f.start < self.end:
                return f

    @cached_property
    def positions(self) -> list[int]:
        """Returns positions within the underlying array to extract signals from."""
        _, lengths = self.signal
        ps = [0]
        for i in range(lengths.shape[0]):
            ps.append(ps[-1] + lengths[i])
        return ps

    def signal_at(self, i: int, resample: int | None) -> np.ndarray:
        """Return (optionally resampled) signal data at a certain position in the file."""
        signal, length = self.signal
        assert signal is not None, f"file {self.path} missing signal data"
        i -= self.start
        signal = signal[self.positions[i]:self.positions[i] + length[i]]

        # Optionally performs resampling.
        if resample > 0 and resample != len(signal):
            signal = np.random.choice(signal, size=resample, replace=True)

        return signal

    @staticmethod
    def _maybe_from_file_path(path: Path) -> Iterator[Fast5]:
        """Attempt parsing the basics of a FAST5 file."""
        try:
            with h5py.File(path) as f:
                yield from Fast5._parse_file(f, path)
        except AssertionError as err:
            log.debug(f"skipped file {path} because {err}")

    @staticmethod
    def _parse_file(f: h5py.File, path: Path) -> Iterator[Fast5]:
        """Get all basic file contents that are useful to us."""
        tpl = f.get("Analyses/RawGenomeCorrected_000/BaseCalled_template")
        assert tpl, "missing BaseCalled_template"
        is_rna = (tpl.attrs.get("rna"))
        assert isinstance(is_rna, np.bool_), "missing BaseCalled_template.rna"
        aln = tpl.get("Alignment")
        assert aln, "missing Alignment"
        strand = aln.attrs.get("mapped_strand")
        assert strand, "missing Alignment.mapped_strand"
        chrom = aln.attrs.get("mapped_chrom")
        assert chrom, "missing Alignment.mapped_chrom"
        start = aln.attrs.get("mapped_start")
        assert isinstance(start, np.int64), "missing Alignment.mapped_start"
        evt = tpl.get("Events")
        assert evt, "missing Events"
        size = evt.shape[0]
        assert size > 0, "empty Events"
        yield Fast5(
            path=path,
            is_rna=is_rna,
            strand=strand,
            chrom=chrom,
            start=start,
            end=start + size)

    @cached_property
    def signal(self) -> tuple[np.ndarray, np.ndarray]:
        """Get signal and signal length data."""
        with h5py.File(self.path) as f:
            tpl = f.get("Analyses/RawGenomeCorrected_000/BaseCalled_template")
            evt = tpl.get("Events")
            read_start_rel_to_raw = evt.attrs.get("read_start_rel_to_raw")
            shift = tpl.attrs.get("shift")
            scale = tpl.attrs.get("scale")
            rds = f.get("Raw/Reads")
            rdn = list(rds.keys())
            sgn = rds[rdn[0]]["Signal"]
            if self.is_rna:
                sgn = np.flip(sgn)
            signal = (sgn[read_start_rel_to_raw:] - shift) / scale
            lengths = evt["length"]
            return signal, lengths

    def __repr__(self) -> str:
        """Get a simple Fast5 representation for debugging."""
        return f"Fast5({self.path}, {self.start}-{self.end}, {'rna' if self.is_rna else 'dna'}, {self.strand}, {self.chrom})"
