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

log = logging.getLogger("nodclust")


@dataclass
class Fast5:
    """A subset of a parsed FAST5 file.

    This is fetched prior to fetching any signal data. Position data is used to
    later on pick which files we actually want to process and in which order.
    """
    path: Path
    type: str
    strand: str
    chrom: str
    start: int
    end: int

    @staticmethod
    def from_path(
            path: Path,
            type: str,
            strand: str,
            chrom: str,
            force_type: bool) \
            -> Iterator[Fast5]:
        """Yield all FAST5 files found under the given path.

        Optionally allows filtering by type, strand or chromosome.
        """
        # Only select one chromosome and strand -- skip others.
        type_selected = None
        chrom_selected = None
        strand_selected = None

        log.debug(f"filtering files in {path}")
        for f in Fast5._from_path(path, forced_type=type if force_type else None):
            satisfies = True
            if type != "autodetect":
                satisfies = satisfies and type == f.type
            if strand is not None:
                satisfies = satisfies and strand == f.strand
            if chrom is not None:
                satisfies = satisfies and re.search(chrom, f.chrom) is not None
            if not satisfies:
                continue

            # Don't compare RNA files to DNA files and vice-versa.
            if type_selected is None:
                type_selected = f.type
                log.debug(f"selected type {type_selected}")
            elif type_selected != f.type:
                log.debug(f"skipped {f.path} due to incompatible type {f.type}")
                continue

            # If we encounter multiple chromosomes, there's a mixed dataset. We
            # currently support only one chromosome at a time so it's up to the
            # user to filter by a specific one.
            if chrom_selected is None:
                chrom_selected = f.chrom
                log.debug(f"selected chromosome {chrom_selected}")
            elif chrom_selected != f.chrom:
                log.debug(f"skipped {f.path} due to incompatible chromosome {f.chrom}")
                continue

            # Similarly for strand -- don't compare across different strands.
            if strand_selected is None:
                strand_selected = f.strand
                log.debug(f"selected strand {strand_selected}")
            elif strand_selected != f.strand:
                log.debug(f"skipped {f.path} due to incompatible strand {f.strand}")
                continue

            log.debug(f"selected {f.path}")
            yield f

    @staticmethod
    def _from_path(path: Path, forced_type: str | None) -> Iterator[Fast5]:
        """Yield all FAST5 files found under the given path."""
        if path.is_file():
            yield from Fast5._maybe_from_file_path(path, forced_type)
            return
        for sub in path.glob("**/*.*"):
            if sub.name.lower().endswith(".fast5") and sub.is_file():
                yield from Fast5._maybe_from_file_path(sub, forced_type)

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

    def signal_at(self, i: int, resample: int) -> np.ndarray:
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
    def _maybe_from_file_path(path: Path, forced_type: str | None) -> Iterator[Fast5]:
        """Attempt parsing the basics of a FAST5 file."""
        try:
            with h5py.File(path) as f:
                yield from Fast5._parse_file(f, path, forced_type)
        except AssertionError as err:
            log.debug(f"skipped {path} due to {err}")

    @staticmethod
    def _parse_file(f: h5py.File, path: Path, forced_type: str | None) -> Iterator[Fast5]:
        """Get all basic file contents that are useful to us."""
        tpl = f.get("Analyses/RawGenomeCorrected_000/BaseCalled_template")
        assert tpl, "missing BaseCalled_template"

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

        # Try to detect whether the file is DNA or RNA by checking a few attributes.
        is_rna = tpl.attrs.get("rna")
        if is_rna is None:
            context_tags = f.get("UniqueGlobalKey/context_tags")
            if context_tags:
                experiment_type = context_tags.attrs.get("experiment_type")
                if experiment_type is not None:
                    is_rna = np.bool_(experiment_type == "rna")

        # If we couldn't autodetect, but are forcing a type, then force the type.
        if not isinstance(is_rna, np.bool_) and forced_type is not None:
            is_rna = np.bool_(forced_type == "rna")

        assert isinstance(is_rna, np.bool_), "missing BaseCalled_template.rna or UniqueGlobalKey.context_tags.experiment_type"

        yield Fast5(
            path=path,
            type="rna" if is_rna else "dna",
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
            sgn = rds[rdn[0]]["Signal"][()]
            if self.type == "rna":
                sgn = np.flip(sgn)
            signal = (sgn[read_start_rel_to_raw:] - shift) / scale
            lengths = evt["length"]
            return signal, lengths

    def __repr__(self) -> str:
        """Get a simple Fast5 representation for debugging."""
        return f"{self.type.upper()}[{self.chrom}{self.strand}]({self.path.name}, {self.start}-{self.end})"
