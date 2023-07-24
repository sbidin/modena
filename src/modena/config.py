"""Provides a configuration object."""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np


@dataclass
class Config:
    """A bundle of all application configuration values."""
    acid: str
    chromosome: Optional[re.Pattern]
    dataset1: Path
    dataset2: Path
    force_acid: bool
    from_position: Optional[int]
    min_coverage: int
    no_distance_sum: bool
    output_bed: Path
    pattern: Optional[re.Pattern]
    random_seed: Optional[int]
    resample_size: int
    strand: Optional[str]
    to_position: Optional[int]
    WINDOW_SIZE: int = 5

    @staticmethod
    def from_cmd_args(**kwargs: Dict[str, Any]) -> "Config":
        """Create a `Config` given `click` commandline arguments."""
        config = Config(**kwargs)
        config._postprocess()
        config._assert()
        return config

    def _postprocess(self) -> None:
        """Apply various transformations to configuration values."""
        try:
            if self.chromosome is not None:
                self.chromosome = re.compile(self.chromosome)
        except re.error as err:
            raise AssertionError("chromosome regex must be valid") from err

        self.dataset1 = Path(self.dataset1)
        self.dataset2 = Path(self.dataset2)
        self.output_bed = Path(self.output_bed)

        try:
            if self.pattern is not None:
                self.pattern = re.compile(self.pattern, re.IGNORECASE)
        except re.error as err:
            raise AssertionError("pattern regex must be valid") from err

        if self.random_seed is not None:
            np.random.seed(self.random_seed)

    def _assert(self) -> None:
        """Assert that configuration is sane."""
        assert self.acid in ("dna", "rna", "autodetect"), "acid must be one of: 'dna', 'rna'"
        assert self.dataset1.exists(), f"no such path exists: {self.dataset1}"
        assert self.dataset2.exists(), f"no such path exists: {self.dataset2}"
        assert not self.force_acid or self.acid != "autodetect", "cannot force-acid without specifying acid"
        assert self.from_position is None or self.from_position >= 0, "from-position must be at least 0"
        assert self.from_position is None or self.to_position is None or self.from_position <= self.to_position, "from-position must be below or equal to to-position"
        assert self.min_coverage >= 1, "min-coverage must be at least 1"
        assert self.resample_size >= 0, "resample-size must be at least 0"
        assert self.strand in (None, "+", "-"), "strand must be one of: '+' or '-'"
        assert self.to_position is None or self.to_position >= 0, "to-position must be at least zero"

