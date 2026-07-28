from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class FixedCosBasisConfig:
    """Fixed COS widths and term counts for each maturity."""

    maturities: tuple[float, ...]
    effective_widths: tuple[float, ...]
    generation_n_cos: int
    estimation_n_cos: int
    maturity_tolerance: float = 1e-10

    def validate(self) -> None:
        maturities = np.asarray(self.maturities, dtype=float)
        widths = np.asarray(self.effective_widths, dtype=float)
        if maturities.ndim != 1 or widths.ndim != 1 or maturities.size == 0:
            raise ValueError("maturities and effective_widths must be non-empty aligned 1D sequences")
        if maturities.shape != widths.shape:
            raise ValueError("maturities and effective_widths must have a one-to-one match")
        if not np.all(np.isfinite(maturities)) or np.any(maturities <= 0.0):
            raise ValueError("fixed COS maturities must be finite and strictly positive")
        if not np.all(np.isfinite(widths)) or np.any(widths <= 0.0):
            raise ValueError("fixed COS effective widths must be finite and strictly positive")
        if self.generation_n_cos <= 0 or self.estimation_n_cos <= 0:
            raise ValueError("generation_n_cos and estimation_n_cos must be positive")
        if not np.isfinite(self.maturity_tolerance) or self.maturity_tolerance <= 0.0:
            raise ValueError("maturity_tolerance must be finite and strictly positive")
        ordered = np.sort(maturities)
        if ordered.size > 1 and np.any(np.diff(ordered) <= self.maturity_tolerance):
            raise ValueError("fixed COS maturity entries must be unique within maturity_tolerance")

    def width_for_maturity(self, maturity: float) -> float:
        distances = np.abs(np.asarray(self.maturities, dtype=float) - float(maturity))
        matches = np.flatnonzero(distances <= self.maturity_tolerance)
        if matches.size == 0:
            raise ValueError(f"maturity {maturity!r} is not present in the fixed COS basis")
        return float(self.effective_widths[int(matches[0])])

    def validate_requested_maturities(self, maturities: Iterable[float]) -> None:
        requested = tuple(float(maturity) for maturity in maturities)
        if not requested:
            raise ValueError("requested maturities must be non-empty")
        for maturity in requested:
            self.width_for_maturity(maturity)
