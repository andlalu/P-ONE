from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np

BASIS_CONVENTION = "fixed_effective_width"
PRICING_IMPLEMENTATION = "heston_cos_fixed_basis_v1"


@dataclass(frozen=True)
class FixedCosBasisConfig:
    """Maturity-aligned fixed effective COS widths.

    Widths are effective log-price truncation half-widths, not multipliers.
    They never depend on simulated or candidate variance.
    """

    maturities: tuple[float, ...]
    effective_widths: tuple[float, ...]
    n_cos: int
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
        if self.n_cos <= 1:
            raise ValueError("n_cos must be > 1")
        if not np.isfinite(self.maturity_tolerance) or self.maturity_tolerance <= 0.0:
            raise ValueError("maturity_tolerance must be finite and strictly positive")
        ordered = np.sort(maturities)
        if ordered.size > 1 and np.any(np.diff(ordered) <= self.maturity_tolerance):
            raise ValueError("fixed COS maturity entries must be unique within maturity_tolerance")

    def width_for_maturity(self, maturity: float) -> float:
        self.validate()
        distances = np.abs(np.asarray(self.maturities, dtype=float) - float(maturity))
        matches = np.flatnonzero(distances <= self.maturity_tolerance)
        if matches.size == 0:
            raise ValueError(f"maturity {maturity!r} is not present in the fixed COS basis")
        if matches.size > 1:
            raise ValueError(f"maturity {maturity!r} matches multiple fixed COS basis entries")
        return float(self.effective_widths[int(matches[0])])

    def validate_requested_maturities(self, maturities: Iterable[float]) -> None:
        for maturity in np.asarray(tuple(maturities), dtype=float).reshape(-1):
            self.width_for_maturity(float(maturity))

    def to_metadata(self) -> dict[str, Any]:
        self.validate()
        return {
            "basis_convention": BASIS_CONVENTION,
            "maturities": [float(value) for value in self.maturities],
            "effective_widths": [float(value) for value in self.effective_widths],
            "n_cos": int(self.n_cos),
            "maturity_tolerance": float(self.maturity_tolerance),
            "pricing_implementation": PRICING_IMPLEMENTATION,
        }

    @classmethod
    def from_metadata(cls, metadata: dict[str, Any]) -> FixedCosBasisConfig:
        if metadata.get("basis_convention") != BASIS_CONVENTION:
            raise ValueError("panel does not use the fixed effective-width COS convention")
        if metadata.get("pricing_implementation") != PRICING_IMPLEMENTATION:
            raise ValueError("panel COS pricing implementation is unsupported")
        config = cls(
            maturities=tuple(float(value) for value in metadata["maturities"]),
            effective_widths=tuple(float(value) for value in metadata["effective_widths"]),
            n_cos=int(metadata["n_cos"]),
            maturity_tolerance=float(metadata.get("maturity_tolerance", 1e-10)),
        )
        config.validate()
        return config

    def assert_matches_metadata(self, metadata: dict[str, Any]) -> None:
        panel_config = self.from_metadata(metadata)
        self.validate()
        if self.n_cos != panel_config.n_cos:
            raise ValueError(f"COS basis n_cos mismatch: panel={panel_config.n_cos}, estimator={self.n_cos}")
        self.validate_requested_maturities(panel_config.maturities)
        panel_config.validate_requested_maturities(self.maturities)
        for maturity in self.maturities:
            expected = self.width_for_maturity(maturity)
            actual = panel_config.width_for_maturity(maturity)
            if not np.isclose(expected, actual, rtol=0.0, atol=max(self.maturity_tolerance, panel_config.maturity_tolerance)):
                raise ValueError(
                    f"COS effective width mismatch at maturity {maturity}: panel={actual}, estimator={expected}"
                )

