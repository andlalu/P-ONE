from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np


@dataclass(frozen=True)
class FixedCosBasisConfig:
    """Maturity-specific COS widths and role-specific term counts."""

    maturities: tuple[float, ...]
    effective_widths: tuple[float, ...]
    generation_n_cos: int
    estimation_n_cos: int
    maturity_tolerance: float = 1e-10

    def validate(self) -> None:
        maturities = np.asarray(self.maturities, dtype=float)
        effective_widths = np.asarray(self.effective_widths, dtype=float)
        if maturities.ndim != 1 or effective_widths.ndim != 1 or maturities.size == 0:
            raise ValueError("maturities and effective_widths must be non-empty aligned 1D sequences")
        if maturities.shape != effective_widths.shape:
            raise ValueError("maturities and effective_widths must have a one-to-one match")
        if not np.all(np.isfinite(maturities)) or np.any(maturities <= 0.0):
            raise ValueError("fixed COS maturities must be finite and strictly positive")
        if not np.all(np.isfinite(effective_widths)) or np.any(effective_widths <= 0.0):
            raise ValueError("fixed COS effective widths must be finite and strictly positive")
        if self.generation_n_cos <= 1 or self.estimation_n_cos <= 1:
            raise ValueError("generation_n_cos and estimation_n_cos must both be > 1")
        if not np.isfinite(self.maturity_tolerance) or self.maturity_tolerance <= 0.0:
            raise ValueError("maturity_tolerance must be finite and strictly positive")
        ordered_maturities = np.sort(maturities)
        if ordered_maturities.size > 1 and np.any(np.diff(ordered_maturities) <= self.maturity_tolerance):
            raise ValueError("fixed COS maturity entries must be unique within maturity_tolerance")

    def width_for_maturity(self, maturity: float) -> float:
        distances = np.abs(np.asarray(self.maturities, dtype=float) - float(maturity))
        matches = np.flatnonzero(distances <= self.maturity_tolerance)
        if matches.size == 0:
            raise ValueError(f"maturity {maturity!r} is not present in the fixed COS basis")
        if matches.size > 1:
            raise ValueError(f"maturity {maturity!r} matches multiple fixed COS basis entries")
        return float(self.effective_widths[int(matches[0])])

    def validate_requested_maturities(self, maturities: Iterable[float]) -> None:
        requested_maturities = np.asarray(tuple(maturities), dtype=float).reshape(-1)
        if requested_maturities.size == 0:
            raise ValueError("requested maturities must be non-empty")
        for maturity in requested_maturities:
            self.width_for_maturity(float(maturity))


def cos_specification_metadata(basis: FixedCosBasisConfig) -> dict[str, Any]:
    """Return the numerical COS specification stored beside generated panels."""

    return {
        "maturities": [float(maturity) for maturity in basis.maturities],
        "effective_widths": [float(width) for width in basis.effective_widths],
        "generation_n_cos": int(basis.generation_n_cos),
        "estimation_n_cos": int(basis.estimation_n_cos),
    }


def validate_panel_cos_compatibility(
    basis: FixedCosBasisConfig,
    panel_metadata: dict[str, Any],
) -> None:
    """Check the generation settings that must match during estimation."""

    panel_maturities = tuple(float(value) for value in panel_metadata.get("maturities", ()))
    panel_widths = tuple(float(value) for value in panel_metadata.get("effective_widths", ()))
    if len(panel_maturities) != len(basis.maturities) or len(panel_widths) != len(basis.effective_widths):
        raise ValueError("panel and estimator COS maturity grids differ")
    if int(panel_metadata.get("generation_n_cos", 0)) != basis.generation_n_cos:
        raise ValueError("panel generation_n_cos differs from the experiment configuration")
    for maturity, width in zip(panel_maturities, panel_widths):
        try:
            configured_width = basis.width_for_maturity(maturity)
        except ValueError as error:
            raise ValueError("panel and estimator COS maturity grids differ") from error
        if not np.isclose(configured_width, width, rtol=0.0, atol=basis.maturity_tolerance):
            raise ValueError(f"COS effective width mismatch at maturity {maturity}")
    basis.validate_requested_maturities(panel_maturities)
