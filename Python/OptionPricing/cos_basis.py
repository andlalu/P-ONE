from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from OptionPricing.config import FixedCosBasisConfig
from OptionPricing.heston_ccf import CoefficientTensor


@dataclass(frozen=True)
class PreparedFixedCosBasis:
    """Everything needed to reuse one fixed COS maturity.

    It keeps the frequency grid, payoff terms, affine ``A`` and ``B``, width
    and term count together. Build it once, then reuse it across dates,
    candidate variances and state-solver calls. The derivative is then cheap:
    ``d exp(A + B v) / d v = B exp(A + B v)``.
    """

    maturity: float
    effective_width: float
    n_cos: int
    u_grid: np.ndarray
    payoff_terms: np.ndarray
    coefficients: CoefficientTensor


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
