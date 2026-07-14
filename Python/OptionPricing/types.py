from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from OptionPricing.base import CcfSolver, OptionPricer


@dataclass(frozen=True)
class VarianceScaledCosConfig:
    """Explicit reference configuration for the legacy variance-scaled width."""

    n_cos: int = 256
    width_multiplier: float = 10.0

    def validate(self) -> None:
        if self.n_cos <= 1:
            raise ValueError("n_cos must be > 1")
        if self.width_multiplier <= 0.0:
            raise ValueError("width_multiplier must be strictly positive")


@dataclass(frozen=True)
class OptionPriceCubeConfig:
    strikes: np.ndarray
    maturities: np.ndarray
    rates: np.ndarray

    def validate(self) -> None:
        if self.strikes.ndim != 1 or np.any(self.strikes <= 0.0):
            raise ValueError("strikes must be positive 1D array")
        if self.maturities.ndim != 1 or np.any(self.maturities < 0.0):
            raise ValueError("maturities must be non-negative 1D array")
        if self.rates.ndim != 1 or self.rates.shape[0] != self.maturities.shape[0]:
            raise ValueError("rates must be 1D and same length as maturities")


@dataclass(frozen=True)
class CoefficientTensor:
    u_grid: np.ndarray
    maturities: np.ndarray
    cf_a: np.ndarray
    cf_b: np.ndarray


@dataclass(frozen=True)
class PreparedFixedCosBasis:
    """Reusable maturity-specific COS grid, payoff terms and affine A/B."""

    maturity: float
    effective_width: float
    n_cos: int
    u_grid: np.ndarray
    payoff_terms: np.ndarray
    coefficients: CoefficientTensor


@dataclass(frozen=True)
class FixedBasisPriceJacobian:
    """Prices and initial-variance derivatives for one fixed-basis maturity."""

    prices: np.ndarray
    initial_variance_jacobian: np.ndarray
    price_clipped: np.ndarray


@dataclass(frozen=True)
class OptionPriceCube:
    """Legacy dense vectorised pricing output; not a market-data panel."""

    prices: np.ndarray
    observation_index: np.ndarray
    strikes: np.ndarray
    maturities: np.ndarray
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class PricingStack:
    ccf_solver: CcfSolver
    option_pricer: OptionPricer
