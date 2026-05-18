from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from OptionPricing.base import CcfSolver, OptionPricer


@dataclass(frozen=True)
class HestonPricingParamsQ:
    kappa: float
    vbar: float
    sigma_v: float
    rho: float
    r: float = 0.0
    q: float = 0.0
    v_risk_premium: float = 0.0

    def validate(self) -> None:
        if self.kappa <= 0.0:
            raise ValueError("kappa must be strictly positive")
        if self.vbar < 0.0:
            raise ValueError("vbar must be non-negative")
        if self.sigma_v < 0.0:
            raise ValueError("sigma_v must be non-negative")
        if not -1.0 <= self.rho <= 1.0:
            raise ValueError("rho must be in [-1, 1]")


@dataclass(frozen=True)
class CosPricingConfig:
    n_cos: int = 256
    truncation_width: float = 10.0

    def validate(self) -> None:
        if self.n_cos <= 1:
            raise ValueError("n_cos must be > 1")
        if self.truncation_width <= 0.0:
            raise ValueError("truncation_width must be strictly positive")


@dataclass(frozen=True)
class OptionPanelConfig:
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
class OptionPanel:
    prices: np.ndarray
    observation_index: np.ndarray
    strikes: np.ndarray
    maturities: np.ndarray
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class PricingStack:
    ccf_solver: CcfSolver
    option_pricer: OptionPricer
