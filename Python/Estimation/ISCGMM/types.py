from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

from DGPSimulation.types import HestonParamsP
from OptionPricing.types import CosPricingConfig
from OptionPricing.types import HestonPricingParamsQ


@dataclass(frozen=True)
class HestonEstimationParams:
    """Current Heston parameter vector estimated by implied-state C-GMM."""

    eta: float
    kappa: float
    vbar: float
    sigma_v: float
    rho: float
    eta_v: float
    r: float = 0.0
    q: float = 0.0

    @property
    def kappa_q(self) -> float:
        return self.kappa - self.eta_v

    @property
    def vbar_q(self) -> float:
        return self.kappa * self.vbar / self.kappa_q

    def validate(self) -> None:
        if not math.isfinite(self.eta):
            raise ValueError("eta must be finite")
        if self.kappa <= 0.0:
            raise ValueError("kappa must be strictly positive")
        if self.vbar <= 0.0:
            raise ValueError("vbar must be strictly positive")
        if self.sigma_v <= 0.0:
            raise ValueError("sigma_v must be strictly positive")
        if not -1.0 < self.rho < 1.0:
            raise ValueError("rho must be in (-1, 1)")
        if self.kappa_q <= 0.0:
            raise ValueError("kappa - eta_v must be strictly positive")

    def to_pricing_q(self) -> HestonPricingParamsQ:
        self.validate()
        return HestonPricingParamsQ(
            kappa=self.kappa_q,
            vbar=self.vbar_q,
            sigma_v=self.sigma_v,
            rho=self.rho,
            r=self.r,
            q=self.q,
        )

    def to_physical_p(self) -> HestonParamsP:
        self.validate()
        return HestonParamsP(
            eta=self.eta,
            kappa=self.kappa,
            vbar=self.vbar,
            sigma_v=self.sigma_v,
            rho=self.rho,
            r=self.r,
            q=self.q,
        )

    def with_rates(self, *, r: float, q: float) -> HestonEstimationParams:
        return replace(self, r=float(r), q=float(q))


def to_free(theta: HestonEstimationParams) -> np.ndarray:
    """Map constrained Heston parameters to a cheap non-periodic free vector."""

    theta.validate()
    return np.array(
        [
            theta.eta,
            math.log(theta.kappa),
            math.log(theta.vbar),
            math.log(theta.sigma_v),
            np.arctanh(theta.rho),
            math.log(theta.kappa_q),
        ],
        dtype=float,
    )


def from_free(free_theta: np.ndarray, *, r: float = 0.0, q: float = 0.0) -> HestonEstimationParams:
    """Map [eta, log kappa, log vbar, log sigma_v, atanh rho, log kappa_q]."""

    free = np.asarray(free_theta, dtype=float)
    if free.shape != (6,):
        raise ValueError("free_theta must have shape (6,)")
    eta = float(free[0])
    kappa = float(np.exp(free[1]))
    vbar = float(np.exp(free[2]))
    sigma_v = float(np.exp(free[3]))
    rho = float(np.tanh(free[4]))
    kappa_q = float(np.exp(free[5]))
    theta = HestonEstimationParams(
        eta=eta,
        kappa=kappa,
        vbar=vbar,
        sigma_v=sigma_v,
        rho=rho,
        eta_v=kappa - kappa_q,
        r=float(r),
        q=float(q),
    )
    theta.validate()
    return theta


@dataclass(frozen=True)
class PanelDate:
    date_index: int
    time: float
    log_spot: float
    spot: float
    strikes: np.ndarray
    maturities: np.ndarray
    option_types: np.ndarray
    observed_iv: np.ndarray
    rates: np.ndarray
    dividend_yields: np.ndarray
    true_variance: float | None = None
    observed_price: np.ndarray | None = None
    clean_iv: np.ndarray | None = None

    @property
    def n_contracts(self) -> int:
        return int(self.observed_iv.shape[0])


@dataclass(frozen=True)
class EstimationPanel:
    dates: tuple[PanelDate, ...]
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if len(self.dates) < 3:
            raise ValueError("at least three panel dates are required for transition moments")

    @property
    def log_spot(self) -> np.ndarray:
        return np.array([item.log_spot for item in self.dates], dtype=float)

    @property
    def spot(self) -> np.ndarray:
        return np.array([item.spot for item in self.dates], dtype=float)

    @property
    def times(self) -> np.ndarray:
        return np.array([item.time for item in self.dates], dtype=float)

    @property
    def true_variance(self) -> np.ndarray | None:
        values = [item.true_variance for item in self.dates]
        if any(value is None for value in values):
            return None
        return np.array(values, dtype=float)

    @property
    def n_contracts(self) -> int:
        return int(sum(item.n_contracts for item in self.dates))

    @property
    def n_dates(self) -> int:
        return len(self.dates)

    def first_rate_pair(self) -> tuple[float, float]:
        first = self.dates[0]
        return float(np.mean(first.rates)), float(np.mean(first.dividend_yields))

    def truncate_dates(self, max_dates: int | None) -> EstimationPanel:
        if max_dates is None or max_dates >= len(self.dates):
            return self
        if max_dates < 3:
            raise ValueError("max_dates must be at least 3")
        return EstimationPanel(dates=self.dates[:max_dates], metadata=dict(self.metadata))


@dataclass(frozen=True)
class ImpliedStateConfig:
    v_min: float = 1e-8
    v_max: float = 1.0
    tol: float = 2e-4
    max_iter: int = 40
    boundary_tol: float = 1e-5
    warm_start_window: float | None = 0.25
    effective_truncation_width: float | None = None
    pricing_config: CosPricingConfig = field(default_factory=lambda: CosPricingConfig(n_cos=256, truncation_width=10.0))

    def validate(self) -> None:
        if self.v_min < 0.0 or self.v_max <= self.v_min:
            raise ValueError("variance bounds must satisfy 0 <= v_min < v_max")
        if self.tol <= 0.0:
            raise ValueError("tol must be strictly positive")
        if self.max_iter <= 0:
            raise ValueError("max_iter must be positive")
        if self.boundary_tol < 0.0:
            raise ValueError("boundary_tol must be non-negative")
        if self.warm_start_window is not None and self.warm_start_window <= 0.0:
            raise ValueError("warm_start_window must be positive when provided")
        if self.effective_truncation_width is not None and self.effective_truncation_width < 0.5:
            raise ValueError("effective_truncation_width must be at least the COS floor of 0.5")
        self.pricing_config.validate()


@dataclass(frozen=True)
class ImpliedStateResult:
    variance: np.ndarray
    objective: np.ndarray
    boundary_hit: np.ndarray
    failed: np.ndarray
    nfev: np.ndarray
    start_values: np.ndarray

    @property
    def success_rate(self) -> float:
        return float(np.mean(~self.failed))

    @property
    def boundary_hit_rate(self) -> float:
        return float(np.mean(self.boundary_hit))


@dataclass(frozen=True)
class QuadratureConfig:
    dim: int = 2
    order: int = 3
    scheme: str = "gauss_hermite"
    scale: float = 1.0

    def validate(self) -> None:
        if self.dim <= 0:
            raise ValueError("dim must be positive")
        if self.order <= 0:
            raise ValueError("order must be positive")
        if self.scheme not in {"gauss_hermite", "gauss_legendre"}:
            raise ValueError("unsupported quadrature scheme")
        if self.scale <= 0.0:
            raise ValueError("scale must be strictly positive")


@dataclass(frozen=True)
class CgmmConfig:
    implied_state: ImpliedStateConfig = field(default_factory=ImpliedStateConfig)
    quadrature: QuadratureConfig = field(default_factory=QuadratureConfig)
    instrument_precision: tuple[float, float] = (10.0, 50.0)
    transition_rk_steps: int = 32
    dt: float | None = None
    second_step: bool = False

    def validate(self) -> None:
        self.implied_state.validate()
        self.quadrature.validate()
        precision = np.asarray(self.instrument_precision, dtype=float)
        if precision.shape != (2,) or np.any(precision <= 0.0):
            raise ValueError("instrument_precision must contain two positive entries")
        if self.transition_rk_steps <= 0:
            raise ValueError("transition_rk_steps must be positive")
        if self.dt is not None and self.dt <= 0.0:
            raise ValueError("dt must be positive when provided")
        if self.second_step:
            raise NotImplementedError("second-step C-GMM is intentionally not implemented yet")


@dataclass(frozen=True)
class CriterionDiagnostics:
    criterion_value: float
    n_dates: int
    n_transitions: int
    n_contracts: int
    implied_state: ImpliedStateResult
    implied_variance_summary: dict[str, float]
    true_variance_rmse: float | None
    true_variance_mae: float | None
    timings: dict[str, float]


# Backwards-compatible aliases for code written against the first draft names.
HestonISParams = HestonEstimationParams
CGMMConfig = CgmmConfig
OptionDateSlice = PanelDate
OptionPanelData = EstimationPanel
