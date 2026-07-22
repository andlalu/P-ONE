from __future__ import annotations

from dataclasses import dataclass, field

from Models.Heston.parameters import HestonParameters
from OptionPricing.cos_basis import FixedCosBasisConfig


@dataclass(frozen=True)
class ImpliedStateConfig:
    cos_basis: FixedCosBasisConfig
    v_min: float = 1e-8
    v_max: float = 1.0
    tol: float = 2e-4
    max_iter: int = 40
    boundary_tol: float = 1e-5
    warm_start_window: float | None = 0.25
    state_solver: str = "golden_section"
    fallback_solver: str | None = "golden_section"
    finite_difference_relative_step: float = 1e-4
    finite_difference_absolute_step: float = 1e-7
    minimum_black_vega: float = 5e-6

    def validate(self) -> None:
        if self.v_min < 0.0 or self.v_max <= self.v_min:
            raise ValueError("variance bounds must satisfy 0 <= v_min < v_max")
        if self.tol <= 0.0 or self.max_iter <= 0:
            raise ValueError("tol and max_iter must be strictly positive")
        if self.boundary_tol < 0.0:
            raise ValueError("boundary_tol must be non-negative")
        if self.warm_start_window is not None and self.warm_start_window <= 0.0:
            raise ValueError("warm_start_window must be positive when provided")
        supported_solvers = {
            "golden_section",
            "bounded_brent",
            "least_squares_finite_difference",
            "least_squares_analytic",
        }
        if self.state_solver not in supported_solvers:
            raise ValueError(f"unsupported state_solver: {self.state_solver!r}")
        if self.fallback_solver not in {None, "golden_section", "bounded_brent"}:
            raise ValueError("fallback_solver must be None, 'golden_section', or 'bounded_brent'")
        if self.fallback_solver == self.state_solver:
            self_fallback_allowed = self.state_solver in {"golden_section", "bounded_brent"}
            if not self_fallback_allowed:
                raise ValueError("fallback_solver must differ from the selected derivative state solver")
        if self.finite_difference_relative_step <= 0.0 or self.finite_difference_absolute_step <= 0.0:
            raise ValueError("finite-difference Jacobian steps must be strictly positive")
        if self.minimum_black_vega <= 0.0:
            raise ValueError("minimum_black_vega must be strictly positive")


@dataclass(frozen=True)
class CcfQuadratureConfig:
    dimension: int = 2
    order: int = 3
    scale: float = 1.0

    def validate(self) -> None:
        if self.dimension <= 0 or self.order <= 0:
            raise ValueError("quadrature dimension and order must be positive")
        if self.scale <= 0.0:
            raise ValueError("quadrature scale must be strictly positive")


@dataclass(frozen=True)
class CgmmConfig:
    implied_state: ImpliedStateConfig
    quadrature: CcfQuadratureConfig = field(default_factory=CcfQuadratureConfig)
    instrument_precision: tuple[float, float] = (10.0, 50.0)
    transition_rk_steps: int = 32
    transition_cf_method: str = "analytic"
    dt: float | None = None
    spacing_tolerance: float = 1e-10

    def validate(self) -> None:
        if len(self.instrument_precision) != 2 or any(value <= 0.0 for value in self.instrument_precision):
            raise ValueError("instrument_precision must contain two positive entries")
        if self.transition_rk_steps <= 0:
            raise ValueError("transition_rk_steps must be positive")
        if self.transition_cf_method not in {"analytic", "rk4"}:
            raise ValueError("transition_cf_method must be 'analytic' or 'rk4'")
        if self.dt is not None and self.dt <= 0.0:
            raise ValueError("dt must be positive when provided")
        if self.spacing_tolerance <= 0.0:
            raise ValueError("spacing_tolerance must be positive")


@dataclass(frozen=True)
class PowellStageConfig:
    max_evaluations: int
    xtol: float
    ftol: float

    def validate(self) -> None:
        if self.max_evaluations <= 0 or self.xtol <= 0.0 or self.ftol <= 0.0:
            raise ValueError("Powell stage limits and tolerances must be strictly positive")


@dataclass(frozen=True)
class OptimizerConfig:
    base_start: HestonParameters
    natural_bounds: tuple[tuple[float, float], ...]
    candidate_relative_perturbations: tuple[tuple[float, ...], ...] = (
        (0.12, -0.12, 0.12, 0.12, 0.12, 0.12),
        (-0.12, 0.12, -0.12, -0.12, -0.12, -0.12),
    )
    stage1: PowellStageConfig = field(default_factory=lambda: PowellStageConfig(120, 2e-2, 2e-3))
    stage2: PowellStageConfig = field(default_factory=lambda: PowellStageConfig(300, 2e-4, 2e-5))
    penalty_value: float = 1e12
    progress_every: int = 10

    def validate(self) -> None:
        if len(self.natural_bounds) != 6:
            raise ValueError("natural_bounds must contain six (lower, upper) pairs")
        if any(len(pair) != 2 or pair[0] >= pair[1] for pair in self.natural_bounds):
            raise ValueError("each natural parameter bound must satisfy lower < upper")
        if len(self.candidate_relative_perturbations) > 2:
            raise ValueError("at most two perturbations are allowed in addition to the base start")
        if any(len(item) != 6 for item in self.candidate_relative_perturbations):
            raise ValueError("each candidate perturbation must contain six entries")
        if self.penalty_value <= 0.0:
            raise ValueError("penalty_value must be positive")
        if self.progress_every <= 0:
            raise ValueError("progress_every must be positive")
