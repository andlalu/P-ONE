from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from Models.Heston.parameters import HestonParameters


def _plain(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


@dataclass(frozen=True)
class ImpliedStateResult:
    variance: np.ndarray
    objective: np.ndarray
    boundary_hit: np.ndarray
    failed: np.ndarray
    nfev: np.ndarray
    start_values: np.ndarray
    coefficient_solve_count: int = 0
    coefficient_cache_hits: int = 0
    fixed_coefficient_count: int = 0

    @property
    def success_rate(self) -> float:
        return float(np.mean(~self.failed))

    @property
    def boundary_hit_rate(self) -> float:
        return float(np.mean(self.boundary_hit))

    def to_dict(self) -> dict[str, Any]:
        return _plain(asdict(self))


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
    metadata: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        return _plain(value)


@dataclass(frozen=True)
class OptimizerStageResult:
    stage: int
    success: bool
    status: int
    message: str
    free_parameters: np.ndarray
    criterion: float
    iterations: int
    function_evaluations: int

    def to_dict(self) -> dict[str, Any]:
        return _plain(asdict(self))


@dataclass(frozen=True)
class FirstStepEstimate:
    candidate_starts: tuple[dict[str, Any], ...]
    selected_start: HestonParameters
    initial_criterion: float
    estimated_parameters: HestonParameters
    free_parameters: np.ndarray
    final_criterion: float
    success: bool
    status: int
    message: str
    iterations: int
    function_evaluations: int
    penalty_evaluations: int
    stage_results: tuple[OptimizerStageResult, ...]
    total_runtime: float
    final_diagnostics: CriterionDiagnostics
    metadata: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return _plain(asdict(self))
