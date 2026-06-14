from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from Estimation.ISCGMM.cgmm_criterion import CgmmFirstStepCriterion
from Estimation.ISCGMM.types import CgmmConfig, EstimationPanel, HestonEstimationParams, from_free, to_free


@dataclass(frozen=True)
class FirstStepEstimate:
    theta: HestonEstimationParams
    free_theta: np.ndarray
    criterion_value: float
    n_evaluations: int
    converged: bool


def estimate_first_step(
    panel: EstimationPanel,
    initial_theta: HestonEstimationParams,
    *,
    config: CgmmConfig | None = None,
    initial_step: np.ndarray | None = None,
    max_iter: int = 8,
    tol: float = 1e-3,
) -> FirstStepEstimate:
    """Deterministic coordinate-search fallback for the first-step objective.

    This is intentionally simple because the first deliverable focuses on a
    stable criterion and implied-state inversion. A SciPy optimizer can be
    plugged in later without changing the criterion object.
    """

    criterion = CgmmFirstStepCriterion(panel, config)
    current = to_free(initial_theta)
    step = np.full_like(current, 0.1, dtype=float) if initial_step is None else np.asarray(initial_step, dtype=float)
    if step.shape != current.shape:
        raise ValueError("initial_step must have shape (6,)")
    rate, dividend_yield = panel.first_rate_pair()

    def objective(free: np.ndarray) -> float:
        return criterion.evaluate_free(free, r=rate, q=dividend_yield)

    value = objective(current)
    n_eval = 1
    converged = False

    for _ in range(max_iter):
        improved = False
        for dim in range(current.shape[0]):
            candidates = []
            for sign in (-1.0, 1.0):
                trial = current.copy()
                trial[dim] += sign * step[dim]
                try:
                    trial_value = objective(trial)
                except ValueError:
                    trial_value = float("inf")
                candidates.append((trial_value, trial))
                n_eval += 1
            best_value, best_free = min(candidates, key=lambda item: item[0])
            if best_value + tol < value:
                current = best_free
                value = best_value
                improved = True
        if not improved:
            step *= 0.5
            if float(np.max(step)) < tol:
                converged = True
                break

    theta = from_free(current, r=rate, q=dividend_yield)
    return FirstStepEstimate(
        theta=theta,
        free_theta=current,
        criterion_value=float(value),
        n_evaluations=n_eval,
        converged=converged,
    )
