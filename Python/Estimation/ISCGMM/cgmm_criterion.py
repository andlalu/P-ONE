from __future__ import annotations

import itertools
import math
import time
from dataclasses import asdict
from typing import Any

import numpy as np

from Estimation.ISCGMM.config import CcfQuadratureConfig, CgmmConfig
from Estimation.ISCGMM.heston_transition_cf import heston_p_transition_cf
from Estimation.ISCGMM.implied_state import imply_heston_variance_path
from Estimation.ISCGMM.parameter_transform import from_free
from Estimation.ISCGMM.results import CriterionDiagnostics
from Models.Heston.parameters import HestonParameters
from OptionData.panel import OptionPanel


def build_estimation_state(log_spot: np.ndarray, implied_variance: np.ndarray) -> np.ndarray:
    """Return X_i = (r_i, V_i), i=1,...,T, using r_i = logS_i - logS_{i-1}."""

    log_s = np.asarray(log_spot, dtype=float)
    variance = np.asarray(implied_variance, dtype=float)
    if log_s.ndim != 1 or variance.ndim != 1 or log_s.shape != variance.shape:
        raise ValueError("log_spot and implied_variance must be same-length 1D arrays")
    if log_s.shape[0] < 3:
        raise ValueError("at least three dates are required")
    return np.column_stack((np.diff(log_s), variance[1:]))


def make_cgmm_quadrature(config: CcfQuadratureConfig | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Create deterministic tensor quadrature nodes for CCF argument integration."""

    cfg = CcfQuadratureConfig() if config is None else config
    cfg.validate()
    one_d_nodes, one_d_weights = np.polynomial.hermite.hermgauss(cfg.order)
    one_d_nodes = math.sqrt(2.0) * cfg.scale * one_d_nodes
    one_d_weights = one_d_weights / math.sqrt(math.pi)

    nodes: list[tuple[float, ...]] = []
    weights: list[float] = []
    for multi_index in itertools.product(range(cfg.order), repeat=cfg.dimension):
        nodes.append(tuple(float(one_d_nodes[idx]) for idx in multi_index))
        weight = 1.0
        for idx in multi_index:
            weight *= float(one_d_weights[idx])
        weights.append(weight)
    return np.array(nodes, dtype=float), np.array(weights, dtype=float)


def _infer_constant_transition_interval(panel: OptionPanel, config: CgmmConfig) -> float:
    diffs = np.diff(panel.times)
    if diffs.size == 0 or np.any(diffs <= 0.0):
        raise ValueError("panel times must contain strictly positive consecutive intervals")
    interval = float(diffs[0])
    if not np.allclose(diffs, interval, rtol=0.0, atol=config.spacing_tolerance):
        raise ValueError("IS-CGMM requires equally spaced panel observations")
    if config.dt is not None and not np.isclose(config.dt, interval, rtol=0.0, atol=config.spacing_tolerance):
        raise ValueError(
            f"configured transition interval dt={config.dt} differs from inferred panel interval {interval}"
        )
    return interval


def _gaussian_instrument_kernel(x_prev: np.ndarray, precision: tuple[float, float]) -> np.ndarray:
    scales = np.asarray(precision, dtype=float)
    diffs = x_prev[:, None, :] - x_prev[None, :, :]
    exponent = -0.5 * np.sum((scales[None, None, :] * diffs) ** 2, axis=2)
    return np.exp(exponent)


def _variance_summary(values: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
    }


class CgmmFirstStepCriterion:
    """First-step Heston implied-state C-GMM criterion.

    The instrument coordinate is analytically integrated under a Gaussian
    weight, leaving numerical quadrature only over the CCF argument s.
    """

    def __init__(self, panel: OptionPanel, config: CgmmConfig):
        self.panel = panel
        self.config = config
        self.config.validate()
        self.s_nodes, self.s_weights = make_cgmm_quadrature(self.config.quadrature)
        if self.s_nodes.shape[1] != 2:
            raise ValueError("Heston IS-CGMM requires two-dimensional s nodes")

    def evaluate(self, theta: HestonParameters, *, return_diagnostics: bool = False) -> float | CriterionDiagnostics:
        panel_rate, panel_dividend_yield = self.panel.first_rate_pair()
        theta = theta.with_rates(r=panel_rate, q=panel_dividend_yield)
        theta.validate()

        t0 = time.perf_counter()
        implied = imply_heston_variance_path(theta, self.panel, self.config.implied_state)
        if np.any(implied.failed) or not np.all(np.isfinite(implied.objective)):
            raise FloatingPointError("implied-state inversion failed")
        t_state = time.perf_counter()

        x_state = build_estimation_state(self.panel.log_spot, implied.variance)
        x_prev = x_state[:-1]
        x_next = x_state[1:]
        dt = _infer_constant_transition_interval(self.panel, self.config)
        n_transitions = x_prev.shape[0]
        if n_transitions <= 0:
            raise ValueError("not enough dates for C-GMM transition residuals")

        phase = np.exp(1j * (self.s_nodes @ x_next.T))
        cond_cf = heston_p_transition_cf(
            self.s_nodes,
            x_prev,
            dt=dt,
            theta=theta,
            rk_steps=self.config.transition_rk_steps,
            method=self.config.transition_cf_method,
        )
        residuals = phase - cond_cf
        t_residuals = time.perf_counter()

        kernel = _gaussian_instrument_kernel(x_prev, self.config.instrument_precision)
        t_kernel = time.perf_counter()

        total = 0.0 + 0.0j
        for weight, residual in zip(self.s_weights, residuals):
            total += weight * np.sum(kernel * residual[:, None] * np.conjugate(residual[None, :]))
        value = float(np.real(total) / float(n_transitions * n_transitions))
        if value < 0.0 and value > -1e-12:
            value = 0.0
        if not np.isfinite(value):
            value = float("inf")
        t_done = time.perf_counter()

        if not return_diagnostics:
            return value

        true_v = self.panel.true_variance
        rmse = None
        mae = None
        if true_v is not None:
            diff = implied.variance - true_v
            rmse = float(np.sqrt(np.mean(diff * diff)))
            mae = float(np.mean(np.abs(diff)))

        return CriterionDiagnostics(
            criterion_value=value,
            n_dates=self.panel.n_dates,
            n_transitions=n_transitions,
            n_contracts=self.panel.n_contracts,
            implied_state=implied,
            implied_variance_summary=_variance_summary(implied.variance),
            true_variance_rmse=rmse,
            true_variance_mae=mae,
            timings={
                "state_inversion": t_state - t0,
                "ccf_residuals": t_residuals - t_state,
                "kernel": t_kernel - t_residuals,
                "quadratic_form": t_done - t_kernel,
                "total": t_done - t0,
            },
            metadata={
                "cos_basis": {
                    **dict(self.panel.metadata["cos_basis"]),
                    "estimation_n_cos": self.config.implied_state.cos_basis.estimation_n_cos,
                }
            },
        )

    def evaluate_free(self, free_theta: np.ndarray, *, r: float, q: float) -> float:
        theta = from_free(free_theta, r=r, q=q)
        return float(self.evaluate(theta))


def criterion_diagnostics_to_dict(theta: HestonParameters, diagnostics: CriterionDiagnostics) -> dict[str, Any]:
    implied = diagnostics.implied_state
    return {
        "theta": asdict(theta),
        "criterion_value": diagnostics.criterion_value,
        "n_dates": diagnostics.n_dates,
        "n_transitions": diagnostics.n_transitions,
        "n_contracts": diagnostics.n_contracts,
        "state_inversion": {
            "success_rate": implied.success_rate,
            "boundary_hit_rate": implied.boundary_hit_rate,
            "nfev_total": int(np.sum(implied.nfev)),
            "objective_mean": float(np.mean(implied.objective)),
            "objective_max": float(np.max(implied.objective)),
            "coefficient_solve_count": int(implied.coefficient_solve_count),
            "coefficient_cache_hits": int(implied.coefficient_cache_hits),
            "fixed_coefficient_count": int(implied.fixed_coefficient_count),
            "solver_name": implied.solver_name,
            "fallback_count": implied.fallback_count,
            "jacobian_evaluation_count": implied.jacobian_evaluation_count,
        },
        "implied_variance": diagnostics.implied_variance_summary,
        "true_variance_rmse": diagnostics.true_variance_rmse,
        "true_variance_mae": diagnostics.true_variance_mae,
        "timings": diagnostics.timings,
        "metadata": diagnostics.metadata,
    }
