from __future__ import annotations

import logging
import math
import time
from dataclasses import asdict
from typing import Any

import numpy as np

from Estimation.ISCGMM.cgmm_criterion import CgmmFirstStepCriterion
from Estimation.ISCGMM.config import CgmmConfig, PowellConfig
from Estimation.ISCGMM.parameter_transform import free_parameter_bounds, from_free, to_free
from Estimation.ISCGMM.results import FirstStepEstimate, PowellPassResult
from Models.Heston.parameters import HestonParameters
from OptionData.panel import OptionPanel

LOGGER = logging.getLogger(__name__)
EXPECTED_NUMERICAL_FAILURES = (ValueError, FloatingPointError, OverflowError, ArithmeticError, np.linalg.LinAlgError)


def _natural_coordinates(theta: HestonParameters) -> np.ndarray:
    return np.array([theta.eta, theta.kappa, theta.vbar, theta.sigma_v, theta.rho, theta.kappa_q], dtype=float)


def _theta_from_natural(values: np.ndarray, *, r: float, q: float) -> HestonParameters:
    eta, kappa, vbar, sigma_v, rho, kappa_q = (float(value) for value in values)
    theta = HestonParameters(
        eta=eta,
        kappa=kappa,
        vbar=vbar,
        sigma_v=sigma_v,
        rho=rho,
        eta_v=kappa - kappa_q,
        r=r,
        q=q,
    )
    theta.validate()
    return theta


def _candidate_starts(config: PowellConfig, *, r: float, q: float) -> tuple[HestonParameters, ...]:
    base = config.base_start.with_rates(r=r, q=q)
    natural = _natural_coordinates(base)
    lower = np.array([pair[0] for pair in config.natural_bounds], dtype=float)
    upper = np.array([pair[1] for pair in config.natural_bounds], dtype=float)
    if np.any(natural < lower) or np.any(natural > upper):
        raise ValueError("configured base start must lie within natural parameter bounds")
    starts = [base]
    for perturbation in config.candidate_relative_perturbations:
        candidate = np.clip(natural * (1.0 + np.asarray(perturbation, dtype=float)), lower, upper)
        starts.append(_theta_from_natural(candidate, r=r, q=q))
    if len(starts) > 3:
        raise ValueError("Powell start screening may contain no more than three candidates")
    return tuple(starts)


def estimate_first_step(
    panel: OptionPanel,
    *,
    criterion_config: CgmmConfig,
    powell_config: PowellConfig,
) -> FirstStepEstimate:
    """Estimate first-step IS-CGMM with coarse and refinement Powell passes."""

    from scipy.optimize import Bounds, minimize  # type: ignore[import-not-found]

    criterion_config.validate()
    powell_config.validate()
    started = time.perf_counter()
    LOGGER.info("estimation job started")
    LOGGER.info(
        "panel loaded sample=%s scenario=%s dates=%d contracts=%d",
        panel.metadata.get("sample_id", "unknown"),
        panel.metadata.get("scenario", "unknown"),
        panel.n_dates,
        panel.n_contracts,
    )

    try:
        criterion = CgmmFirstStepCriterion(panel, criterion_config)
        rate, dividend_yield = panel.first_rate_pair()
        free_bounds = free_parameter_bounds(powell_config.natural_bounds)
        lower = np.array([pair[0] for pair in free_bounds])
        upper = np.array([pair[1] for pair in free_bounds])
        evaluation_count = 0
        penalty_count = 0
        best_value = math.inf

        def objective(free: np.ndarray) -> float:
            nonlocal evaluation_count, penalty_count, best_value
            evaluation_count += 1
            try:
                value = float(criterion.evaluate_free(np.asarray(free), r=rate, q=dividend_yield))
                if not math.isfinite(value):
                    raise FloatingPointError("non-finite C-GMM criterion")
            except EXPECTED_NUMERICAL_FAILURES:
                penalty_count += 1
                value = float(powell_config.penalty_value)
            best_value = min(best_value, value)
            if evaluation_count % powell_config.progress_every == 0:
                LOGGER.info(
                    "Powell progress evaluations=%d best=%.8g elapsed=%.1fs penalties=%d",
                    evaluation_count,
                    best_value,
                    time.perf_counter() - started,
                    penalty_count,
                )
            return value

        starts = _candidate_starts(powell_config, r=rate, q=dividend_yield)
        screened: list[dict[str, Any]] = []
        for index, theta in enumerate(starts):
            value = objective(to_free(theta))
            screened.append({"index": index, "parameters": asdict(theta), "criterion": value})
        selected_index = int(np.argmin([item["criterion"] for item in screened]))
        selected = starts[selected_index]
        initial_criterion = float(screened[selected_index]["criterion"])
        LOGGER.info(
            "candidate-start screening completed criteria=%s selected=%d",
            [float(item["criterion"]) for item in screened],
            selected_index,
        )

        # These are two searches over the same first-step criterion.
        powell_passes: list[PowellPassResult] = []
        current = to_free(selected)
        scipy_bounds = Bounds(lower, upper)
        scipy_result = None
        for pass_name, pass_config in (
            ("coarse", powell_config.coarse_pass),
            ("refinement", powell_config.refinement_pass),
        ):
            LOGGER.info("Powell %s pass started", pass_name)
            before = evaluation_count
            scipy_result = minimize(
                objective,
                current,
                method="Powell",
                bounds=scipy_bounds,
                options={
                    "maxfev": pass_config.max_evaluations,
                    "xtol": pass_config.xtol,
                    "ftol": pass_config.ftol,
                    "disp": False,
                },
            )
            current = np.asarray(scipy_result.x, dtype=float)
            pass_result = PowellPassResult(
                pass_name=pass_name,
                success=bool(scipy_result.success),
                status=int(scipy_result.status),
                message=str(scipy_result.message),
                free_parameters=current.copy(),
                criterion=float(scipy_result.fun),
                iterations=int(scipy_result.nit),
                function_evaluations=evaluation_count - before,
            )
            powell_passes.append(pass_result)
            LOGGER.info(
                "Powell %s pass completed success=%s criterion=%.8g evaluations=%d",
                pass_name,
                pass_result.success,
                pass_result.criterion,
                pass_result.function_evaluations,
            )
            LOGGER.debug(
                "Powell %s pass best free parameters=%s",
                pass_name,
                current.tolist(),
            )

        assert scipy_result is not None
        estimated = from_free(current, r=rate, q=dividend_yield)
        final_diagnostics = criterion.evaluate(estimated, return_diagnostics=True)
        evaluation_count += 1
        LOGGER.info("final diagnostic evaluation completed criterion=%.8g", final_diagnostics.criterion_value)
        LOGGER.debug(
            "implied-state success=%.3f boundary-hit=%.3f timings=%s",
            final_diagnostics.implied_state.success_rate,
            final_diagnostics.implied_state.boundary_hit_rate,
            final_diagnostics.timings,
        )
        estimate = FirstStepEstimate(
            candidate_starts=tuple(screened),
            selected_start=selected,
            initial_criterion=initial_criterion,
            estimated_parameters=estimated,
            free_parameters=current,
            final_criterion=float(final_diagnostics.criterion_value),
            success=bool(scipy_result.success),
            status=int(scipy_result.status),
            message=str(scipy_result.message),
            iterations=sum(result.iterations for result in powell_passes),
            function_evaluations=evaluation_count,
            penalty_evaluations=penalty_count,
            powell_passes=tuple(powell_passes),
            total_runtime=time.perf_counter() - started,
            final_diagnostics=final_diagnostics,
            metadata=dict(final_diagnostics.metadata),
        )
        LOGGER.info("job completed success=%s runtime=%.1fs", estimate.success, estimate.total_runtime)
        return estimate
    except Exception:
        LOGGER.exception("estimation job failed")
        raise
