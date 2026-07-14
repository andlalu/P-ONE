from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from Estimation.ISCGMM.config import ImpliedStateConfig
from Estimation.ISCGMM.implied_state_jacobian import (
    finite_difference_iv_variance_jacobian,
    price_iv_and_initial_variance_jacobian_fixed_basis,
)
from Estimation.ISCGMM.results import ImpliedStateResult
from ImpliedVolatility.black_iv import implied_vol_black76
from Models.Heston.parameters import HestonParameters, HestonRiskNeutralParameters
from OptionData.panel import OptionPanel, OptionPanelDate
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.types import PreparedFixedCosBasis


@dataclass(frozen=True)
class _ImpliedVarianceSolveResult:
    x: float
    fun: float
    nfev: int
    nit: int
    success: bool
    njev: int = 0
    solver_name: str = "golden_section"
    fallback_used: bool = False


@dataclass
class _PricingCacheDiagnostics:
    coefficient_solve_count: int = 0
    coefficient_cache_hits: int = 0
    fixed_coefficient_count: int = 0


@dataclass(frozen=True)
class _MaturityPricingBatch:
    maturity: float
    indices: np.ndarray
    rate: float
    dividend_yield: float
    strikes: np.ndarray
    option_types: np.ndarray
    option_types_for_pricing: np.ndarray
    forwards: np.ndarray
    discounts: np.ndarray
    basis: PreparedFixedCosBasis


@dataclass(frozen=True)
class _DateImpliedStatePricingContext:
    date: OptionPanelDate
    maturity_batches: tuple[_MaturityPricingBatch, ...]


def imply_heston_variance_path(
    theta: HestonParameters,
    panel: OptionPanel,
    config: ImpliedStateConfig,
) -> ImpliedStateResult:
    """Imply the scalar variance path from cross-sectional IV fit losses.

    Each date is grouped by maturity so one affine solve and one fixed COS
    payoff expansion serve all strikes at that maturity. The fixed basis also
    makes a future semi-analytical initial-variance sensitivity well defined:
    its affine dependence is ``exp(A + B v)`` and does not hide a derivative
    of the truncation rule.
    """

    config.validate()
    theta.validate()
    if panel.n_dates < 3:
        raise ValueError("at least three panel dates are required for transition moments")
    basis_metadata = panel.metadata.get("cos_basis")
    if not isinstance(basis_metadata, dict):
        raise ValueError("option panel is missing required fixed COS-basis metadata")
    config.cos_basis.assert_panel_compatible(basis_metadata)

    pricer = CosOptionPricer()
    params_q = theta.to_risk_neutral()
    cache_diagnostics = _PricingCacheDiagnostics()
    prepared_by_maturity = _prepare_fixed_bases(
        panel=panel,
        params_q=params_q,
        pricer=pricer,
        config=config,
        diagnostics=cache_diagnostics,
    )
    contexts = tuple(
        _prepare_date_pricing_context(
            date,
            prepared_by_maturity=prepared_by_maturity,
            config=config,
            diagnostics=cache_diagnostics,
        )
        for date in panel.dates
    )

    variance = np.empty(panel.n_dates, dtype=float)
    objective = np.empty(panel.n_dates, dtype=float)
    boundary_hit = np.zeros(panel.n_dates, dtype=bool)
    failed = np.zeros(panel.n_dates, dtype=bool)
    nfev = np.zeros(panel.n_dates, dtype=int)
    start_values = np.empty(panel.n_dates, dtype=float)
    jacobian_evaluation_count = 0
    fallback_count = 0

    previous = min(max(theta.vbar, config.v_min), config.v_max)
    for date_index, context in enumerate(contexts):
        start_values[date_index] = previous
        result = _imply_variance_for_date(
            context=context,
            start_value=previous,
            pricer=pricer,
            config=config,
        )
        variance[date_index] = min(max(result.x, config.v_min), config.v_max)
        objective[date_index] = result.fun
        nfev[date_index] = result.nfev
        failed[date_index] = not result.success
        jacobian_evaluation_count += result.njev
        fallback_count += int(result.fallback_used)
        boundary_hit[date_index] = (
            variance[date_index] <= config.v_min + config.boundary_tol
            or variance[date_index] >= config.v_max - config.boundary_tol
        )
        previous = variance[date_index]

    return ImpliedStateResult(
        variance=variance,
        objective=objective,
        boundary_hit=boundary_hit,
        failed=failed,
        nfev=nfev,
        start_values=start_values,
        coefficient_solve_count=cache_diagnostics.coefficient_solve_count,
        coefficient_cache_hits=cache_diagnostics.coefficient_cache_hits,
        fixed_coefficient_count=cache_diagnostics.fixed_coefficient_count,
        solver_name=config.state_solver,
        fallback_count=fallback_count,
        jacobian_evaluation_count=jacobian_evaluation_count,
    )


def _prepare_fixed_bases(
    *,
    panel: OptionPanel,
    params_q: HestonRiskNeutralParameters,
    pricer: CosOptionPricer,
    config: ImpliedStateConfig,
    diagnostics: _PricingCacheDiagnostics,
) -> dict[float, PreparedFixedCosBasis]:
    requested = np.unique(np.concatenate([date.maturities for date in panel.dates]))
    config.cos_basis.validate_requested_maturities(requested)
    prepared: dict[float, PreparedFixedCosBasis] = {}
    for maturity in requested:
        maturity_value = float(maturity)
        prepared[maturity_value] = pricer.prepare_fixed_basis(
            maturity=maturity_value,
            effective_width=config.cos_basis.width_for_maturity(maturity_value),
            n_cos=config.cos_basis.estimation_n_cos,
            model_params=params_q,
        )
        diagnostics.coefficient_solve_count += 1
        diagnostics.fixed_coefficient_count += 1
    return prepared


def _prepare_date_pricing_context(
    date: OptionPanelDate,
    *,
    prepared_by_maturity: dict[float, PreparedFixedCosBasis],
    config: ImpliedStateConfig,
    diagnostics: _PricingCacheDiagnostics,
) -> _DateImpliedStatePricingContext:
    """Group a date by maturity for vectorised strikes and basis reuse."""

    batches: list[_MaturityPricingBatch] = []
    for maturity in np.unique(date.maturities):
        indices = np.flatnonzero(np.isclose(date.maturities, maturity, rtol=0.0, atol=config.cos_basis.maturity_tolerance))
        rates = date.rates[indices]
        dividend_yields = date.dividend_yields[indices]
        if not np.allclose(rates, rates[0], rtol=0.0, atol=1e-14):
            raise ValueError("rates must be constant within each date/maturity batch")
        if not np.allclose(dividend_yields, dividend_yields[0], rtol=0.0, atol=1e-14):
            raise ValueError("dividend yields must be constant within each date/maturity batch")
        basis_key = min(prepared_by_maturity, key=lambda value: abs(value - float(maturity)))
        basis = prepared_by_maturity[basis_key]
        diagnostics.coefficient_cache_hits += 1
        rate = float(rates[0])
        dividend_yield = float(dividend_yields[0])
        maturities = date.maturities[indices].astype(float)
        batches.append(
            _MaturityPricingBatch(
                maturity=float(maturity),
                indices=indices,
                rate=rate,
                dividend_yield=dividend_yield,
                strikes=date.strikes[indices].astype(float),
                option_types=date.option_types[indices].astype(str),
                option_types_for_pricing=date.option_types[indices].astype(str).reshape(1, -1, 1),
                forwards=date.spot * np.exp((rate - dividend_yield) * maturities),
                discounts=np.exp(-rate * maturities),
                basis=basis,
            )
        )
    return _DateImpliedStatePricingContext(date=date, maturity_batches=tuple(batches))


def _compute_model_implied_volatilities(
    *,
    variance: float,
    context: _DateImpliedStatePricingContext,
    pricer: CosOptionPricer,
) -> np.ndarray:
    model_iv = np.empty(context.date.n_contracts, dtype=float)
    for batch in context.maturity_batches:
        prices = pricer.price_matrix_fixed_basis(
            log_s=np.array([context.date.log_spot]),
            variance=np.array([variance]),
            strike_grid=batch.strikes,
            rate=batch.rate,
            dividend_yield=batch.dividend_yield,
            basis=batch.basis,
            option_type=batch.option_types_for_pricing,
        )[0, :, 0]
        for local_position, price, forward, strike, discount, option_type in zip(
            batch.indices,
            prices,
            batch.forwards,
            batch.strikes,
            batch.discounts,
            batch.option_types,
        ):
            try:
                model_iv[local_position] = implied_vol_black76(
                    price=float(price),
                    forward=float(forward),
                    strike=float(strike),
                    tau=batch.maturity,
                    discount_factor=float(discount),
                    option_type=str(option_type),
                    on_bounds="clip",
                )
            except (ValueError, FloatingPointError):
                model_iv[local_position] = np.nan
    return model_iv


def _compute_model_iv_and_initial_variance_jacobian(
    *,
    variance: float,
    context: _DateImpliedStatePricingContext,
    config: ImpliedStateConfig,
    jacobian_method: str,
) -> tuple[np.ndarray, np.ndarray]:
    model_iv = np.empty(context.date.n_contracts, dtype=float)
    initial_variance_jacobian = np.empty(context.date.n_contracts, dtype=float)
    for maturity_batch in context.maturity_batches:
        if jacobian_method == "analytic":
            derivative_output = price_iv_and_initial_variance_jacobian_fixed_basis(
                log_spot=context.date.log_spot,
                implied_variance=variance,
                strikes=maturity_batch.strikes,
                option_types=maturity_batch.option_types,
                rate=maturity_batch.rate,
                dividend_yield=maturity_batch.dividend_yield,
                fixed_cos_basis=maturity_batch.basis,
                minimum_black_vega=config.minimum_black_vega,
            )
            if np.any(derivative_output.failed_contracts):
                raise FloatingPointError("semi-analytical IV Jacobian failed for one or more contracts")
            batch_iv = derivative_output.model_iv
            batch_jacobian = derivative_output.initial_variance_jacobian
        elif jacobian_method == "finite_difference":
            derivative_output = finite_difference_iv_variance_jacobian(
                log_spot=context.date.log_spot,
                implied_variance=variance,
                strikes=maturity_batch.strikes,
                option_types=maturity_batch.option_types,
                rate=maturity_batch.rate,
                dividend_yield=maturity_batch.dividend_yield,
                fixed_cos_basis=maturity_batch.basis,
                variance_bounds=(config.v_min, config.v_max),
                relative_step=config.finite_difference_relative_step,
                absolute_minimum_step=config.finite_difference_absolute_step,
            )
            batch_iv = derivative_output.model_iv
            batch_jacobian = derivative_output.initial_variance_jacobian
        else:
            raise ValueError(f"unsupported Jacobian method: {jacobian_method}")
        model_iv[maturity_batch.indices] = batch_iv
        initial_variance_jacobian[maturity_batch.indices] = batch_jacobian
    if not np.all(np.isfinite(model_iv)) or not np.all(np.isfinite(initial_variance_jacobian)):
        raise FloatingPointError("date-level IV or initial-variance Jacobian contains non-finite entries")
    return model_iv, initial_variance_jacobian


def _implied_variance_fit_loss(
    *,
    variance: float,
    context: _DateImpliedStatePricingContext,
    pricer: CosOptionPricer,
    config: ImpliedStateConfig,
) -> float:
    if variance < config.v_min or variance > config.v_max or not math.isfinite(variance):
        return float("inf")
    try:
        model_iv = _compute_model_implied_volatilities(variance=variance, context=context, pricer=pricer)
    except (ValueError, FloatingPointError, OverflowError):
        return float("inf")
    residual = context.date.observed_iv - model_iv
    finite = np.isfinite(residual)
    if not np.any(finite):
        return float("inf")
    return float(np.dot(residual[finite], residual[finite]))


def _local_variance_search_bounds(previous: float, config: ImpliedStateConfig) -> tuple[float, float]:
    if config.warm_start_window is None:
        return config.v_min, config.v_max
    lower = max(config.v_min, previous - config.warm_start_window)
    upper = min(config.v_max, previous + config.warm_start_window)
    return (config.v_min, config.v_max) if upper <= lower else (lower, upper)


def _imply_variance_for_date(
    *,
    context: _DateImpliedStatePricingContext,
    start_value: float,
    pricer: CosOptionPricer,
    config: ImpliedStateConfig,
) -> _ImpliedVarianceSolveResult:
    objective = lambda value: _implied_variance_fit_loss(
        variance=float(value), context=context, pricer=pricer, config=config
    )
    lower, upper = _local_variance_search_bounds(start_value, config)
    try:
        result = _solve_implied_variance_with_method(
            objective=objective,
            context=context,
            start_value=start_value,
            lower=lower,
            upper=upper,
            config=config,
            method=config.state_solver,
        )
    except (ValueError, FloatingPointError, OverflowError):
        result = _fallback_implied_variance_solve(objective=objective, config=config)
    hit_artificial_lower = lower > config.v_min and result.x <= lower + config.boundary_tol
    hit_artificial_upper = upper < config.v_max and result.x >= upper - config.boundary_tol
    if hit_artificial_lower or hit_artificial_upper:
        try:
            full = _solve_implied_variance_with_method(
                objective=objective,
                context=context,
                start_value=result.x,
                lower=config.v_min,
                upper=config.v_max,
                config=config,
                method=config.state_solver,
            )
        except (ValueError, FloatingPointError, OverflowError):
            full = _fallback_implied_variance_solve(objective=objective, config=config)
        return _ImpliedVarianceSolveResult(
            x=full.x,
            fun=full.fun,
            nfev=result.nfev + full.nfev,
            nit=result.nit + full.nit,
            success=full.success,
            njev=result.njev + full.njev,
            solver_name=full.solver_name,
            fallback_used=result.fallback_used or full.fallback_used,
        )
    return result


def _fallback_implied_variance_solve(*, objective, config: ImpliedStateConfig) -> _ImpliedVarianceSolveResult:
    if config.fallback_solver is None:
        raise FloatingPointError("selected implied-state solver failed and no fallback is configured")
    fallback_result = _solve_scalar_objective(
        objective,
        lower=config.v_min,
        upper=config.v_max,
        xatol=config.tol,
        max_iter=config.max_iter,
        method=config.fallback_solver,
    )
    return _ImpliedVarianceSolveResult(
        x=fallback_result.x,
        fun=fallback_result.fun,
        nfev=fallback_result.nfev,
        nit=fallback_result.nit,
        success=fallback_result.success,
        solver_name=fallback_result.solver_name,
        fallback_used=True,
    )


def _solve_implied_variance_with_method(
    *,
    objective,
    context: _DateImpliedStatePricingContext,
    start_value: float,
    lower: float,
    upper: float,
    config: ImpliedStateConfig,
    method: str,
) -> _ImpliedVarianceSolveResult:
    if method in {"golden_section", "bounded_brent"}:
        return _solve_scalar_objective(
            objective,
            lower=lower,
            upper=upper,
            xatol=config.tol,
            max_iter=config.max_iter,
            method=method,
        )

    from scipy.optimize import least_squares  # type: ignore[import-not-found]

    jacobian_method = "analytic" if method == "least_squares_analytic" else "finite_difference"
    residual_evaluation_count = 0
    jacobian_evaluation_count = 0

    def residual_function(candidate: np.ndarray) -> np.ndarray:
        nonlocal residual_evaluation_count
        residual_evaluation_count += 1
        model_iv = _compute_model_implied_volatilities(
            variance=float(candidate[0]), context=context, pricer=CosOptionPricer()
        )
        if not np.all(np.isfinite(model_iv)):
            raise FloatingPointError("least-squares residual contains non-finite model IVs")
        return context.date.observed_iv - model_iv

    def jacobian_function(candidate: np.ndarray) -> np.ndarray:
        nonlocal jacobian_evaluation_count
        jacobian_evaluation_count += 1
        _, model_iv_jacobian = _compute_model_iv_and_initial_variance_jacobian(
            variance=float(candidate[0]),
            context=context,
            config=config,
            jacobian_method=jacobian_method,
        )
        return -model_iv_jacobian[:, None]

    initial_value = min(max(start_value, lower), upper)
    least_squares_result = least_squares(
        residual_function,
        x0=np.array([initial_value]),
        jac=jacobian_function,
        bounds=(np.array([lower]), np.array([upper])),
        xtol=config.tol,
        ftol=config.tol,
        gtol=config.tol,
        max_nfev=config.max_iter,
    )
    loss = float(2.0 * least_squares_result.cost)
    return _ImpliedVarianceSolveResult(
        x=float(least_squares_result.x[0]),
        fun=loss,
        nfev=residual_evaluation_count,
        nit=int(least_squares_result.nfev),
        success=bool(least_squares_result.success and np.isfinite(loss)),
        njev=jacobian_evaluation_count,
        solver_name=method,
    )


def _solve_scalar_objective(
    objective,
    *,
    lower: float,
    upper: float,
    xatol: float,
    max_iter: int,
    method: str,
) -> _ImpliedVarianceSolveResult:
    if method == "golden_section":
        return _bounded_golden_section_search(
            objective, lower=lower, upper=upper, xatol=xatol, max_iter=max_iter
        )
    if method != "bounded_brent":
        raise ValueError(f"unsupported scalar implied-state solver: {method}")
    from scipy.optimize import minimize_scalar  # type: ignore[import-not-found]

    scipy_result = minimize_scalar(
        objective,
        method="bounded",
        bounds=(lower, upper),
        options={"xatol": xatol, "maxiter": max_iter},
    )
    return _ImpliedVarianceSolveResult(
        x=float(scipy_result.x),
        fun=float(scipy_result.fun),
        nfev=int(scipy_result.nfev),
        nit=int(scipy_result.nit),
        success=bool(scipy_result.success and np.isfinite(scipy_result.fun)),
        solver_name="bounded_brent",
    )


def _bounded_golden_section_search(
    func,
    *,
    lower: float,
    upper: float,
    xatol: float,
    max_iter: int,
) -> _ImpliedVarianceSolveResult:
    if upper <= lower:
        value = float(func(lower))
        return _ImpliedVarianceSolveResult(lower, value, 1, 0, math.isfinite(value))
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    inv_phi_sq = (3.0 - math.sqrt(5.0)) / 2.0
    a, b = float(lower), float(upper)
    h = b - a
    c, d = a + inv_phi_sq * h, a + inv_phi * h
    fc, fd = float(func(c)), float(func(d))
    nfev, nit = 2, 0
    for nit in range(1, max_iter + 1):
        if abs(b - a) <= xatol:
            break
        if fc < fd:
            b, d, fd = d, c, fc
            h = inv_phi * h
            c = a + inv_phi_sq * h
            fc = float(func(c))
        else:
            a, c, fc = c, d, fd
            h = inv_phi * h
            d = a + inv_phi * h
            fd = float(func(d))
        nfev += 1
    x, value = (c, fc) if fc < fd else (d, fd)
    return _ImpliedVarianceSolveResult(float(x), float(value), nfev, nit, math.isfinite(value))
