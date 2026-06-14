from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ImpliedVolatility.black_iv import implied_vol_black76
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver
from OptionPricing.types import CosPricingConfig

from Estimation.ISCGMM.types import EstimationPanel, HestonEstimationParams, ImpliedStateConfig, ImpliedStateResult, PanelDate


@dataclass(frozen=True)
class _ScalarMinResult:
    x: float
    fun: float
    nfev: int
    nit: int
    success: bool


@dataclass(frozen=True)
class _ContractGroup:
    maturity: float
    indices: np.ndarray
    rate: float
    dividend_yield: float


@dataclass(frozen=True)
class _DatePricingCache:
    date: PanelDate
    groups: tuple[_ContractGroup, ...]


def _bounded_minimize_scalar(
    func,
    *,
    lower: float,
    upper: float,
    xatol: float,
    max_iter: int,
) -> _ScalarMinResult:
    if upper <= lower:
        value = float(func(lower))
        return _ScalarMinResult(x=lower, fun=value, nfev=1, nit=0, success=math.isfinite(value))

    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    inv_phi_sq = (3.0 - math.sqrt(5.0)) / 2.0

    a = float(lower)
    b = float(upper)
    h = b - a
    c = a + inv_phi_sq * h
    d = a + inv_phi * h
    fc = float(func(c))
    fd = float(func(d))
    nfev = 2
    nit = 0

    for nit in range(1, max_iter + 1):
        if abs(b - a) <= xatol:
            break
        if fc < fd:
            b = d
            d = c
            fd = fc
            h = inv_phi * h
            c = a + inv_phi_sq * h
            fc = float(func(c))
        else:
            a = c
            c = d
            fc = fd
            h = inv_phi * h
            d = a + inv_phi * h
            fd = float(func(d))
        nfev += 1

    if fc < fd:
        x = c
        fun = fc
    else:
        x = d
        fun = fd
    return _ScalarMinResult(x=float(x), fun=float(fun), nfev=nfev, nit=nit, success=math.isfinite(fun))


def _make_date_cache(date: PanelDate) -> _DatePricingCache:
    groups: list[_ContractGroup] = []
    for maturity in np.unique(date.maturities):
        idx = np.flatnonzero(np.isclose(date.maturities, maturity, rtol=0.0, atol=1e-14))
        rates = date.rates[idx]
        dividend_yields = date.dividend_yields[idx]
        if not np.allclose(rates, rates[0]):
            raise ValueError("rates must be constant within each maturity/date group")
        if not np.allclose(dividend_yields, dividend_yields[0]):
            raise ValueError("dividend yields must be constant within each maturity/date group")
        groups.append(
            _ContractGroup(
                maturity=float(maturity),
                indices=idx,
                rate=float(rates[0]),
                dividend_yield=float(dividend_yields[0]),
            )
        )
    return _DatePricingCache(date=date, groups=tuple(groups))


def _pricing_config_for_group(
    *,
    variance: float,
    maturity: float,
    base_config: ImpliedStateConfig,
) -> tuple[CosPricingConfig, float]:
    if base_config.effective_truncation_width is None:
        maturity_grid = np.array([maturity], dtype=float)
        width = CosOptionPricer.effective_truncation_width(
            np.array([variance], dtype=float),
            maturity_grid,
            base_config.pricing_config,
        )
        return base_config.pricing_config, width

    target_width = float(base_config.effective_truncation_width)
    vol_time_scale = math.sqrt(max(float(variance), 1e-12) * max(float(maturity), 1e-12))
    local_config = CosPricingConfig(
        n_cos=base_config.pricing_config.n_cos,
        truncation_width=target_width / vol_time_scale,
    )
    return local_config, target_width


def _model_iv_for_date(
    *,
    variance: float,
    cache: _DatePricingCache,
    theta: HestonEstimationParams,
    pricer: CosOptionPricer,
    solver: HestonAnalyticCcfSolver,
    config: ImpliedStateConfig,
) -> np.ndarray:
    date = cache.date
    model_iv = np.empty(date.n_contracts, dtype=float)
    params_q = theta.to_pricing_q()
    for group in cache.groups:
        idx = group.indices
        maturity_grid = np.array([group.maturity], dtype=float)
        pricing_config, truncation_width = _pricing_config_for_group(
            variance=variance,
            maturity=group.maturity,
            base_config=config,
        )
        u_grid, _ = pricer._get_static_terms(pricing_config.n_cos, truncation_width)
        coefficients = solver.solve_coefficients(
            u_grid,
            maturity_grid,
            model_params=params_q,
        )
        option_types = date.option_types[idx].reshape(1, -1, 1)
        prices = pricer.price_matrix(
            log_s=np.array([date.log_spot], dtype=float),
            variance=np.array([variance], dtype=float),
            strike_grid=date.strikes[idx],
            maturity_grid=maturity_grid,
            rate_grid=np.array([group.rate], dtype=float),
            dividend_yield_grid=np.array([group.dividend_yield], dtype=float),
            coefficients=coefficients,
            pricing_config=pricing_config,
            option_type=option_types,
        )[0, :, 0]
        for local_pos, price in zip(idx, prices):
            tau = float(date.maturities[local_pos])
            rate = float(date.rates[local_pos])
            q = float(date.dividend_yields[local_pos])
            forward = date.spot * math.exp((rate - q) * tau)
            discount = math.exp(-rate * tau)
            try:
                model_iv[local_pos] = implied_vol_black76(
                    price=float(price),
                    forward=forward,
                    strike=float(date.strikes[local_pos]),
                    tau=tau,
                    discount_factor=discount,
                    option_type=str(date.option_types[local_pos]),
                    on_bounds="clip",
                )
            except (ValueError, FloatingPointError):
                model_iv[local_pos] = np.nan
    return model_iv


def _date_objective(
    *,
    variance: float,
    cache: _DatePricingCache,
    theta: HestonEstimationParams,
    pricer: CosOptionPricer,
    solver: HestonAnalyticCcfSolver,
    config: ImpliedStateConfig,
) -> float:
    if variance < config.v_min or variance > config.v_max or not math.isfinite(variance):
        return float("inf")
    try:
        model_iv = _model_iv_for_date(
            variance=variance,
            cache=cache,
            theta=theta,
            pricer=pricer,
            solver=solver,
            config=config,
        )
    except (ValueError, FloatingPointError, OverflowError):
        return float("inf")
    residual = cache.date.observed_iv - model_iv
    finite = np.isfinite(residual)
    if not np.any(finite):
        return float("inf")
    return float(np.dot(residual[finite], residual[finite]))


def _warm_bounds(previous: float, config: ImpliedStateConfig) -> tuple[float, float]:
    if config.warm_start_window is None:
        return config.v_min, config.v_max
    lower = max(config.v_min, previous - config.warm_start_window)
    upper = min(config.v_max, previous + config.warm_start_window)
    if upper <= lower:
        return config.v_min, config.v_max
    return lower, upper


def _minimize_date_variance(
    *,
    cache: _DatePricingCache,
    theta: HestonEstimationParams,
    start_value: float,
    pricer: CosOptionPricer,
    solver: HestonAnalyticCcfSolver,
    config: ImpliedStateConfig,
) -> _ScalarMinResult:
    objective = lambda value: _date_objective(
        variance=float(value),
        cache=cache,
        theta=theta,
        pricer=pricer,
        solver=solver,
        config=config,
    )
    lower, upper = _warm_bounds(start_value, config)
    result = _bounded_minimize_scalar(
        objective,
        lower=lower,
        upper=upper,
        xatol=config.tol,
        max_iter=config.max_iter,
    )
    interior_lower = lower > config.v_min and result.x <= lower + config.boundary_tol
    interior_upper = upper < config.v_max and result.x >= upper - config.boundary_tol
    if interior_lower or interior_upper:
        global_result = _bounded_minimize_scalar(
            objective,
            lower=config.v_min,
            upper=config.v_max,
            xatol=config.tol,
            max_iter=config.max_iter,
        )
        return _ScalarMinResult(
            x=global_result.x,
            fun=global_result.fun,
            nfev=result.nfev + global_result.nfev,
            nit=result.nit + global_result.nit,
            success=global_result.success,
        )
    return result


def imply_heston_variance_path(
    theta: HestonEstimationParams,
    panel: EstimationPanel,
    config: ImpliedStateConfig | None = None,
) -> ImpliedStateResult:
    """Imply the scalar Heston variance path from date-by-date IV panels."""

    cfg = ImpliedStateConfig() if config is None else config
    cfg.validate()
    theta.validate()
    caches = tuple(_make_date_cache(date) for date in panel.dates)
    pricer = CosOptionPricer()
    solver = HestonAnalyticCcfSolver()

    n_dates = panel.n_dates
    variance = np.empty(n_dates, dtype=float)
    objective = np.empty(n_dates, dtype=float)
    boundary_hit = np.zeros(n_dates, dtype=bool)
    failed = np.zeros(n_dates, dtype=bool)
    nfev = np.zeros(n_dates, dtype=int)
    start_values = np.empty(n_dates, dtype=float)

    previous = min(max(theta.vbar, cfg.v_min), cfg.v_max)
    for date_idx, cache in enumerate(caches):
        start_values[date_idx] = previous
        result = _minimize_date_variance(
            cache=cache,
            theta=theta,
            start_value=previous,
            pricer=pricer,
            solver=solver,
            config=cfg,
        )
        variance[date_idx] = min(max(result.x, cfg.v_min), cfg.v_max)
        objective[date_idx] = result.fun
        nfev[date_idx] = result.nfev
        failed[date_idx] = not result.success
        boundary_hit[date_idx] = (
            variance[date_idx] <= cfg.v_min + cfg.boundary_tol
            or variance[date_idx] >= cfg.v_max - cfg.boundary_tol
        )
        previous = variance[date_idx]

    return ImpliedStateResult(
        variance=variance,
        objective=objective,
        boundary_hit=boundary_hit,
        failed=failed,
        nfev=nfev,
        start_values=start_values,
    )
