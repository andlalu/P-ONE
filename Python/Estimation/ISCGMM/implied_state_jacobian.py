from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_vega
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.types import PreparedFixedCosBasis


@dataclass(frozen=True)
class FiniteDifferenceJacobianDiagnostics:
    step_used: float
    lower_evaluation_point: float
    upper_evaluation_point: float
    stencil_type: str
    failed_contracts: np.ndarray
    low_vega_contracts: np.ndarray
    function_evaluation_count: int


@dataclass(frozen=True)
class FiniteDifferenceIvVarianceJacobian:
    model_iv: np.ndarray
    initial_variance_jacobian: np.ndarray
    diagnostics: FiniteDifferenceJacobianDiagnostics


@dataclass(frozen=True)
class PriceIvVarianceJacobian:
    prices: np.ndarray
    model_iv: np.ndarray
    price_jacobian: np.ndarray
    initial_variance_jacobian: np.ndarray
    vegas: np.ndarray
    failed_contracts: np.ndarray
    low_vega_contracts: np.ndarray
    price_or_iv_boundary_contracts: np.ndarray
    pricing_call_count: int


def _model_iv_fixed_basis(
    *,
    log_spot: float,
    implied_variance: float,
    strikes: np.ndarray,
    option_types: np.ndarray,
    rate: float,
    dividend_yield: float,
    fixed_cos_basis: PreparedFixedCosBasis,
) -> np.ndarray:
    pricer = CosOptionPricer()
    strike_array = np.asarray(strikes, dtype=float)
    option_type_array = np.asarray(option_types).astype(str)
    prices = pricer.price_matrix_fixed_basis(
        log_s=np.array([log_spot]),
        variance=np.array([implied_variance]),
        strike_grid=strike_array,
        rate=rate,
        dividend_yield=dividend_yield,
        basis=fixed_cos_basis,
        option_type=option_type_array.reshape(1, -1, 1),
    )[0, :, 0]
    maturity = fixed_cos_basis.maturity
    spot = float(np.exp(log_spot))
    forward = spot * float(np.exp((rate - dividend_yield) * maturity))
    discount = float(np.exp(-rate * maturity))
    model_iv = np.empty(strike_array.size, dtype=float)
    for contract_index, (price, strike, option_type) in enumerate(
        zip(prices, strike_array, option_type_array)
    ):
        model_iv[contract_index] = implied_vol_black76(
            price=float(price),
            forward=forward,
            strike=float(strike),
            tau=maturity,
            discount_factor=discount,
            option_type=str(option_type),
            on_bounds="raise",
        )
    return model_iv


def finite_difference_iv_variance_jacobian(
    *,
    log_spot: float,
    implied_variance: float,
    strikes: np.ndarray,
    option_types: np.ndarray,
    rate: float,
    dividend_yield: float,
    fixed_cos_basis: PreparedFixedCosBasis,
    variance_bounds: tuple[float, float],
    relative_step: float = 1e-4,
    absolute_minimum_step: float = 1e-7,
    minimum_black_vega: float = 1e-8,
) -> FiniteDifferenceIvVarianceJacobian:
    """Finite-difference ``d model IV / d initial variance`` for one maturity.

    Central differences are used in the interior and deterministic one-sided
    differences are used when the requested perturbation reaches a state bound.
    """

    variance_lower, variance_upper = (float(bound) for bound in variance_bounds)
    if variance_lower < 0.0 or variance_upper <= variance_lower:
        raise ValueError("variance_bounds must satisfy 0 <= lower < upper")
    if not variance_lower <= implied_variance <= variance_upper:
        raise ValueError("implied_variance lies outside variance_bounds")
    if relative_step <= 0.0 or absolute_minimum_step <= 0.0:
        raise ValueError("finite-difference steps must be strictly positive")
    if minimum_black_vega <= 0.0:
        raise ValueError("minimum_black_vega must be strictly positive")
    requested_step = max(abs(implied_variance) * relative_step, absolute_minimum_step)
    lower_point = max(variance_lower, implied_variance - requested_step)
    upper_point = min(variance_upper, implied_variance + requested_step)
    if lower_point == upper_point:
        raise FloatingPointError("variance bounds leave no usable finite-difference perturbation")

    base_iv = _model_iv_fixed_basis(
        log_spot=log_spot,
        implied_variance=implied_variance,
        strikes=strikes,
        option_types=option_types,
        rate=rate,
        dividend_yield=dividend_yield,
        fixed_cos_basis=fixed_cos_basis,
    )
    if lower_point < implied_variance and upper_point > implied_variance:
        lower_iv = _model_iv_fixed_basis(
            log_spot=log_spot,
            implied_variance=lower_point,
            strikes=strikes,
            option_types=option_types,
            rate=rate,
            dividend_yield=dividend_yield,
            fixed_cos_basis=fixed_cos_basis,
        )
        upper_iv = _model_iv_fixed_basis(
            log_spot=log_spot,
            implied_variance=upper_point,
            strikes=strikes,
            option_types=option_types,
            rate=rate,
            dividend_yield=dividend_yield,
            fixed_cos_basis=fixed_cos_basis,
        )
        jacobian = (upper_iv - lower_iv) / (upper_point - lower_point)
        stencil = "central"
        function_evaluation_count = 3
    elif upper_point > implied_variance:
        upper_iv = _model_iv_fixed_basis(
            log_spot=log_spot,
            implied_variance=upper_point,
            strikes=strikes,
            option_types=option_types,
            rate=rate,
            dividend_yield=dividend_yield,
            fixed_cos_basis=fixed_cos_basis,
        )
        jacobian = (upper_iv - base_iv) / (upper_point - implied_variance)
        stencil = "forward"
        function_evaluation_count = 2
    else:
        lower_iv = _model_iv_fixed_basis(
            log_spot=log_spot,
            implied_variance=lower_point,
            strikes=strikes,
            option_types=option_types,
            rate=rate,
            dividend_yield=dividend_yield,
            fixed_cos_basis=fixed_cos_basis,
        )
        jacobian = (base_iv - lower_iv) / (implied_variance - lower_point)
        stencil = "backward"
        function_evaluation_count = 2

    failed_contracts = ~np.isfinite(jacobian)
    if np.any(failed_contracts):
        raise FloatingPointError("finite-difference IV Jacobian contains failed contracts")
    maturity = fixed_cos_basis.maturity
    forward = float(np.exp(log_spot + (rate - dividend_yield) * maturity))
    discount = float(np.exp(-rate * maturity))
    vegas = np.array(
        [
            black76_vega(
                forward=forward,
                strike=float(strike),
                tau=maturity,
                vol=float(volatility),
                discount_factor=discount,
            )
            for strike, volatility in zip(np.asarray(strikes, dtype=float), base_iv)
        ]
    )
    low_vega_contracts = vegas < minimum_black_vega
    return FiniteDifferenceIvVarianceJacobian(
        model_iv=base_iv,
        initial_variance_jacobian=jacobian,
        diagnostics=FiniteDifferenceJacobianDiagnostics(
            step_used=max(implied_variance - lower_point, upper_point - implied_variance),
            lower_evaluation_point=lower_point,
            upper_evaluation_point=upper_point,
            stencil_type=stencil,
            failed_contracts=failed_contracts,
            low_vega_contracts=low_vega_contracts,
            function_evaluation_count=function_evaluation_count,
        ),
    )


def price_iv_and_initial_variance_jacobian_fixed_basis(
    *,
    log_spot: float,
    implied_variance: float,
    strikes: np.ndarray,
    option_types: np.ndarray,
    rate: float,
    dividend_yield: float,
    fixed_cos_basis: PreparedFixedCosBasis,
    minimum_black_vega: float = 1e-8,
) -> PriceIvVarianceJacobian:
    """Evaluate fixed-basis price, IV and their initial-variance Jacobians."""

    if minimum_black_vega <= 0.0:
        raise ValueError("minimum_black_vega must be strictly positive")
    strike_array = np.asarray(strikes, dtype=float)
    option_type_array = np.char.lower(np.asarray(option_types).astype(str))
    price_output = CosOptionPricer().price_and_initial_variance_jacobian_fixed_basis(
        log_s=log_spot,
        variance=implied_variance,
        strikes=strike_array,
        rate=rate,
        dividend_yield=dividend_yield,
        basis=fixed_cos_basis,
        option_types=option_type_array,
    )
    maturity = fixed_cos_basis.maturity
    spot = float(np.exp(log_spot))
    forward = spot * float(np.exp((rate - dividend_yield) * maturity))
    discount = float(np.exp(-rate * maturity))
    model_iv = np.full(strike_array.size, np.nan)
    vegas = np.full(strike_array.size, np.nan)
    iv_jacobian = np.full(strike_array.size, np.nan)
    low_vega = np.zeros(strike_array.size, dtype=bool)
    boundary = price_output.price_clipped.copy()
    failed = np.zeros(strike_array.size, dtype=bool)

    for contract_index, (price, price_jacobian, strike, option_type) in enumerate(
        zip(price_output.prices, price_output.initial_variance_jacobian, strike_array, option_type_array)
    ):
        intrinsic = discount * max((forward - strike) if option_type == "call" else (strike - forward), 0.0)
        upper_bound = discount * (forward if option_type == "call" else strike)
        bound_epsilon = 1e-12 * max(1.0, upper_bound)
        if price <= intrinsic + bound_epsilon or price >= upper_bound - bound_epsilon:
            boundary[contract_index] = True
        try:
            model_iv[contract_index] = implied_vol_black76(
                price=float(price),
                forward=forward,
                strike=float(strike),
                tau=maturity,
                discount_factor=discount,
                option_type=str(option_type),
                on_bounds="raise",
            )
        except (ValueError, FloatingPointError):
            failed[contract_index] = True
            continue
        vegas[contract_index] = black76_vega(
            forward=forward,
            strike=float(strike),
            tau=maturity,
            vol=float(model_iv[contract_index]),
            discount_factor=discount,
        )
        low_vega[contract_index] = vegas[contract_index] < minimum_black_vega
        if boundary[contract_index] or low_vega[contract_index] or not np.isfinite(price_jacobian):
            failed[contract_index] = True
            continue
        iv_jacobian[contract_index] = price_jacobian / vegas[contract_index]

    return PriceIvVarianceJacobian(
        prices=price_output.prices,
        model_iv=model_iv,
        price_jacobian=price_output.initial_variance_jacobian,
        initial_variance_jacobian=iv_jacobian,
        vegas=vegas,
        failed_contracts=failed,
        low_vega_contracts=low_vega,
        price_or_iv_boundary_contracts=boundary,
        pricing_call_count=1,
    )


def implied_variance_loss_and_jacobian(
    observed_iv: np.ndarray,
    model_iv: np.ndarray,
    initial_variance_jacobian: np.ndarray,
) -> tuple[float, float]:
    """Return ``L(v)`` and ``dL/dv = -2 J(v)' residual(v)``."""

    observed = np.asarray(observed_iv, dtype=float)
    model = np.asarray(model_iv, dtype=float)
    jacobian = np.asarray(initial_variance_jacobian, dtype=float)
    if observed.ndim != 1 or observed.shape != model.shape or observed.shape != jacobian.shape:
        raise ValueError("observed IV, model IV and Jacobian must be aligned 1D arrays")
    if not np.all(np.isfinite(observed)) or not np.all(np.isfinite(model)) or not np.all(np.isfinite(jacobian)):
        raise FloatingPointError("loss derivative inputs must all be finite")
    residual = observed - model
    return float(residual @ residual), float(-2.0 * jacobian @ residual)
