from __future__ import annotations

import math

import numpy as np

from Models.Heston.parameters import HestonRiskNeutralParameters


def _heston_log_spot_characteristic_function(
    transform_argument: complex,
    *,
    log_spot: float,
    initial_variance: float,
    maturity: float,
    parameters: HestonRiskNeutralParameters,
) -> complex:
    """Heston characteristic function of terminal log spot for complex arguments."""

    u = complex(transform_argument)
    iu = 1j * u
    kappa = parameters.kappa
    sigma = parameters.sigma_v
    rho = parameters.rho
    sigma_squared = sigma * sigma
    discriminant = np.sqrt((rho * sigma * iu - kappa) ** 2 + sigma_squared * (iu + u * u))
    numerator = kappa - rho * sigma * iu - discriminant
    denominator = kappa - rho * sigma * iu + discriminant
    ratio = numerator / denominator
    exponential = np.exp(-discriminant * maturity)
    affine_a = (kappa * parameters.vbar / sigma_squared) * (
        numerator * maturity - 2.0 * np.log((1.0 - ratio * exponential) / (1.0 - ratio))
    )
    affine_b = (numerator / sigma_squared) * ((1.0 - exponential) / (1.0 - ratio * exponential))
    drifted_log_spot = log_spot + (parameters.r - parameters.q) * maturity
    return complex(np.exp(iu * drifted_log_spot + affine_a + affine_b * initial_variance))


def heston_option_price_fourier_reference(
    *,
    spot: float,
    initial_variance: float,
    maturity: float,
    strike: float,
    option_type: str,
    parameters: HestonRiskNeutralParameters,
    integration_limit: float = 200.0,
    absolute_tolerance: float = 2e-10,
    relative_tolerance: float = 2e-10,
    integration_subinterval_limit: int = 400,
) -> float:
    """Independent Heston price from Carr--Madan damped Fourier inversion."""

    from scipy.integrate import quad  # type: ignore[import-not-found]

    parameters.validate()
    if spot <= 0.0 or strike <= 0.0 or initial_variance < 0.0 or maturity <= 0.0:
        raise ValueError("spot, strike and maturity must be positive and initial_variance non-negative")
    if integration_limit <= 0.0 or absolute_tolerance <= 0.0 or relative_tolerance <= 0.0:
        raise ValueError("integration settings must be strictly positive")
    kind = option_type.lower()
    if kind not in {"call", "put"}:
        raise ValueError("option_type must be 'call' or 'put'")
    # Heston prices are homogeneous in spot and strike. Normalising spot to
    # one keeps the damped integrand near unit scale and materially improves
    # quadrature accuracy for very small short-maturity time values.
    normalised_strike = strike / spot
    log_spot = 0.0
    log_strike = math.log(normalised_strike)
    damping = 1.5

    def damped_call_integrand(frequency: float) -> float:
        characteristic_value = _heston_log_spot_characteristic_function(
            frequency - 1j * (damping + 1.0),
            log_spot=log_spot,
            initial_variance=initial_variance,
            maturity=maturity,
            parameters=parameters,
        )
        denominator = (
            damping * damping
            + damping
            - frequency * frequency
            + 1j * (2.0 * damping + 1.0) * frequency
        )
        return float(
            np.real(
                np.exp(-1j * frequency * log_strike)
                * math.exp(-parameters.r * maturity)
                * characteristic_value
                / denominator
            )
        )

    integral, _ = quad(
        damped_call_integrand,
        0.0,
        integration_limit,
        epsabs=absolute_tolerance,
        epsrel=relative_tolerance,
        limit=integration_subinterval_limit,
    )
    call_price = spot * math.exp(-damping * log_strike) * integral / math.pi
    discounted_spot = spot * math.exp(-parameters.q * maturity)
    discounted_strike = strike * math.exp(-parameters.r * maturity)
    if kind == "call":
        return float(call_price)
    return float(call_price - discounted_spot + discounted_strike)
