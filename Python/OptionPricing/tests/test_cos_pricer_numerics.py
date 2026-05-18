import math

import numpy as np
import pytest

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_price
from OptionPricing.tests.pricing_test_helpers import cos_heston_price, forward_strike, mc_heston_price
from OptionPricing.types import CosPricingConfig, HestonPricingParamsQ


def _baseline_params() -> HestonPricingParamsQ:
    return HestonPricingParamsQ(kappa=3.0, vbar=0.04, sigma_v=0.4, rho=-0.7, r=0.02, q=0.0)


def test_cos_put_call_parity():
    params = _baseline_params()
    config = CosPricingConfig(n_cos=256, truncation_width=10.0)
    spot = 100.0
    variance = 0.04
    tau = 60.0 / 252.0
    strike = 103.0

    call = cos_heston_price(
        spot=spot,
        variance=variance,
        tau=tau,
        strike=strike,
        option_type="call",
        params=params,
        config=config,
    )
    put = cos_heston_price(
        spot=spot,
        variance=variance,
        tau=tau,
        strike=strike,
        option_type="put",
        params=params,
        config=config,
    )
    forward = spot * math.exp((params.r - params.q) * tau)
    parity = math.exp(-params.r * tau) * (forward - strike)

    assert call - put == pytest.approx(parity, abs=2e-10)


def test_cos_otm_prices_are_finite_and_bounded():
    params = _baseline_params()
    config = CosPricingConfig(n_cos=256, truncation_width=10.0)
    spot = 100.0
    variance = 0.04
    tau = 60.0 / 252.0
    discount = math.exp(-params.r * tau)
    forward = spot * math.exp((params.r - params.q) * tau)

    cases = [
        ("put", forward_strike(spot=spot, tau=tau, log_moneyness=-0.1, r=params.r, q=params.q), discount * forward),
        ("call", forward_strike(spot=spot, tau=tau, log_moneyness=0.1, r=params.r, q=params.q), discount * forward),
    ]
    for option_type, strike, upper_bound in cases:
        price = cos_heston_price(
            spot=spot,
            variance=variance,
            tau=tau,
            strike=strike,
            option_type=option_type,
            params=params,
            config=config,
        )
        assert math.isfinite(price)
        assert 0.0 <= price <= upper_bound


def test_cos_price_stabilizes_as_terms_increase():
    params = _baseline_params()
    spot = 100.0
    variance = 0.04
    tau = 60.0 / 252.0
    strike = forward_strike(spot=spot, tau=tau, log_moneyness=0.075, r=params.r, q=params.q)

    p128 = cos_heston_price(
        spot=spot,
        variance=variance,
        tau=tau,
        strike=strike,
        option_type="call",
        params=params,
        config=CosPricingConfig(n_cos=128, truncation_width=10.0),
    )
    p256 = cos_heston_price(
        spot=spot,
        variance=variance,
        tau=tau,
        strike=strike,
        option_type="call",
        params=params,
        config=CosPricingConfig(n_cos=256, truncation_width=10.0),
    )

    assert p256 == pytest.approx(p128, abs=2e-3)


def test_black_scholes_pricing_and_iv_edges():
    cases = [
        ("call", 100.0, 100.0, 1.0, 0.2),
        ("call", 100.0, 130.0, 0.25, 0.35),
        ("put", 100.0, 70.0, 0.25, 0.35),
        ("call", 100.0, 105.0, 5.0 / 252.0, 0.08),
        ("put", 100.0, 99.0, 2.0, 0.75),
    ]
    discount = 0.997
    for option_type, forward, strike, tau, vol in cases:
        price = black76_price(
            forward=forward,
            strike=strike,
            tau=tau,
            vol=vol,
            discount_factor=discount,
            option_type=option_type,
        )
        implied = implied_vol_black76(
            price=price,
            forward=forward,
            strike=strike,
            tau=tau,
            discount_factor=discount,
            option_type=option_type,
        )
        assert implied == pytest.approx(vol, abs=1e-11)


@pytest.mark.slow
@pytest.mark.parametrize(
    ("n_weeks", "log_moneyness", "option_type"),
    [
        (12, 0.1, "call"),
        (12, -0.1, "put"),
        (26, 0.1, "call"),
        (26, -0.1, "put"),
    ],
)
def test_cos_matches_q_measure_monte_carlo_for_otm_options(n_weeks, log_moneyness, option_type):
    params = _baseline_params()
    config = CosPricingConfig(n_cos=256, truncation_width=10.0)
    spot = 100.0
    variance = 0.04
    tau = n_weeks * 5.0 / 252.0
    strike = forward_strike(spot=spot, tau=tau, log_moneyness=log_moneyness, r=params.r, q=params.q)

    cos_price = cos_heston_price(
        spot=spot,
        variance=variance,
        tau=tau,
        strike=strike,
        option_type=option_type,
        params=params,
        config=config,
    )
    estimate = mc_heston_price(
        spot=spot,
        variance=variance,
        n_weeks=n_weeks,
        strike=strike,
        option_type=option_type,
        params=params,
        n_paths=1500,
    )

    assert abs(cos_price - estimate.mean) <= max(0.15, 4.0 * estimate.standard_error)


@pytest.mark.slow
@pytest.mark.parametrize(
    "params",
    [
        HestonPricingParamsQ(kappa=3.0, vbar=0.04, sigma_v=0.05, rho=-0.2, r=0.02, q=0.0),
        HestonPricingParamsQ(kappa=1.2, vbar=0.04, sigma_v=0.9, rho=-0.85, r=0.02, q=0.0),
        HestonPricingParamsQ(kappa=2.5, vbar=0.04, sigma_v=0.4, rho=0.0, r=0.02, q=0.0),
    ],
)
def test_cos_monte_carlo_stress_cases_are_plausible(params):
    config = CosPricingConfig(n_cos=256, truncation_width=10.0)
    spot = 100.0
    variance = 0.04
    n_weeks = 12
    tau = n_weeks * 5.0 / 252.0
    strike = forward_strike(spot=spot, tau=tau, log_moneyness=0.1, r=params.r, q=params.q)

    cos_price = cos_heston_price(
        spot=spot,
        variance=variance,
        tau=tau,
        strike=strike,
        option_type="call",
        params=params,
        config=config,
    )
    estimate = mc_heston_price(
        spot=spot,
        variance=variance,
        n_weeks=n_weeks,
        strike=strike,
        option_type="call",
        params=params,
        n_paths=1200,
        base_seed=51_000,
    )

    assert math.isfinite(cos_price)
    assert abs(cos_price - estimate.mean) <= max(0.25, 5.0 * estimate.standard_error)
