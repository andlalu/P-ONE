import math

import numpy as np
import pytest

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_price, black76_vega
from OptionData.add_noise import generate_noisy_panel_rows
from OptionData.noise_common import NoiseSettings, apply_price_mechanics, marginal_scale
from OptionData.noise_persistent_factor import (
    compute_q_diag_from_stationary_std,
    persistent_factor_noise,
    persistent_factor_residual_scale,
    validate_persistent_factor_settings,
)


def _noise_settings() -> NoiseSettings:
    return NoiseSettings(
        base_seed=9900000,
        sigma_min=0.0001,
        price_epsilon=1e-10,
        tick_size=0.01,
        scenarios={
            "low_iid": {
                "alpha_0": 0.0025,
                "alpha_m": 1.0,
                "alpha_tau": 0.05,
                "tau_min": 0.01984126984126984,
            },
            "spatial_corr": {
                "alpha_0": 0.0050,
                "alpha_m": 1.0,
                "alpha_tau": 0.05,
                "tau_min": 0.01984126984126984,
                "ell_m": 0.10,
                "ell_tau": 0.25,
                "correlation_jitter": 1e-10,
                "max_correlation_jitter": 1e-6,
            },
            "persistent_factor": {
                "alpha_0": 0.0015,
                "alpha_m": 2.0,
                "alpha_tau": 0.05,
                "tau_min": 0.01984126984126984,
                "a_diag": [0.85, 0.65, 0.65],
                "stationary_factor_std": [0.0007, 0.0025, 0.0008],
                "q_diag": [1.35975e-7, 3.609375e-6, 3.696e-7],
                "residual_policy": "match_total_marginal_scale",
                "factor_initialization": "zero",
            },
        },
    )


def _clean_rows(sample_id=0):
    rows = []
    for week_index, spot in enumerate((100.0, 101.0)):
        for tau in (0.25, 0.5):
            forward = spot * math.exp(0.02 * tau)
            discount = math.exp(-0.02 * tau)
            for lm in (-0.05, 0.0, 0.05):
                strike = forward * math.exp(lm)
                option_type = "put" if lm < 0.0 else "call"
                vol = 0.20 + 0.02 * tau + 0.01 * abs(lm)
                price = black76_price(
                    forward=forward,
                    strike=strike,
                    tau=tau,
                    vol=vol,
                    discount_factor=discount,
                    option_type=option_type,
                )
                rows.append(
                    {
                        "run_id": "test_run",
                        "sample_id": sample_id,
                        "week_index": week_index,
                        "t": week_index / 52.0,
                        "S": spot,
                        "logS": math.log(spot),
                        "V": 0.04,
                        "r": 0.02,
                        "q": 0.0,
                        "maturity_years": tau,
                        "expiry_time": week_index / 52.0 + tau,
                        "forward": forward,
                        "log_moneyness": lm,
                        "strike": strike,
                        "option_type": option_type,
                        "is_otm": True,
                        "pricing_method": "COS",
                        "model_price": price,
                        "model_iv": implied_vol_black76(
                            price=price,
                            forward=forward,
                            strike=strike,
                            tau=tau,
                            discount_factor=discount,
                            option_type=option_type,
                        ),
                        "model_vega": black76_vega(
                            forward=forward,
                            strike=strike,
                            tau=tau,
                            vol=vol,
                            discount_factor=discount,
                        ),
                        "iv_method": "lets_be_rational",
                    }
                )
    return rows


def test_low_iid_noise_is_deterministic():
    config = _noise_settings()
    rows = _clean_rows()
    first, _ = generate_noisy_panel_rows(rows, scenario="low_iid", seed=123, config=config)
    second, _ = generate_noisy_panel_rows(rows, scenario="low_iid", seed=123, config=config)
    assert [row["raw_noisy_iv"] for row in first] == [row["raw_noisy_iv"] for row in second]
    assert all(row["noise_scenario"] == "low_iid" for row in first)
    assert all(row["observed_iv"] >= config.sigma_min for row in first)


def test_spatial_corr_noise_has_contract_level_draws():
    config = _noise_settings()
    rows, _ = generate_noisy_panel_rows(_clean_rows(), scenario="spatial_corr", seed=456, config=config)
    draws = np.array([row["noise_draw"] for row in rows])
    assert np.std(draws) > 0.0
    assert all(row["noise_scenario"] == "spatial_corr" for row in rows)


def test_persistent_factor_q_diag_is_computed_from_stationary_std():
    a_diag = np.array([0.85, 0.65, 0.65])
    stationary_factor_std = np.array([0.0007, 0.0025, 0.0008])
    expected = np.array([1.35975e-7, 3.609375e-6, 3.696e-7])

    assert np.allclose(compute_q_diag_from_stationary_std(a_diag, stationary_factor_std), expected)

    config = _noise_settings().scenarios["persistent_factor"]
    assert np.allclose(np.array(config["q_diag"]), expected)


def test_persistent_factor_q_diag_consistency_validation_rejects_stationary_std_values():
    factor = dict(_noise_settings().scenarios["persistent_factor"])
    factor["q_diag"] = [0.0007, 0.0025, 0.0008]
    with pytest.raises(ValueError, match="innovation covariance diagonal"):
        validate_persistent_factor_settings(factor)


@pytest.mark.parametrize(
    ("bad_fields", "message"),
    [
        ({"a_diag": [0.85, 0.65]}, "a_diag must have length 3"),
        ({"a_diag": [1.0, 0.65, 0.65]}, "absolute value strictly below 1"),
        ({"stationary_factor_std": [0.0007, -0.0025, 0.0008]}, "stationary_factor_std entries must be non-negative"),
        ({"q_diag": [1.0e-7, -2.0e-7, 3.0e-7]}, "q_diag entries must be non-negative"),
    ],
)
def test_persistent_factor_config_validation_rejects_invalid_values(bad_fields, message):
    factor_config = dict(_noise_settings().scenarios["persistent_factor"])
    factor_config.update(bad_fields)

    with pytest.raises(ValueError, match=message):
        validate_persistent_factor_settings(factor_config)


def test_persistent_factor_innovation_draw_uses_sqrt_q_diag_scale():
    config = _noise_settings()

    class RecordingRng:
        def __init__(self):
            self.normal_scales = []

        def normal(self, *, loc, scale, size):
            self.normal_scales.append(np.array(scale, dtype=float))
            return np.zeros(size, dtype=float) + loc

        def standard_normal(self, size):
            return np.zeros(size, dtype=float)

    rng = RecordingRng()
    rows = [{"week_index": 0, "model_iv": 0.2, "log_moneyness": 0.0, "maturity_years": 0.25}]
    persistent_factor_noise(rows, rng, config)  # type: ignore[arg-type]

    assert len(rng.normal_scales) == 1
    assert np.allclose(
        rng.normal_scales[0],
        np.sqrt(np.array(config.scenarios["persistent_factor"]["q_diag"])),
    )


def test_persistent_factor_residual_policy_matches_total_marginal_scale():
    config = _noise_settings().scenarios["persistent_factor"]
    lm = np.array([-0.15, -0.075, 0.0, 0.075, 0.15])
    tau = np.array([1.0 / 12.0, 1.0 / 4.0, 1.0 / 2.0])
    grid_lm, grid_tau = np.meshgrid(lm, tau, indexing="ij")
    grid_lm = grid_lm.ravel()
    grid_tau = grid_tau.ravel()

    residual_scale = persistent_factor_residual_scale(grid_lm, grid_tau, config)
    stationary_std = np.array(config["stationary_factor_std"])
    factor_variance = (
        stationary_std[0] * stationary_std[0]
        + grid_lm * grid_lm * stationary_std[1] * stationary_std[1]
        + grid_tau * grid_tau * stationary_std[2] * stationary_std[2]
    )
    total_scale = np.sqrt(factor_variance + residual_scale * residual_scale)

    assert np.allclose(total_scale, marginal_scale(grid_lm, grid_tau, config))


def test_tick_rounding_and_capping_are_recorded():
    config = _noise_settings()
    rows = _clean_rows()
    noisy, _ = generate_noisy_panel_rows(rows, scenario="low_iid", seed=1, config=config)
    assert all(abs(row["price_after_rounding"] / config.tick_size - round(row["price_after_rounding"] / config.tick_size)) < 1e-9 for row in noisy)

    _, _, was_capped, cap_direction = apply_price_mechanics(
        raw_price=10_000.0,
        spot=100.0,
        strike=100.0,
        tau=0.25,
        rate=0.02,
        dividend_yield=0.0,
        option_type="call",
        tick_size=config.tick_size,
        price_epsilon=config.price_epsilon,
    )
    assert was_capped
    assert cap_direction == "upper"
