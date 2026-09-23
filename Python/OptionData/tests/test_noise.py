import math
from dataclasses import replace

import numpy as np
import pytest

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_price, black76_vega
from Models.Heston.parameters import HestonPhysicalParameters
from OptionData.add_noise import generate_noisy_panel_rows, validate_noise_settings
from OptionData.noise_common import (
    NoiseSettings,
    apply_price_mechanics,
    apply_price_projection,
    marginal_scale,
    price_bounds,
    scenario_seed,
)
from OptionData.noise_persistent_factor import (
    compute_q_diag_from_stationary_std,
    persistent_factor_noise,
    persistent_factor_residual_scale,
    validate_persistent_factor_settings,
)
from OptionData.noise_variance_linked import (
    variance_linked_gamma,
    variance_linked_residual_scale,
)


def _noise_settings() -> NoiseSettings:
    return NoiseSettings(
        base_seed=9900000,
        sigma_min=0.0001,
        price_epsilon=1e-10,
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
            "variance_linked_factor": {
                "alpha_0": 0.0015,
                "alpha_m": 2.0,
                "alpha_tau": 0.05,
                "tau_min": 0.01984126984126984,
                "a_diag": [0.85, 0.65, 0.65],
                "stationary_factor_std": [0.0007, 0.0025, 0.0008],
                "q_diag": [1.35975e-7, 3.609375e-6, 3.696e-7],
                "residual_policy": "match_total_marginal_scale",
                "factor_initialization": "zero",
                "latent_variance_share": 0.5,
            },
        },
    )


def _physical_parameters() -> HestonPhysicalParameters:
    return HestonPhysicalParameters(
        eta=1.5,
        kappa=3.0,
        vbar=0.04,
        sigma_v=0.4,
        rho=-0.7,
        r=0.02,
        q=0.0,
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
                        "V": 0.035 + 0.01 * week_index,
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


def test_run_001_tick_rounding_is_optional_and_precedes_projection():
    clean = _clean_rows()
    unrounded, _ = generate_noisy_panel_rows(
        clean, scenario="low_iid", seed=123, config=_noise_settings()
    )
    tick_config = replace(_noise_settings(), tick_size=0.01)
    rounded, _ = generate_noisy_panel_rows(
        clean, scenario="low_iid", seed=123, config=tick_config
    )
    assert [row["raw_noisy_iv"] for row in rounded] == [
        row["raw_noisy_iv"] for row in unrounded
    ]
    assert any(
        row["observed_price"] != plain["observed_price"]
        for row, plain in zip(rounded, unrounded)
    )
    for row in rounded:
        expected = apply_price_mechanics(
            raw_price=row["raw_price_before_rounding"],
            spot=row["S"],
            strike=row["strike"],
            tau=row["maturity_years"],
            rate=row["r"],
            dividend_yield=row["q"],
            option_type=row["option_type"],
            tick_size=0.01,
            price_epsilon=tick_config.price_epsilon,
        )
        assert row["raw_noisy_price"] == row["raw_price_before_rounding"]
        assert row["price_after_rounding"] == expected[0]
        assert row["observed_price"] == expected[1]
        assert row["was_price_capped"] == expected[2]
        assert row["cap_direction"] == expected[3]


def test_invalid_tick_size_is_rejected():
    with pytest.raises(ValueError, match="tick_size"):
        validate_noise_settings(replace(_noise_settings(), tick_size=0.0))


def test_spatial_corr_noise_has_contract_level_draws():
    config = _noise_settings()
    rows, _ = generate_noisy_panel_rows(_clean_rows(), scenario="spatial_corr", seed=456, config=config)
    draws = np.array([row["noise_draw"] for row in rows])
    assert np.std(draws) > 0.0
    assert all(row["noise_scenario"] == "spatial_corr" for row in rows)


@pytest.mark.parametrize(
    ("scenario", "seed"),
    [("low_iid", 101), ("spatial_corr", 202), ("persistent_factor", 303)],
)
def test_existing_noise_streams_do_not_depend_on_design_d_configuration(
    scenario,
    seed,
):
    config = _noise_settings()
    without_design_d = NoiseSettings(
        base_seed=config.base_seed,
        sigma_min=config.sigma_min,
        price_epsilon=config.price_epsilon,
        scenarios={
            name: settings
            for name, settings in config.scenarios.items()
            if name != "variance_linked_factor"
        },
    )

    with_d, _ = generate_noisy_panel_rows(
        _clean_rows(),
        scenario=scenario,
        seed=seed,
        config=config,
    )
    without_d, _ = generate_noisy_panel_rows(
        _clean_rows(),
        scenario=scenario,
        seed=seed,
        config=without_design_d,
    )

    np.testing.assert_array_equal(
        [row["noise_draw"] for row in with_d],
        [row["noise_draw"] for row in without_d],
    )
    np.testing.assert_array_equal(
        [row["raw_noisy_iv"] for row in with_d],
        [row["raw_noisy_iv"] for row in without_d],
    )


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


def test_price_projection_preserves_valid_price_and_caps_bounds():
    arguments = {
        "spot": 100.0,
        "strike": 100.0,
        "tau": 0.25,
        "rate": 0.02,
        "dividend_yield": 0.0,
        "option_type": "call",
        "price_epsilon": 1e-10,
    }
    lower, upper = price_bounds(
        arguments["spot"],
        arguments["strike"],
        arguments["tau"],
        arguments["rate"],
        arguments["dividend_yield"],
        arguments["option_type"],
    )

    observed, was_capped, direction = apply_price_projection(
        raw_price=5.123456789,
        **arguments,
    )
    assert observed == 5.123456789
    assert not was_capped
    assert direction == "none"

    observed, was_capped, direction = apply_price_projection(
        raw_price=10_000.0,
        **arguments,
    )
    assert observed == pytest.approx(upper - arguments["price_epsilon"])
    assert was_capped
    assert direction == "upper"

    observed, was_capped, direction = apply_price_projection(
        raw_price=-1.0,
        **arguments,
    )
    assert observed == pytest.approx(lower + arguments["price_epsilon"])
    assert was_capped
    assert direction == "lower"


def test_variance_linked_factor_is_deterministic_and_matches_true_variance():
    config = _noise_settings()
    rows = _clean_rows()
    params = _physical_parameters()
    first, first_factors = generate_noisy_panel_rows(
        rows,
        scenario="variance_linked_factor",
        seed=789,
        config=config,
        params_p=params,
    )
    second, second_factors = generate_noisy_panel_rows(
        rows,
        scenario="variance_linked_factor",
        seed=789,
        config=config,
        params_p=params,
    )
    assert [row["noise_draw"] for row in first] == [
        row["noise_draw"] for row in second
    ]
    assert first_factors == second_factors

    factor_config = config.scenarios["variance_linked_factor"]
    gamma_v = variance_linked_gamma(
        np.array([row["log_moneyness"] for row in rows]),
        np.array([row["maturity_years"] for row in rows]),
        factor_config,
    )
    variance_state_std = math.sqrt(
        params.vbar * params.sigma_v**2 / (2.0 * params.kappa)
    )
    for factor in first_factors:
        week = int(factor["week_index"])
        variance = next(row["V"] for row in rows if row["week_index"] == week)
        assert factor["variance_linked_factor"] == pytest.approx(
            gamma_v * (variance - params.vbar) / variance_state_std,
            rel=1e-14,
            abs=1e-16,
        )


@pytest.mark.parametrize("share", [float("nan"), float("inf"), -0.1, 0.0, 1.0, 1.1])
def test_variance_linked_share_must_be_finite_and_strictly_between_zero_and_one(
    share,
):
    config = _noise_settings()
    config.scenarios["variance_linked_factor"]["latent_variance_share"] = share
    with pytest.raises(ValueError, match="latent_variance_share"):
        validate_noise_settings(config)


def test_variance_linked_calibration_and_residual_scale_match_production_target():
    config = _noise_settings().scenarios["variance_linked_factor"]
    log_moneyness, maturities = np.meshgrid(
        np.array([-0.15, -0.075, 0.0, 0.075, 0.15]),
        np.array([1.0 / 12.0, 0.25, 0.5]),
        indexing="ij",
    )
    log_moneyness = log_moneyness.ravel()
    maturities = maturities.ravel()
    gamma_v = variance_linked_gamma(log_moneyness, maturities, config)
    total_scale = marginal_scale(log_moneyness, maturities, config)
    shares = gamma_v * gamma_v / (total_scale * total_scale)
    assert gamma_v == pytest.approx(0.00135463, abs=5e-9)
    assert float(np.mean(shares)) == pytest.approx(0.5, abs=1e-15)

    residual_scale = variance_linked_residual_scale(
        log_moneyness,
        maturities,
        config,
        gamma_v,
    )
    stationary_std = np.array(config["stationary_factor_std"])
    factor_variance = (
        stationary_std[0] ** 2
        + log_moneyness**2 * stationary_std[1] ** 2
        + maturities**2 * stationary_std[2] ** 2
    )
    assert np.all(residual_scale >= 0.0)
    np.testing.assert_allclose(
        factor_variance + gamma_v**2 + residual_scale**2,
        total_scale**2,
        rtol=1e-14,
        atol=1e-18,
    )


def test_variance_linked_residual_scale_rejects_infeasible_calibration():
    config = _noise_settings().scenarios["variance_linked_factor"]
    with pytest.raises(ValueError, match="negative residual variance"):
        variance_linked_residual_scale(
            np.array([0.0]),
            np.array([0.25]),
            config,
            gamma_v=0.01,
        )


def test_variance_linked_factor_requires_physical_parameters():
    with pytest.raises(ValueError, match="requires physical Heston parameters"):
        generate_noisy_panel_rows(
            _clean_rows(),
            scenario="variance_linked_factor",
            seed=789,
            config=_noise_settings(),
        )


def test_scenario_seed_offsets_preserve_existing_streams_and_add_design_d():
    assert scenario_seed(9_900_000, 2, "low_iid") == 9_902_101
    assert scenario_seed(9_900_000, 2, "spatial_corr") == 9_902_202
    assert scenario_seed(9_900_000, 2, "persistent_factor") == 9_902_303
    assert scenario_seed(9_900_000, 2, "variance_linked_factor") == 9_902_404
