import math

import numpy as np
from scipy.integrate import quad

from Scripts.make_is_cgmm_report_assets import (
    HORIZON_YEARS,
    _build_noise_mechanics_frame,
    _validate_estimate,
    average_variance,
    paired_path_metrics,
    risk_premium_paths,
)


def test_noise_mechanics_keep_raw_error_and_lower_adjustments_separate():
    cell = (1.0 / 12.0, 0.15)
    raw = {
        ("low_iid", *cell): [np.array([0.001, -0.002])],
        ("spatial_corr", *cell): [np.array([0.003, 0.0])],
        ("persistent_factor", *cell): [np.array([-0.004, 0.001])],
    }
    lower = {
        ("low_iid", *cell): [np.array([True, False])],
        ("spatial_corr", *cell): [np.array([True, True])],
        ("persistent_factor", *cell): [np.array([False, False])],
    }

    diagnostics = _build_noise_mechanics_frame(raw, lower)
    pooled = diagnostics.loc[diagnostics["scenario"] == "pooled_noisy"].iloc[0]
    raw_bp = np.array([10.0, -20.0, 30.0, 0.0, -40.0, 10.0])
    assert pooled["q95_absolute_raw_error_bp"] == np.quantile(np.abs(raw_bp), 0.95)
    assert pooled["lower_bound_adjustment_percent"] == 50.0
    assert pooled["observations"] == 6


def test_average_variance_limit_and_direct_numerical_integration():
    variance = np.array([0.01, 0.0225, 0.08])
    kappa = 7.0
    vbar = 0.0225

    np.testing.assert_allclose(
        average_variance(variance, kappa=kappa, vbar=vbar, h=0.0),
        variance,
        rtol=0.0,
        atol=0.0,
    )
    near_zero = average_variance(variance, kappa=kappa, vbar=vbar, h=1.0e-10)
    np.testing.assert_allclose(near_zero, variance, rtol=1.0e-9, atol=1.0e-12)

    analytic = average_variance(variance, kappa=kappa, vbar=vbar, h=HORIZON_YEARS)
    numerical = np.array(
        [
            quad(lambda s: vbar + (value - vbar) * math.exp(-kappa * s), 0.0, HORIZON_YEARS)[0]
            / HORIZON_YEARS
            for value in variance
        ]
    )
    np.testing.assert_allclose(analytic, numerical, rtol=2.0e-13, atol=2.0e-15)


def test_risk_premium_paths_use_physical_minus_risk_neutral_variance():
    variance = np.array([0.01, 0.0225, 0.05])
    erp, vrp = risk_premium_paths(
        variance,
        eta=5.0,
        kappa=7.0,
        vbar=0.0225,
        kappa_q=2.0,
    )
    a_p = average_variance(variance, kappa=7.0, vbar=0.0225, h=HORIZON_YEARS)
    a_q = average_variance(
        variance,
        kappa=2.0,
        vbar=7.0 * 0.0225 / 2.0,
        h=HORIZON_YEARS,
    )
    np.testing.assert_allclose(erp, 5.0 * a_p)
    np.testing.assert_allclose(vrp, a_p - a_q)


def test_paired_path_metrics_are_estimated_minus_true():
    truth = np.array([1.0, 2.0, 4.0])
    estimate = np.array([1.5, 1.0, 5.5])
    metrics = paired_path_metrics(estimate, truth)
    error = estimate - truth
    assert metrics["mean_error"] == np.mean(error)
    assert metrics["mae"] == np.mean(np.abs(error))
    assert metrics["rmse"] == np.sqrt(np.mean(error * error))
    assert metrics["path_average_truth"] == np.mean(truth)
    assert metrics["path_average_estimate"] == np.mean(estimate)


def test_complete_finite_estimate_is_included_without_using_internal_search_flag():
    implied = {
        "variance": np.full(526, 0.02).tolist(),
        "objective": np.zeros(526).tolist(),
    }
    result = {
        "success": False,
        "final_criterion": 1.0e-6,
        "estimated_parameters": {
            "eta": 5.0,
            "kappa": 7.0,
            "vbar": 0.0225,
            "sigma_v": 0.4,
            "rho": -0.5,
            "eta_v": 5.0,
            "r": 0.02,
            "q": 0.0,
        },
        "final_diagnostics": {
            "n_dates": 526,
            "n_transitions": 524,
            "n_contracts": 7890,
            "implied_state": implied,
        },
    }
    _validate_estimate(result, sample_id=0, scenario="clean")
