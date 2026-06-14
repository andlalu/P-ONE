import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from DGPSimulation.heston_simulator import HestonPathSimulator
from DGPSimulation.types import HestonParamsP, HestonSimConfig
from DGPSimulation.variance_drawers import AndersenQeVarianceDrawer
from Estimation.ISCGMM.cgmm_criterion import CgmmFirstStepCriterion
from Estimation.ISCGMM.implied_state import imply_heston_variance_path
from Estimation.ISCGMM.panel import load_option_panel_data
from Estimation.ISCGMM.types import CgmmConfig, HestonEstimationParams, ImpliedStateConfig, QuadratureConfig
from OptionPricing.clean_panel import generate_clean_option_panel_rows
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.types import CosPricingConfig


def _true_params() -> tuple[HestonParamsP, HestonEstimationParams]:
    params_p = HestonParamsP(
        eta=5.0,
        kappa=7.0,
        vbar=0.0225,
        sigma_v=0.4,
        rho=-0.5,
        r=0.02,
        q=0.0,
    )
    theta = HestonEstimationParams(
        eta=params_p.eta,
        kappa=params_p.kappa,
        vbar=params_p.vbar,
        sigma_v=params_p.sigma_v,
        rho=params_p.rho,
        eta_v=5.0,
        r=params_p.r,
        q=params_p.q,
    )
    return params_p, theta


def _write_tiny_panel(tmp_path: Path, *, n_weeks: int = 7) -> Path:
    params_p, _ = _true_params()
    path = HestonPathSimulator(
        params=params_p,
        config=HestonSimConfig(
            delta=1.0 / 252.0,
            m_week=5,
            t_week=n_weeks,
            burnin_days=60,
            s0=100.0,
            seed=321,
        ),
        variance_drawer=AndersenQeVarianceDrawer(),
    ).simulate()
    rows = generate_clean_option_panel_rows(
        run_id="unit",
        sample_id=0,
        path=path,
        params_p=params_p,
        eta_v=5.0,
        maturities_years=(1.0 / 12.0, 0.25),
        log_moneyness=(-0.10, 0.0, 0.10),
        atm_option_type="call",
        pricing_method="COS",
        iv_method="lets_be_rational",
        pricer=CosOptionPricer(),
        pricing_config=CosPricingConfig(n_cos=64, truncation_width=10.0),
    )
    target = tmp_path / "panel.parquet"
    pd.DataFrame(rows).to_parquet(target, index=False)
    return target


def _state_config(effective_truncation_width: float | None = None) -> ImpliedStateConfig:
    return ImpliedStateConfig(
        tol=2e-5,
        max_iter=80,
        effective_truncation_width=effective_truncation_width,
        pricing_config=CosPricingConfig(n_cos=64, truncation_width=10.0),
    )


def _generation_effective_width(panel) -> float:
    true_v = panel.true_variance
    assert true_v is not None
    maturities = np.unique(np.concatenate([date.maturities for date in panel.dates]))
    return CosOptionPricer.effective_truncation_width(
        true_v,
        maturities,
        CosPricingConfig(n_cos=64, truncation_width=10.0),
    )


def test_implied_state_recovers_clean_variance_and_reports_warm_starts(tmp_path):
    _, theta = _true_params()
    panel = load_option_panel_data(_write_tiny_panel(tmp_path))
    result = imply_heston_variance_path(theta, panel, _state_config(_generation_effective_width(panel)))
    true_v = panel.true_variance
    assert true_v is not None
    rmse = math.sqrt(float(np.mean((result.variance - true_v) ** 2)))
    assert rmse < 5e-5
    assert np.max(np.abs(result.variance - true_v)) < 1e-4
    assert result.success_rate == pytest.approx(1.0)
    np.testing.assert_allclose(result.start_values[0], theta.vbar)
    np.testing.assert_allclose(result.start_values[1:], result.variance[:-1])


def test_implied_state_reports_boundary_hits(tmp_path):
    _, theta = _true_params()
    panel = load_option_panel_data(_write_tiny_panel(tmp_path), max_dates=4)
    cfg = ImpliedStateConfig(
        v_min=1e-8,
        v_max=0.01,
        tol=2e-4,
        boundary_tol=1e-4,
        max_iter=40,
        pricing_config=CosPricingConfig(n_cos=48, truncation_width=10.0),
    )
    result = imply_heston_variance_path(theta, panel, cfg)
    assert np.any(result.boundary_hit)


def test_first_step_criterion_is_deterministic_finite_and_penalizes_strong_perturbation(tmp_path):
    _, theta = _true_params()
    panel = load_option_panel_data(_write_tiny_panel(tmp_path), max_dates=6)
    config = CgmmConfig(
        implied_state=_state_config(_generation_effective_width(panel)),
        quadrature=QuadratureConfig(order=2, scale=1.0),
        instrument_precision=(10.0, 50.0),
        transition_rk_steps=16,
    )
    criterion = CgmmFirstStepCriterion(panel, config)
    q_true_1 = criterion.evaluate(theta)
    q_true_2 = criterion.evaluate(theta)
    assert q_true_1 == pytest.approx(q_true_2, rel=0.0, abs=1e-14)
    assert math.isfinite(q_true_1)

    bad_kappa = theta.kappa * 0.55
    bad_kappa_q = theta.kappa_q * 1.60
    perturbed = HestonEstimationParams(
        eta=theta.eta * 1.4,
        kappa=bad_kappa,
        vbar=theta.vbar * 1.6,
        sigma_v=theta.sigma_v * 0.55,
        rho=-0.9,
        eta_v=bad_kappa - bad_kappa_q,
        r=theta.r,
        q=theta.q,
    )
    q_bad = criterion.evaluate(perturbed)
    assert math.isfinite(q_bad)
    assert q_true_1 < q_bad

    before = panel.log_spot.copy()
    _ = criterion.evaluate(theta)
    np.testing.assert_allclose(panel.log_spot, before)
