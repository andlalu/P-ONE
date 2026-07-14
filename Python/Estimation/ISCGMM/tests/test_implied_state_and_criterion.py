import math
from pathlib import Path

import numpy as np
import pytest

from DGPSimulation.heston_simulator import HestonPathSimulator
from DGPSimulation.types import HestonSimConfig
from DGPSimulation.variance_drawers import AndersenQeVarianceDrawer
from Estimation.ISCGMM.cgmm_criterion import CgmmFirstStepCriterion, criterion_diagnostics_to_dict
from Estimation.ISCGMM.config import CcfQuadratureConfig, CgmmConfig, ImpliedStateConfig
from Estimation.ISCGMM.implied_state import imply_heston_variance_path
from Models.Heston.parameters import HestonParameters, HestonPhysicalParameters
from OptionData.io import load_option_panel
from OptionData.panel import OptionPanel
from OptionPricing.clean_panel import generate_clean_option_panel_rows, write_panel
from OptionPricing.cos_basis import FixedCosBasisConfig
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver


def _basis() -> FixedCosBasisConfig:
    return FixedCosBasisConfig(
        maturities=(1.0 / 12.0, 0.25),
        effective_widths=(1.0, 1.5),
        generation_n_cos=64,
        estimation_n_cos=64,
    )


def _true_params() -> tuple[HestonPhysicalParameters, HestonParameters]:
    params_p = HestonPhysicalParameters(
        eta=5.0, kappa=7.0, vbar=0.0225, sigma_v=0.4, rho=-0.5, r=0.02, q=0.0
    )
    theta = HestonParameters(
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
        config=HestonSimConfig(t_week=n_weeks, burnin_days=60, s0=100.0, seed=321),
        variance_drawer=AndersenQeVarianceDrawer(),
    ).simulate()
    basis = _basis()
    rows = generate_clean_option_panel_rows(
        run_id="unit",
        sample_id=0,
        path=path,
        params_p=params_p,
        eta_v=5.0,
        maturities_years=basis.maturities,
        log_moneyness=(-0.10, 0.0, 0.10),
        atm_option_type="call",
        pricing_method="COS",
        iv_method="lets_be_rational",
        pricer=CosOptionPricer(),
        cos_basis=basis,
    )
    return write_panel(
        rows,
        tmp_path / "panel",
        metadata={"sample_id": 0, "scenario": "clean", "cos_basis": basis.generation_metadata()},
        panel_format="csv",
    )


def _state_config(*, v_max: float = 1.0, boundary_tol: float = 1e-5) -> ImpliedStateConfig:
    return ImpliedStateConfig(
        cos_basis=_basis(),
        v_min=1e-8,
        v_max=v_max,
        tol=2e-5,
        max_iter=80,
        boundary_tol=boundary_tol,
    )


def test_implied_state_recovers_clean_variance_and_reports_warm_starts(tmp_path):
    _, theta = _true_params()
    panel = load_option_panel(_write_tiny_panel(tmp_path))
    result = imply_heston_variance_path(theta, panel, _state_config())
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
    panel = load_option_panel(_write_tiny_panel(tmp_path), max_dates=4)
    result = imply_heston_variance_path(theta, panel, _state_config(v_max=0.01, boundary_tol=1e-4))
    assert np.any(result.boundary_hit)


def test_implied_state_rejects_panel_configuration_cos_basis_mismatch(tmp_path):
    _, theta = _true_params()
    panel = load_option_panel(_write_tiny_panel(tmp_path), max_dates=3)
    mismatched = FixedCosBasisConfig(_basis().maturities, (1.1, 1.5), 64, 64)
    panel = OptionPanel(
        panel.dates,
        metadata={**panel.metadata, "cos_basis": mismatched.generation_metadata()},
    )
    with pytest.raises(ValueError, match="width mismatch"):
        imply_heston_variance_path(theta, panel, _state_config())


def test_fixed_width_inversion_reuses_heston_coefficients_across_dates(tmp_path, monkeypatch):
    _, theta = _true_params()
    panel = load_option_panel(_write_tiny_panel(tmp_path), max_dates=4)
    calls = []
    original = HestonAnalyticCcfSolver.solve_coefficients

    def counting_solve(self, *args, **kwargs):
        calls.append(1)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(HestonAnalyticCcfSolver, "solve_coefficients", counting_solve)
    result = imply_heston_variance_path(theta, panel, _state_config())
    n_maturities = len(np.unique(np.concatenate([date.maturities for date in panel.dates])))
    assert result.success_rate == pytest.approx(1.0)
    assert len(calls) == n_maturities
    assert result.coefficient_solve_count == n_maturities


def test_first_step_criterion_is_deterministic_and_preserves_ordering(tmp_path):
    _, theta = _true_params()
    panel = load_option_panel(_write_tiny_panel(tmp_path), max_dates=6)
    config = CgmmConfig(
        implied_state=_state_config(),
        quadrature=CcfQuadratureConfig(order=2, scale=1.0),
        instrument_precision=(10.0, 50.0),
        transition_rk_steps=16,
    )
    criterion = CgmmFirstStepCriterion(panel, config)
    diagnostics = criterion.evaluate(theta, return_diagnostics=True)
    q_true_1 = diagnostics.criterion_value
    q_true_2 = criterion.evaluate(theta)
    assert q_true_1 == pytest.approx(q_true_2, rel=0.0, abs=1e-14)
    assert math.isfinite(q_true_1)
    state_info = criterion_diagnostics_to_dict(theta, diagnostics)["state_inversion"]
    assert state_info["coefficient_solve_count"] > 0
    assert state_info["coefficient_cache_hits"] > 0
    assert state_info["fixed_coefficient_count"] == state_info["coefficient_solve_count"]

    bad_kappa = theta.kappa * 0.55
    bad_kappa_q = theta.kappa_q * 1.60
    perturbed = HestonParameters(
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
    criterion.evaluate(theta)
    np.testing.assert_allclose(panel.log_spot, before)
