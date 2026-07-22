import json
from dataclasses import asdict
from pathlib import Path

import pytest

from Scripts.experiment_config import load_experiment_config

CONFIG_DIRECTORY = Path(__file__).resolve().parents[1] / "configs"
PRODUCTION_CONFIG = CONFIG_DIRECTORY / "heston_experiment_run_001.json"


def test_frozen_production_experiment_values():
    experiment = load_experiment_config(PRODUCTION_CONFIG)

    assert experiment.run_id == "run_001"
    assert experiment.n_samples == 100
    assert experiment.base_seed == 1234500
    assert experiment.panel_format == "parquet"
    assert experiment.workers is None
    assert asdict(experiment.dgp) == {
        "eta": 5.0,
        "kappa": 7.0,
        "vbar": 0.0225,
        "sigma_v": 0.4,
        "rho": -0.5,
        "r": 0.02,
        "q": 0.0,
    }
    assert experiment.simulation.delta == 1.0 / 252.0
    assert experiment.simulation.m_week == 5
    assert experiment.simulation.t_week == 525
    assert experiment.simulation.burnin_days == 756
    assert experiment.simulation.s0 == 100.0
    assert experiment.simulation.return_daily is True
    assert experiment.cos_basis.maturities == (1.0 / 12.0, 0.25, 0.5)
    assert experiment.cos_basis.effective_widths == (0.75, 1.25, 2.0)
    assert experiment.cos_basis.generation_n_cos == 576
    assert experiment.cos_basis.estimation_n_cos == 576
    assert experiment.noise is not None
    assert experiment.noise.scenario_names() == ("low_iid", "spatial_corr", "persistent_factor")

    state = experiment.criterion_config.implied_state
    assert state.state_solver == "bounded_brent"
    assert state.fallback_solver == "golden_section"
    assert (state.v_min, state.v_max, state.tol, state.max_iter) == (1e-6, 0.15, 1e-7, 100)

    optimizer = experiment.optimizer_config
    assert asdict(optimizer.base_start) == {
        "eta": 4.8,
        "kappa": 6.8,
        "vbar": 0.024,
        "sigma_v": 0.42,
        "rho": -0.52,
        "eta_v": 4.8,
        "r": 0.02,
        "q": 0.0,
    }
    assert optimizer.natural_bounds == (
        (2.0, 8.0),
        (4.0, 9.0),
        (0.012, 0.04),
        (0.2, 0.65),
        (-0.8, -0.2),
        (1.5, 4.0),
    )
    assert optimizer.stage1.max_evaluations == 120
    assert (optimizer.stage1.xtol, optimizer.stage1.ftol) == (0.02, 0.002)
    assert optimizer.stage2.max_evaluations == 300
    assert (optimizer.stage2.xtol, optimizer.stage2.ftol) == (0.0002, 0.00002)
    assert optimizer.progress_every == 10


def test_mismatched_cos_maturities_and_widths_are_rejected(tmp_path):
    payload = json.loads(PRODUCTION_CONFIG.read_text())
    payload["cos"]["effective_widths"] = [0.75, 1.25]
    config_path = tmp_path / "invalid_experiment.json"
    config_path.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="one-to-one"):
        load_experiment_config(config_path)
