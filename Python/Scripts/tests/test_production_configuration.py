import json
from pathlib import Path

from OptionPricing.cos_basis import FixedCosBasisConfig
from Scripts.generation import load_config


CONFIG_DIRECTORY = Path(__file__).resolve().parents[1] / "configs"
CALIBRATION_DIRECTORY = Path(__file__).resolve().parents[3] / "outputs" / "calibration" / "fixed_cos"


def test_frozen_production_solver_matches_deterministic_benchmark():
    profile = json.loads((CONFIG_DIRECTORY / "is_cgmm_production_profile.json").read_text())
    benchmark = json.loads(
        (CALIBRATION_DIRECTORY / "state_solver_benchmark_summary.json").read_text()
    )
    assert benchmark["deterministic_seeds"] == {
        "low_iid": 301,
        "spatial_corr": 302,
        "persistent_factor": 303,
    }
    assert profile["implied_state"]["state_solver"] == benchmark["selected_production_solver"]
    assert profile["implied_state"]["fallback_solver"] == benchmark["fallback_solver"]
    assert benchmark["aggregate"]["bounded_brent"]["failure_count"] == 0


def test_generation_and_estimation_load_one_authoritative_basis():
    generation = load_config(CONFIG_DIRECTORY / "generation_run_001.json")
    profile = json.loads((CONFIG_DIRECTORY / "is_cgmm_production_profile.json").read_text())
    estimation_basis = FixedCosBasisConfig.load(CONFIG_DIRECTORY / profile["cos_basis_path"])
    assert generation.cos_basis_path == str(
        (CONFIG_DIRECTORY / "heston_cos_basis_production.json").resolve()
    )
    assert generation.cos_basis == estimation_basis
    assert generation.cos_basis.stable_hash() == estimation_basis.stable_hash()
