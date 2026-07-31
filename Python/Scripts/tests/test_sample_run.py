import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from DGPSimulation.heston_simulator import HestonPathSimulator
from DGPSimulation.io import load_heston_path_npz
from DGPSimulation.variance_drawers import AndersenQeVarianceDrawer
from Estimation.ISCGMM.estimate import estimate_first_step
from OptionData.add_noise import generate_noisy_panel_rows
from OptionData.clean_panel import generate_clean_option_panel_rows
from OptionData.io import load_option_panel, write_records
from OptionData.noise_common import NOISE_SCENARIOS, scenario_seed
from OptionPricing.cos_basis import cos_specification_metadata
from OptionPricing.cos_pricer import CosOptionPricer
from Scripts.experiment_config import load_experiment_config
from Scripts.sample_run import (
    COMBINED_PANEL_COLUMNS,
    SCENARIO_ORDER,
    _estimate_is_complete,
    atomic_json,
    initialise_or_verify_run,
    run_sample,
    sha256_file,
)

PRODUCTION_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "heston_experiment_run_001.json"


def _mini_config(tmp_path: Path, *, n_samples: int = 2):
    payload = json.loads(PRODUCTION_CONFIG.read_text())
    payload["run"].update(
        run_id="test_run",
        n_samples=n_samples,
        base_seed=9000,
        output_root="configured_output",
        panel_format="parquet",
        workers=1,
    )
    payload["dgp"] = {
        "eta": 1.5,
        "kappa": 3.0,
        "vbar": 0.04,
        "sigma_v": 0.4,
        "rho": -0.7,
        "r": 0.02,
        "q": 0.0,
    }
    payload["simulation"].update(t_week=4, burnin_days=2, return_daily=True)
    payload["q_measure"]["eta_v"] = 0.0
    payload["panel"]["log_moneyness"] = [-0.05, 0.0, 0.05]
    payload["cos"].update(
        maturities_years=[0.25],
        effective_widths=[1.5],
        generation_n_cos=64,
        estimation_n_cos=32,
    )
    payload["estimation"]["optimizer"]["coarse_pass"]["max_evaluations"] = 12
    payload["estimation"]["optimizer"]["refinement_pass"]["max_evaluations"] = 16
    payload["estimation"]["optimizer"]["progress_every"] = 1000
    config_path = tmp_path / "experiment.json"
    config_path.write_text(json.dumps(payload))
    return load_experiment_config(config_path)


def _fake_estimate(_panel_file, scenario, _config, *, max_dates):
    return {
        "scenario": scenario,
        "max_dates": max_dates,
        "success": True,
        "final_criterion": 0.25,
        "estimated_parameters": {
            "eta": 4.8,
            "kappa": 6.8,
            "vbar": 0.024,
            "sigma_v": 0.42,
            "rho": -0.52,
            "eta_v": 4.8,
            "r": 0.02,
            "q": 0.0,
        },
        "free_parameters": [4.8, 1.9, -3.7, -0.87, -0.57, 0.69],
        "function_evaluations": 12,
        "penalty_evaluations": 0,
        "powell_passes": [
            {"pass_name": "coarse"},
            {"pass_name": "refinement"},
        ],
        "final_diagnostics": {"criterion_value": 0.25},
        "total_runtime_seconds": 1.0,
    }


def test_atomic_json_and_checksum(tmp_path):
    target = tmp_path / "record.json"
    atomic_json(target, {"status": "initialising", "value": 3})
    assert json.loads(target.read_text()) == {"status": "initialising", "value": 3}
    assert sha256_file(target) == sha256_file(target)
    assert not list(tmp_path.glob("*.tmp.json"))


def test_run_record_creation_and_verification(tmp_path):
    config = _mini_config(tmp_path)
    root = tmp_path / "run"
    first = initialise_or_verify_run(root, config)
    second = initialise_or_verify_run(root, config)
    assert first == second
    assert first["run_id"] == "test_run"
    assert first["configuration_hash"] == config.experiment_config_hash
    assert first["configuration"] == config.raw_config
    assert set(first["environment"]) == {"python", "numpy", "scipy", "pandas", "pyarrow"}
    assert set(path.name for path in root.iterdir()) == {"run.json"}

    changed = json.loads(Path(config.source_path).read_text())
    changed["run"]["base_seed"] += 1
    changed_path = tmp_path / "changed.json"
    changed_path.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="configuration hash"):
        initialise_or_verify_run(root, load_experiment_config(changed_path))


def test_generation_combines_scenarios_validates_and_resumes(tmp_path):
    config = _mini_config(tmp_path)
    root = tmp_path / "run"
    record = run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        overwrite=True,
        generation_only=True,
        log_level="ERROR",
    )
    sample = root / "sample_000"
    assert {path.name for path in root.iterdir()} == {"run.json", "sample_000"}
    assert {path.name for path in sample.iterdir()} == {
        "path.npz",
        "panels.parquet",
        "record.json",
        "sample.log",
    }
    assert record["status"] == "generated"
    assert record["current_stage"] == "validated"
    assert record["validation"]["passed"]
    assert set(record["validation"]["panels"]["scenarios"]) == set(SCENARIO_ORDER)
    assert record["artifacts"]["panels"]["rows"] == 60
    assert record["artifacts"]["panels"]["rows_by_scenario"] == {
        name: 15 for name in SCENARIO_ORDER
    }
    assert not list(root.rglob("*.tmp*"))

    path_hash = record["artifacts"]["path"]["sha256"]
    panel_hash = record["artifacts"]["panels"]["sha256"]
    resumed = run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        resume=True,
        generation_only=True,
        log_level="ERROR",
    )
    assert resumed["artifacts"]["path"]["sha256"] == path_hash
    assert resumed["artifacts"]["panels"]["sha256"] == panel_hash
    assert resumed["status"] == "generated"


def test_estimate_completion_requires_a_complete_result_not_optimizer_success():
    complete = _fake_estimate(None, "clean", None, max_dates=None)
    complete["success"] = False
    assert _estimate_is_complete(complete)

    for missing in (
        "final_criterion",
        "estimated_parameters",
        "free_parameters",
        "function_evaluations",
        "penalty_evaluations",
        "final_diagnostics",
        "total_runtime_seconds",
    ):
        incomplete = dict(complete)
        incomplete.pop(missing)
        assert not _estimate_is_complete(incomplete)

    incomplete = dict(complete)
    incomplete["powell_passes"] = [{"pass_name": "coarse"}]
    assert not _estimate_is_complete(incomplete)


def test_combined_panel_schema_values_factors_and_scenario_loading(tmp_path):
    config = _mini_config(tmp_path)
    root = tmp_path / "run"
    run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        overwrite=True,
        generation_only=True,
        log_level="ERROR",
    )
    target = root / "sample_000" / "panels.parquet"
    frame = pd.read_parquet(target)
    assert tuple(frame.columns) == COMBINED_PANEL_COLUMNS
    assert tuple(frame["scenario"].drop_duplicates()) == SCENARIO_ORDER
    clean = frame[frame["scenario"] == "clean"]
    assert np.array_equal(clean["estimation_price"], clean["model_price"])
    assert np.array_equal(clean["estimation_iv"], clean["model_iv"])
    assert clean["noise_draw"].isna().all()
    assert clean["noise_seed"].isna().all()
    assert not clean["was_price_capped"].any()
    assert (clean["cap_direction"] == "none").all()

    persistent = frame[frame["scenario"] == "persistent_factor"]
    factor_columns = [
        "persistent_factor_level",
        "persistent_factor_moneyness",
        "persistent_factor_maturity",
    ]
    assert persistent[factor_columns].notna().all().all()
    assert (
        persistent.groupby("week_index")[factor_columns].nunique().to_numpy() == 1
    ).all()
    for scenario in SCENARIO_ORDER:
        panel = load_option_panel(target, scenario=scenario)
        assert panel.n_dates == 5
        assert panel.n_contracts == 15
        assert panel.metadata["scenario"] == scenario
        assert panel.metadata["sample_id"] == 0
        selected = frame[frame["scenario"] == scenario]
        assert np.array_equal(
            np.concatenate([date.observed_iv for date in panel.dates]),
            selected["estimation_iv"].to_numpy(),
        )
    with pytest.raises(ValueError, match="unknown panel scenario"):
        load_option_panel(target, scenario="missing")
    with pytest.raises(ValueError, match="scenario must be selected"):
        load_option_panel(target)


def test_failure_recording_and_overwrite_isolation(tmp_path, monkeypatch):
    config = _mini_config(tmp_path)
    root = tmp_path / "run"
    other = root / "sample_001"
    other.mkdir(parents=True)
    sentinel = other / "user-owned.txt"
    sentinel.write_text("preserve")

    def fail_build(*_args, **_kwargs):
        raise RuntimeError("deliberate panel failure")

    monkeypatch.setattr("Scripts.sample_run.build_combined_panel", fail_build)
    with pytest.raises(RuntimeError, match="deliberate panel failure"):
        run_sample(
            config=config,
            sample_id=0,
            output_root=root,
            overwrite=True,
            generation_only=True,
            log_level="ERROR",
        )
    record = json.loads((root / "sample_000" / "record.json").read_text())
    assert record["status"] == "failed"
    assert record["current_stage"] == "generating"
    assert record["errors"][-1]["exception_type"] == "RuntimeError"
    assert "deliberate panel failure" in record["errors"][-1]["traceback"]
    assert sentinel.read_text() == "preserve"


def test_complete_sample_skip_and_partial_estimation_resume(tmp_path, monkeypatch):
    config = _mini_config(tmp_path)
    root = tmp_path / "run"
    calls: list[str] = []

    def recording_estimate(panel_file, scenario, config, *, max_dates):
        calls.append(scenario)
        result = _fake_estimate(panel_file, scenario, config, max_dates=max_dates)
        result["success"] = scenario != "clean"
        result["status"] = 1 if scenario == "clean" else 0
        result["message"] = (
            "Maximum number of function evaluations has been exceeded."
            if scenario == "clean"
            else "Optimization terminated successfully."
        )
        return result

    monkeypatch.setattr("Scripts.sample_run.estimate_scenario", recording_estimate)
    partial = run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        overwrite=True,
        scenario="clean",
        log_level="ERROR",
    )
    assert partial["status"] == "estimating"
    assert calls == ["clean"]

    second_partial = run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        resume=True,
        scenario="low_iid",
        log_level="ERROR",
    )
    assert second_partial["status"] == "estimating"
    assert calls == ["clean", "low_iid"]

    complete = run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        resume=True,
        log_level="ERROR",
    )
    assert complete["status"] == "complete"
    assert calls == ["clean", "low_iid", "spatial_corr", "persistent_factor"]
    assert tuple(complete["estimation"]) == SCENARIO_ORDER
    assert complete["estimation"]["clean"]["success"] is False
    assert complete["errors"] == []
    assert all(len(item["powell_passes"]) == 2 for item in complete["estimation"].values())

    monkeypatch.setattr(
        "Scripts.sample_run.estimate_scenario",
        lambda *_args, **_kwargs: pytest.fail("complete sample should not be estimated again"),
    )
    skipped = run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        resume=True,
        log_level="ERROR",
    )
    assert skipped["status"] == "complete"
    assert {path.name for path in (root / "sample_000").iterdir()} == {
        "path.npz",
        "panels.parquet",
        "record.json",
        "sample.log",
    }


def test_resume_rejects_checksum_mismatch_and_state_only_remains_diagnostic(tmp_path):
    config = _mini_config(tmp_path)
    root = tmp_path / "run"
    run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        overwrite=True,
        generation_only=True,
        log_level="ERROR",
    )
    diagnostic = run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        resume=True,
        state_only=True,
        scenario="clean",
        max_dates=3,
        log_level="ERROR",
    )
    assert diagnostic["status"] == "estimating"
    assert "clean" in diagnostic["diagnostics"]["state_only"]
    assert diagnostic["estimation"] == {}

    panel_path = root / "sample_000" / "panels.parquet"
    with panel_path.open("ab") as file_handle:
        file_handle.write(b"corruption")
    with pytest.raises(ValueError, match="checksum"):
        run_sample(
            config=config,
            sample_id=0,
            output_root=root,
            resume=True,
            generation_only=True,
            log_level="ERROR",
        )


def _reference_rows(config, sample_id):
    simulation = replace(config.simulation, seed=config.base_seed + sample_id)
    path = HestonPathSimulator(
        params=config.dgp,
        config=simulation,
        variance_drawer=AndersenQeVarianceDrawer(),
    ).simulate()
    clean = generate_clean_option_panel_rows(
        run_id=config.run_id,
        sample_id=sample_id,
        path=path,
        params_p=config.dgp,
        eta_v=config.eta_v,
        maturities_years=config.cos_basis.maturities,
        log_moneyness=config.log_moneyness,
        atm_option_type=config.atm_option_type,
        pricing_method="COS",
        iv_method="lets_be_rational",
        pricer=CosOptionPricer(),
        cos_basis=config.cos_basis,
    )
    scenarios = {"clean": clean}
    factors = {}
    assert config.noise is not None
    for scenario in NOISE_SCENARIOS:
        rows, factor_rows = generate_noisy_panel_rows(
            clean,
            scenario=scenario,
            seed=scenario_seed(config.noise.base_seed, sample_id, scenario),
            config=config.noise,
        )
        scenarios[scenario] = rows
        factors[scenario] = factor_rows
    return path, scenarios, factors


def test_short_profile_path_panels_factors_and_estimates_are_numerically_equivalent(tmp_path):
    config = _mini_config(tmp_path)
    reference_path, reference_rows, reference_factors = _reference_rows(config, 0)
    root = tmp_path / "run"
    run_sample(
        config=config,
        sample_id=0,
        output_root=root,
        overwrite=True,
        generation_only=True,
        log_level="ERROR",
    )

    generated_path, _, _ = load_heston_path_npz(root / "sample_000" / "path.npz")
    for name in ("t_week", "logS_week", "V_week", "dlogS_week", "logS_daily", "V_daily"):
        np.testing.assert_array_equal(getattr(generated_path, name), getattr(reference_path, name))

    combined_path = root / "sample_000" / "panels.parquet"
    combined = pd.read_parquet(combined_path)
    numeric_columns = [
        "t",
        "S",
        "logS",
        "V",
        "maturity_years",
        "strike",
        "forward",
        "model_price",
        "model_iv",
    ]
    noise_numeric = [
        "noise_draw",
        "raw_noisy_iv",
        "raw_price_before_rounding",
        "price_after_rounding",
        "observed_price",
        "observed_iv",
    ]
    legacy_paths = {}
    metadata = {
        "sample_id": 0,
        "cos_basis": cos_specification_metadata(config.cos_basis),
    }
    for scenario in SCENARIO_ORDER:
        selected = combined[combined["scenario"] == scenario].reset_index(drop=True)
        reference = pd.DataFrame(reference_rows[scenario]).reset_index(drop=True)
        for column in numeric_columns:
            np.testing.assert_array_equal(selected[column].to_numpy(), reference[column].to_numpy())
        if scenario == "clean":
            np.testing.assert_array_equal(
                selected["estimation_iv"].to_numpy(),
                reference["model_iv"].to_numpy(),
            )
        else:
            for column in noise_numeric:
                np.testing.assert_array_equal(
                    selected[column].to_numpy(),
                    reference[column].to_numpy(),
                )
            np.testing.assert_array_equal(
                selected["estimation_iv"].to_numpy(),
                reference["observed_iv"].to_numpy(),
            )
        legacy_paths[scenario] = write_records(
            reference_rows[scenario],
            tmp_path / "legacy" / scenario,
            panel_format="parquet",
            metadata={**metadata, "scenario": scenario},
        )

    factor_frame = combined[combined["scenario"] == "persistent_factor"]
    factor_once = (
        factor_frame[
            [
                "week_index",
                "persistent_factor_level",
                "persistent_factor_moneyness",
                "persistent_factor_maturity",
            ]
        ]
        .drop_duplicates("week_index")
        .reset_index(drop=True)
    )
    expected_factors = pd.DataFrame(reference_factors["persistent_factor"])
    np.testing.assert_array_equal(
        factor_once["persistent_factor_level"].to_numpy(),
        expected_factors["factor_0"].to_numpy(),
    )
    np.testing.assert_array_equal(
        factor_once["persistent_factor_moneyness"].to_numpy(),
        expected_factors["factor_1"].to_numpy(),
    )
    np.testing.assert_array_equal(
        factor_once["persistent_factor_maturity"].to_numpy(),
        expected_factors["factor_2"].to_numpy(),
    )

    for scenario in SCENARIO_ORDER:
        legacy_panel = load_option_panel(legacy_paths[scenario])
        combined_panel = load_option_panel(combined_path, scenario=scenario)
        for legacy_date, combined_date in zip(legacy_panel.dates, combined_panel.dates):
            np.testing.assert_array_equal(legacy_date.observed_iv, combined_date.observed_iv)
            np.testing.assert_array_equal(legacy_date.observed_price, combined_date.observed_price)
            assert legacy_date.true_variance == combined_date.true_variance
        legacy_estimate = estimate_first_step(
            legacy_panel,
            criterion_config=config.criterion_config,
            powell_config=config.powell_config,
        )
        combined_estimate = estimate_first_step(
            combined_panel,
            criterion_config=config.criterion_config,
            powell_config=config.powell_config,
        )
        assert combined_estimate.final_criterion == pytest.approx(
            legacy_estimate.final_criterion,
            rel=0.0,
            abs=0.0,
        )
        assert asdict(combined_estimate.estimated_parameters) == pytest.approx(
            asdict(legacy_estimate.estimated_parameters),
            rel=0.0,
            abs=0.0,
        )
        np.testing.assert_array_equal(
            combined_estimate.free_parameters,
            legacy_estimate.free_parameters,
        )
        assert [item.criterion for item in combined_estimate.powell_passes] == pytest.approx(
            [item.criterion for item in legacy_estimate.powell_passes],
            rel=0.0,
            abs=0.0,
        )
