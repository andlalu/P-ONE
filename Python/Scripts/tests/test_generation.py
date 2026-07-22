import csv
import json
from pathlib import Path

import pytest

from DGPSimulation.io import load_heston_path_npz
from OptionData.io import load_option_panel, panel_metadata_path
from OptionPricing.clean_panel import write_panel
from OptionPricing.noisy_panel import read_table
from Scripts.experiment_config import load_experiment_config
from Scripts.generation import run_samples, with_overrides, write_run_metadata
from Scripts.validate_generation import validate_run

PRODUCTION_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "heston_experiment_run_001.json"


def _mini_config(tmp_path, *, noise=True, skip_existing=False):
    payload = json.loads(PRODUCTION_CONFIG.read_text())
    payload["run"].update(
        run_id="test_run",
        n_samples=1,
        base_seed=9000,
        output_root="generation_test",
        panel_format="csv",
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
    if not noise:
        payload["noise"] = None
    config_path = tmp_path / "experiment.json"
    config_path.write_text(json.dumps(payload))
    return with_overrides(load_experiment_config(config_path), skip_existing=skip_existing)


def test_clean_and_noisy_generation_workflow_and_skip_existing(tmp_path):
    config = _mini_config(tmp_path)
    results = run_samples(config=config)
    write_run_metadata(config, results)

    assert [result.status for result in results] == ["ok"]
    result = results[0]
    panel = load_option_panel(result.panel_file)
    assert panel.metadata["cos_basis"]["effective_widths"] == [1.5]
    assert panel.metadata["cos_basis"]["generation_n_cos"] == 64
    assert panel.metadata["experiment_config_hash"] == config.experiment_config_hash
    assert len(read_table(result.panel_file)) == 15
    assert {item.noise_scenario for item in result.noisy_results} == {
        "low_iid",
        "spatial_corr",
        "persistent_factor",
    }
    assert all(item.n_rows == 15 and item.status == "ok" for item in result.noisy_results)
    assert (Path(config.output_root) / "noise_factors" / "persistent_factor" / "sample_000.csv").exists()
    validate_run(run_root=config.output_root, expected_samples=1)

    manifest = Path(config.output_root) / "config" / "manifest_generation.csv"
    with manifest.open(newline="") as file_handle:
        assert len(list(csv.DictReader(file_handle))) == 1
    assert (Path(config.output_root) / "config" / "manifest_noisy_panels.csv").exists()

    skipped = run_samples(config=_mini_config(tmp_path, skip_existing=True))
    assert [item.status for item in skipped] == ["skipped"]


def test_null_noise_and_relative_output_path(tmp_path):
    config = _mini_config(tmp_path, noise=False)
    assert config.noise is None
    assert config.output_root == str((tmp_path / "generation_test").resolve())
    result = run_samples(config=config)[0]
    assert result.noisy_results == ()
    path, _, stored_config = load_heston_path_npz(result.path_file)
    assert stored_config.return_daily is True
    assert path.logS_daily is not None and path.V_daily is not None


def test_run_metadata_and_completeness_aware_skip(tmp_path):
    config = _mini_config(tmp_path)
    result = run_samples(config=config)[0]
    write_run_metadata(config, [result])

    run_metadata = json.loads((Path(config.output_root) / "config" / "run_metadata.json").read_text())
    assert run_metadata["git_sha"] != ""
    assert run_metadata["requested_workers"] == 1
    assert run_metadata["resolved_workers"] == 1
    assert run_metadata["experiment_config_hash"] == config.experiment_config_hash
    assert run_metadata["sample_range"] == [0, 1]
    copied = json.loads((Path(config.output_root) / "config" / "experiment_config.json").read_text())
    assert copied == config.raw_config
    assert not list(Path(config.output_root).rglob("*.tmp*"))

    panel_metadata_path(result.panel_file).unlink()
    regenerated = run_samples(config=_mini_config(tmp_path, skip_existing=True))[0]
    assert regenerated.status == "ok"


def test_requested_parquet_format_does_not_fall_back_to_csv(tmp_path, monkeypatch):
    monkeypatch.setattr("OptionPricing.clean_panel.parquet_available", lambda: False)
    with pytest.raises(RuntimeError, match="requires pandas"):
        write_panel([], tmp_path / "panel", metadata={}, panel_format="parquet")
    assert not (tmp_path / "panel.csv").exists()
