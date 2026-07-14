import csv
import json
import pytest

from DGPSimulation.types import HestonSimConfig
from Models.Heston.parameters import HestonPhysicalParameters
from OptionData.io import load_option_panel
from OptionPricing.cos_basis import FixedCosBasisConfig
from OptionPricing.clean_panel import write_panel
from DGPSimulation.io import load_heston_path_npz
from OptionData.io import panel_metadata_path
from Scripts.generation import (
    GenerationConfig,
    PanelConfig,
    load_config,
    resolved_config_dict,
    run_samples,
    write_run_metadata,
)
from Scripts.validate_generation import validate_run


def _mini_config(tmp_path, skip_existing=False):
    return GenerationConfig(
        run_id="test_run",
        n_samples=2,
        base_seed=9000,
        output_root=str(tmp_path / "generation_test"),
        dgp=HestonPhysicalParameters(eta=1.5, kappa=3.0, vbar=0.04, sigma_v=0.4, rho=-0.7, r=0.02, q=0.0),
        simulation=HestonSimConfig(t_week=4, burnin_days=2, return_daily=True),
        eta_v=0.0,
        panel=PanelConfig(
            maturities_years=(0.25,),
            log_moneyness=(-0.05, 0.0, 0.05),
            option_type_rule="otm_by_forward_moneyness",
            atm_option_type="call",
            pricing_method="COS",
            iv_method="lets_be_rational",
        ),
        panel_format="csv",
        cos_basis=FixedCosBasisConfig((0.25,), (1.5,), 64, 32),
        cos_basis_path=str(tmp_path / "test_basis.json"),
        configuration_hash="test-configuration-hash",
        workers=2,
        skip_existing=skip_existing,
    )


def test_clean_generation_workflow_and_skip_existing(tmp_path):
    config = _mini_config(tmp_path)
    results = run_samples(config=config)
    write_run_metadata(config, results)

    assert {result.status for result in results} == {"ok"}
    panel = load_option_panel(results[0].panel_file)
    assert panel.metadata["cos_basis"]["effective_widths"] == [1.5]
    assert panel.metadata["cos_basis"]["generation_n_cos"] == 64
    validate_run(run_root=config.output_root, expected_samples=2)

    manifest = tmp_path / "generation_test" / "config" / "manifest_generation.csv"
    with manifest.open(newline="") as fh:
        assert len(list(csv.DictReader(fh))) == 2

    skipped = run_samples(config=_mini_config(tmp_path, skip_existing=True))
    assert {result.status for result in skipped} == {"skipped"}


def test_null_noise_and_relative_paths_are_explicit(tmp_path):
    basis = FixedCosBasisConfig((0.25,), (1.5,), 64, 32)
    (tmp_path / "basis.json").write_text(json.dumps(basis.to_dict()))
    payload = {
        "run_id": "null_noise",
        "n_samples": 1,
        "base_seed": 10,
        "output_root": "relative-output",
        "panel_format": "csv",
        "cos_basis_path": "basis.json",
        "dgp": {"eta": 1.5, "kappa": 3.0, "vbar": 0.04, "sigma_v": 0.4, "rho": -0.7, "r": 0.02, "q": 0.0},
        "simulation": {"t_week": 2, "burnin_days": 0, "return_daily": False},
        "q_measure": {"eta_v": 0.0},
        "panel": {
            "maturities_years": [0.25],
            "log_moneyness": [-0.05, 0.0, 0.05],
            "option_type_rule": "otm_by_forward_moneyness",
            "atm_option_type": "call",
            "pricing_method": "COS",
            "iv_method": "lets_be_rational"
        },
        "parallel": {"workers": 1},
        "noise": None,
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(payload))
    config = load_config(config_path)
    assert config.noise is None
    assert config.output_root == str((tmp_path / "relative-output").resolve())
    results = run_samples(config=config)
    assert results[0].noisy_results == ()
    path, _, stored_config = load_heston_path_npz(results[0].path_file)
    assert stored_config.return_daily is False
    assert path.logS_daily is None and path.V_daily is None


def test_execution_metadata_atomic_outputs_and_completeness_aware_skip(tmp_path):
    config = _mini_config(tmp_path)
    result = run_samples(config=config, sample_end=1)[0]
    write_run_metadata(config, [result])
    execution = resolved_config_dict(config, requested_workers=2, resolved_workers=1)["execution"]
    assert execution["git_sha"] != ""
    assert execution["requested_worker_count"] == 2
    assert execution["resolved_worker_count"] == 1
    assert not list((tmp_path / "generation_test").rglob("*.tmp*"))

    panel_metadata_path(result.panel_file).unlink()
    regenerated = run_samples(config=_mini_config(tmp_path, skip_existing=True), sample_end=1)[0]
    assert regenerated.status == "ok"


def test_requested_parquet_format_does_not_fall_back_to_csv(tmp_path, monkeypatch):
    monkeypatch.setattr("OptionPricing.clean_panel.parquet_available", lambda: False)
    with pytest.raises(RuntimeError, match="requires pandas"):
        write_panel([], tmp_path / "panel", metadata={}, panel_format="parquet")
    assert not (tmp_path / "panel.csv").exists()
