import csv

from DGPSimulation.types import HestonSimConfig
from Models.Heston.parameters import HestonPhysicalParameters
from OptionData.io import load_option_panel
from Scripts.generation import GenerationConfig, PanelConfig, run_samples, write_run_metadata
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
            cos_n_terms=64,
            cos_effective_widths=(1.5,),
        ),
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
    validate_run(run_root=config.output_root, expected_samples=2)

    manifest = tmp_path / "generation_test" / "config" / "manifest_generation.csv"
    with manifest.open(newline="") as fh:
        assert len(list(csv.DictReader(fh))) == 2

    skipped = run_samples(config=_mini_config(tmp_path, skip_existing=True))
    assert {result.status for result in skipped} == {"skipped"}
