from pathlib import Path


PYTHON_ROOT = Path(__file__).resolve().parents[2]


def _production_sources() -> dict[Path, str]:
    return {
        path: path.read_text()
        for path in PYTHON_ROOT.rglob("*.py")
        if "tests" not in path.parts and not path.name.startswith(".")
    }


def test_removed_interfaces_and_legacy_cos_paths_have_no_production_references():
    sources = _production_sources()
    forbidden = (
        "DGPSimulation.base",
        "OptionPricing.base",
        "VarianceScaledCosConfig",
        "variance_scaled_effective_width",
        "price_matrix_variance_scaled_reference",
        "price_one_variance_scaled_reference",
    )
    for path, source in sources.items():
        for name in forbidden:
            assert name not in source, f"{name} remains in {path}"


def test_panel_persistence_is_owned_by_option_data_io():
    clean_source = (PYTHON_ROOT / "OptionData" / "clean_panel.py").read_text()
    noise_sources = [
        (PYTHON_ROOT / "OptionData" / name).read_text()
        for name in (
            "noise_common.py",
            "noise_low_iid.py",
            "noise_spatial.py",
            "noise_persistent_factor.py",
            "add_noise.py",
        )
    ]
    assert "write_panel" not in clean_source
    assert "heston_q_from_p" not in clean_source
    assert all("csv." not in source and "to_parquet" not in source for source in noise_sources)
