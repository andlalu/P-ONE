import numpy as np
import pytest
from pathlib import Path

from Models.Heston.parameters import HestonRiskNeutralParameters
from OptionPricing.cos_basis import FixedCosBasisConfig
from OptionPricing.cos_pricer import CosOptionPricer


def test_fixed_cos_basis_validation_and_unambiguous_matching():
    config = FixedCosBasisConfig((0.25, 0.5), (1.2, 1.8), 32, 16, maturity_tolerance=1e-8)
    config.validate()
    assert config.width_for_maturity(0.25 + 1e-9) == 1.2
    with pytest.raises(ValueError, match="not present"):
        config.width_for_maturity(1.0)
    with pytest.raises(ValueError, match="unique"):
        FixedCosBasisConfig((0.25, 0.25 + 1e-11), (1.2, 1.3), 32, 16).validate()
    with pytest.raises(ValueError, match="one-to-one"):
        FixedCosBasisConfig((0.25,), (1.2, 1.3), 32, 16).validate()


def test_fixed_basis_grid_is_independent_of_candidate_variance(monkeypatch):
    pricer = CosOptionPricer()
    params = HestonRiskNeutralParameters(kappa=2.0, vbar=0.04, sigma_v=0.4, rho=-0.5)
    basis = pricer.prepare_fixed_basis(maturity=0.25, effective_width=1.5, n_cos=48, model_params=params)
    u_before = basis.u_grid.copy()

    def dynamic_path_forbidden(*args, **kwargs):
        raise AssertionError("variance-dependent width path was invoked")

    monkeypatch.setattr(pricer, "variance_scaled_effective_width", dynamic_path_forbidden)
    prices = pricer.price_matrix_fixed_basis(
        log_s=np.log(np.array([100.0, 100.0])),
        variance=np.array([0.01, 0.09]),
        strike_grid=np.array([100.0]),
        rate=0.0,
        dividend_yield=0.0,
        basis=basis,
    )
    assert prices.shape == (2, 1, 1)
    np.testing.assert_array_equal(basis.u_grid, u_before)


def test_cos_basis_metadata_roundtrip_and_mismatch_rejection():
    config = FixedCosBasisConfig((0.25, 0.5), (1.2, 1.8), 64, 32)
    restored = FixedCosBasisConfig.from_dict(config.to_dict())
    assert restored == config
    config.assert_panel_compatible(restored.generation_metadata())
    with pytest.raises(ValueError, match="width mismatch"):
        FixedCosBasisConfig((0.25, 0.5), (1.2, 1.9), 64, 32).assert_panel_compatible(
            config.generation_metadata()
        )


def test_generation_and_estimation_term_counts_are_role_specific_but_compatible():
    panel_basis = FixedCosBasisConfig((0.25,), (1.5,), 128, 64)
    estimator_basis = FixedCosBasisConfig((0.25,), (1.5,), 128, 32)
    with pytest.raises(ValueError, match="configuration hash"):
        estimator_basis.assert_panel_compatible(panel_basis.generation_metadata())

    shared_basis = FixedCosBasisConfig((0.25,), (1.5,), 128, 32)
    shared_basis.assert_panel_compatible(shared_basis.generation_metadata())


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("basis_convention", "variance_scaled", "conventions differ"),
        ("pricing_implementation", "another_version", "implementations differ"),
        ("maturities", [0.25, 1.0], "maturity grids differ"),
        ("effective_widths", [1.2, 1.9], "width mismatch"),
    ],
)
def test_panel_basis_metadata_rejects_exact_architecture_mismatches(field, value, message):
    config = FixedCosBasisConfig((0.25, 0.5), (1.2, 1.8), 64, 32)
    metadata = config.generation_metadata()
    metadata[field] = value
    with pytest.raises(ValueError, match=message):
        config.assert_panel_compatible(metadata)


def test_production_basis_load_and_hash_are_stable():
    basis_path = (
        Path(__file__).resolve().parents[2] / "Scripts" / "configs" / "heston_cos_basis_production.json"
    )
    config = FixedCosBasisConfig.load(basis_path)
    assert config.stable_hash() == FixedCosBasisConfig.from_dict(config.to_dict()).stable_hash()
    assert config.generation_n_cos >= config.estimation_n_cos
