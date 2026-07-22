import numpy as np
import pytest

from Models.Heston.parameters import HestonRiskNeutralParameters
from OptionPricing.cos_basis import (
    FixedCosBasisConfig,
    cos_specification_metadata,
    validate_panel_cos_compatibility,
)
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


def test_panel_compatibility_uses_generation_settings_only():
    panel_basis = FixedCosBasisConfig((0.25,), (1.5,), 128, 64)
    estimator_basis = FixedCosBasisConfig((0.25,), (1.5,), 128, 32)
    validate_panel_cos_compatibility(estimator_basis, cos_specification_metadata(panel_basis))

    mismatched_terms = FixedCosBasisConfig((0.25,), (1.5,), 64, 32)
    with pytest.raises(ValueError, match="generation_n_cos"):
        validate_panel_cos_compatibility(mismatched_terms, cos_specification_metadata(panel_basis))


def test_panel_compatibility_rejects_maturity_and_width_mismatches():
    config = FixedCosBasisConfig((0.25, 0.5), (1.2, 1.8), 64, 32)
    metadata = cos_specification_metadata(config)
    metadata["effective_widths"] = [1.2, 1.9]
    with pytest.raises(ValueError, match="width mismatch"):
        validate_panel_cos_compatibility(config, metadata)

    metadata = cos_specification_metadata(config)
    metadata["maturities"] = [0.25, 1.0]
    with pytest.raises(ValueError, match="maturity grids differ"):
        validate_panel_cos_compatibility(config, metadata)
