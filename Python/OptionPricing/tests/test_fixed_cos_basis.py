import pytest

from OptionPricing.config import FixedCosBasisConfig
from OptionPricing.cos_basis import cos_specification_metadata, validate_panel_cos_compatibility


def test_valid_fixed_cos_basis_lookup():
    config = FixedCosBasisConfig((0.25, 0.5), (1.2, 1.8), 32, 16, maturity_tolerance=1e-8)
    config.validate()
    assert config.width_for_maturity(0.25 + 1e-9) == 1.2


def test_unknown_fixed_cos_maturity_is_rejected():
    config = FixedCosBasisConfig((0.25, 0.5), (1.2, 1.8), 32, 16)
    with pytest.raises(ValueError, match="not present"):
        config.width_for_maturity(1.0)


def test_mismatched_fixed_cos_configuration_is_rejected():
    with pytest.raises(ValueError, match="unique"):
        FixedCosBasisConfig((0.25, 0.25 + 1e-11), (1.2, 1.3), 32, 16).validate()
    with pytest.raises(ValueError, match="one-to-one"):
        FixedCosBasisConfig((0.25,), (1.2, 1.3), 32, 16).validate()


def test_panel_compatibility_uses_generation_settings_only_and_rejects_mismatch():
    panel_basis = FixedCosBasisConfig((0.25,), (1.5,), 128, 64)
    estimator_basis = FixedCosBasisConfig((0.25,), (1.5,), 128, 32)
    validate_panel_cos_compatibility(estimator_basis, cos_specification_metadata(panel_basis))

    mismatched_terms = FixedCosBasisConfig((0.25,), (1.5,), 64, 32)
    with pytest.raises(ValueError, match="generation_n_cos"):
        validate_panel_cos_compatibility(mismatched_terms, cos_specification_metadata(panel_basis))
    metadata = cos_specification_metadata(panel_basis)
    metadata["effective_widths"] = [1.6]
    with pytest.raises(ValueError, match="width mismatch"):
        validate_panel_cos_compatibility(estimator_basis, metadata)
