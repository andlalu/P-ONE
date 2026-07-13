import pytest

from Models.Heston.parameters import HestonParameters, HestonPhysicalParameters, HestonRiskNeutralParameters


def test_joint_heston_parameter_projections_are_canonical():
    theta = HestonParameters(eta=5.0, kappa=7.0, vbar=0.0225, sigma_v=0.4, rho=-0.5, eta_v=5.0, r=0.02)
    physical = theta.to_physical()
    risk_neutral = theta.to_risk_neutral()
    assert physical == HestonPhysicalParameters(5.0, 7.0, 0.0225, 0.4, -0.5, 0.02, 0.0)
    assert risk_neutral.kappa == pytest.approx(2.0)
    assert risk_neutral.vbar == pytest.approx(7.0 * 0.0225 / 2.0)


@pytest.mark.parametrize(
    "parameters",
    [
        HestonPhysicalParameters(0.0, 1.0, 0.0, 0.4, -0.5),
        HestonRiskNeutralParameters(1.0, 0.04, 0.0, -0.5),
        HestonRiskNeutralParameters(1.0, 0.04, 0.4, 1.0),
    ],
)
def test_physical_and_risk_neutral_validation_rejects_invalid_values(parameters):
    with pytest.raises(ValueError):
        parameters.validate()
