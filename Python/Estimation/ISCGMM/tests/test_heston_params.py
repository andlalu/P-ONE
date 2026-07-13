import numpy as np

from Estimation.ISCGMM.parameter_transform import from_free, to_free
from Models.Heston.parameters import HestonParameters


def test_parameter_transform_roundtrip_and_constraints():
    theta = HestonParameters(
        eta=5.0,
        kappa=7.0,
        vbar=0.0225,
        sigma_v=0.4,
        rho=-0.5,
        eta_v=5.0,
        r=0.02,
        q=0.0,
    )
    free = to_free(theta)
    recovered = from_free(free, r=theta.r, q=theta.q)
    np.testing.assert_allclose(
        [recovered.eta, recovered.kappa, recovered.vbar, recovered.sigma_v, recovered.rho, recovered.eta_v],
        [theta.eta, theta.kappa, theta.vbar, theta.sigma_v, theta.rho, theta.eta_v],
        rtol=1e-12,
        atol=1e-12,
    )
    assert recovered.kappa > 0.0
    assert recovered.vbar > 0.0
    assert recovered.sigma_v > 0.0
    assert -1.0 < recovered.rho < 1.0
    assert recovered.kappa - recovered.eta_v > 0.0
