import numpy as np

from Estimation.ISCGMM.heston_transition_cf import heston_p_transition_cf
from Estimation.ISCGMM.types import HestonEstimationParams


def test_heston_transition_cf_zero_and_conjugacy():
    theta = HestonEstimationParams(
        eta=5.0,
        kappa=7.0,
        vbar=0.0225,
        sigma_v=0.4,
        rho=-0.5,
        eta_v=5.0,
        r=0.02,
        q=0.0,
    )
    x_prev = np.array([[0.01, 0.02], [-0.02, 0.04]], dtype=float)
    nodes = np.array([[0.0, 0.0], [0.4, -0.25], [-0.4, 0.25]], dtype=float)
    values = heston_p_transition_cf(nodes, x_prev, dt=5.0 / 252.0, theta=theta, rk_steps=24)
    assert values.shape == (3, 2)
    np.testing.assert_allclose(values[0], np.ones(2), atol=1e-12)
    np.testing.assert_allclose(values[2], np.conjugate(values[1]), rtol=1e-11, atol=1e-11)
    assert np.all(np.isfinite(values.real))
    assert np.all(np.isfinite(values.imag))
