import math

import numpy as np

from sim.types import HestonParamsP
from sim.variance_drawers import AndersenQeVarianceDrawer


def test_qe_stepper_matches_conditional_moments():
    params = HestonParamsP(eta=1.2, kappa=2.2, vbar=0.05, sigma_v=0.5, rho=-0.4, r=0.0, q=0.0)
    drawer = AndersenQeVarianceDrawer()

    delta = 1.0 / 252.0
    v_n = 0.06

    e = math.exp(-params.kappa * delta)
    m = params.vbar + (v_n - params.vbar) * e
    s2 = (
        v_n * (params.sigma_v ** 2) * e * (1.0 - e) / params.kappa
        + params.vbar * (params.sigma_v ** 2) * ((1.0 - e) ** 2) / (2.0 * params.kappa)
    )

    rng = np.random.default_rng(123)
    draws = np.array([drawer.draw_next_variance(v_n, delta, rng, params) for _ in range(200_000)])

    sample_mean = float(np.mean(draws))
    sample_var = float(np.var(draws))

    assert abs(sample_mean - m) < 3.0e-4
    assert abs(sample_var - s2) < 8.0e-5
