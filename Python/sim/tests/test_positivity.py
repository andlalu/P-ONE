import numpy as np

from sim.heston_simulator import HestonPathSimulator
from sim.types import HestonParamsP, HestonSimConfig
from sim.variance_drawers import AndersenQeVarianceDrawer, EulerVarianceDrawer


STRESS_PARAM_SETS = [
    HestonParamsP(eta=1.0, kappa=1.5, vbar=0.04, sigma_v=0.7, rho=-0.8, r=0.01, q=0.0),
    HestonParamsP(eta=1.0, kappa=0.3, vbar=0.02, sigma_v=1.0, rho=0.85, r=0.01, q=0.0),
]


def test_variance_nonnegative_qe_stress_cases():
    config = HestonSimConfig(seed=100, return_daily=True)

    for params in STRESS_PARAM_SETS:
        path = HestonPathSimulator(params, config, AndersenQeVarianceDrawer()).simulate()
        assert path.V_daily is not None
        assert np.all(path.V_daily >= 0.0)


def test_variance_nonnegative_euler_stress_cases():
    config = HestonSimConfig(seed=100, return_daily=True)

    for params in STRESS_PARAM_SETS:
        path = HestonPathSimulator(params, config, EulerVarianceDrawer()).simulate()
        assert path.V_daily is not None
        assert np.all(path.V_daily >= 0.0)
