import numpy as np

from sim.heston_simulator import HestonPathSimulator
from sim.types import HestonParamsP, HestonSimConfig
from sim.variance_drawers import AndersenQeVarianceDrawer, EulerVarianceDrawer


def _base_params() -> HestonParamsP:
    return HestonParamsP(eta=1.2, kappa=2.0, vbar=0.04, sigma_v=0.5, rho=-0.6, r=0.01, q=0.0)


def test_determinism_qe_drawer():
    params = _base_params()
    config = HestonSimConfig(seed=42)

    sim_1 = HestonPathSimulator(params, config, AndersenQeVarianceDrawer())
    sim_2 = HestonPathSimulator(params, config, AndersenQeVarianceDrawer())

    path_1 = sim_1.simulate()
    path_2 = sim_2.simulate()

    np.testing.assert_allclose(path_1.logS_week, path_2.logS_week)
    np.testing.assert_allclose(path_1.V_week, path_2.V_week)


def test_determinism_euler_drawer():
    params = _base_params()
    config = HestonSimConfig(seed=42)

    sim_1 = HestonPathSimulator(params, config, EulerVarianceDrawer())
    sim_2 = HestonPathSimulator(params, config, EulerVarianceDrawer())

    path_1 = sim_1.simulate()
    path_2 = sim_2.simulate()

    np.testing.assert_allclose(path_1.logS_week, path_2.logS_week)
    np.testing.assert_allclose(path_1.V_week, path_2.V_week)
