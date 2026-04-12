import numpy as np

from sim.heston_simulator import HestonPathSimulator
from sim.io import load_heston_path_npz, save_heston_path_npz
from sim.types import HestonParamsP, HestonSimConfig
from sim.variance_drawers import AndersenQeVarianceDrawer


def _params() -> HestonParamsP:
    return HestonParamsP(eta=1.1, kappa=1.8, vbar=0.03, sigma_v=0.4, rho=-0.5, r=0.01, q=0.0)


def test_npz_roundtrip_with_daily_series(tmp_path):
    params = _params()
    config = HestonSimConfig(seed=202, return_daily=True, t_week=32, burnin_days=25)
    sim = HestonPathSimulator(params=params, config=config, variance_drawer=AndersenQeVarianceDrawer())
    path = sim.simulate()

    out_file = tmp_path / "heston_path_daily.npz"
    save_heston_path_npz(out_file, path=path, params=params, config=config)

    loaded_path, loaded_params, loaded_config = load_heston_path_npz(out_file)

    assert loaded_path.seed == path.seed
    np.testing.assert_allclose(loaded_path.t_week, path.t_week)
    np.testing.assert_allclose(loaded_path.logS_week, path.logS_week)
    np.testing.assert_allclose(loaded_path.V_week, path.V_week)
    np.testing.assert_allclose(loaded_path.dlogS_week, path.dlogS_week)
    np.testing.assert_allclose(loaded_path.logS_daily, path.logS_daily)
    np.testing.assert_allclose(loaded_path.V_daily, path.V_daily)
    assert loaded_params == params
    assert loaded_config == config


def test_npz_roundtrip_without_daily_series(tmp_path):
    params = _params()
    config = HestonSimConfig(seed=303, return_daily=False, t_week=16, burnin_days=10)
    sim = HestonPathSimulator(params=params, config=config, variance_drawer=AndersenQeVarianceDrawer())
    path = sim.simulate()

    out_file = tmp_path / "heston_path_weekly_only.npz"
    save_heston_path_npz(out_file, path=path, params=params, config=config, compressed=False)

    loaded_path, loaded_params, loaded_config = load_heston_path_npz(out_file)

    assert loaded_path.seed == config.seed
    assert loaded_path.logS_daily is None
    assert loaded_path.V_daily is None
    np.testing.assert_allclose(loaded_path.logS_week, path.logS_week)
    np.testing.assert_allclose(loaded_path.V_week, path.V_week)
    assert loaded_params == params
    assert loaded_config == config
