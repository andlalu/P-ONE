import numpy as np

from sim.heston_simulator import HestonPathSimulator
from sim.types import HestonParamsP, HestonSimConfig
from sim.variance_drawers import AndersenQeVarianceDrawer


def _estimate_corr_for_rho(rho: float) -> float:
    corrs = []
    for seed in range(30):
        params = HestonParamsP(eta=1.0, kappa=2.0, vbar=0.04, sigma_v=0.45, rho=rho, r=0.0, q=0.0)
        config = HestonSimConfig(seed=seed, burnin_days=126, t_week=120, return_daily=True)
        path = HestonPathSimulator(params, config, AndersenQeVarianceDrawer()).simulate()

        dy = np.diff(path.logS_daily)
        dv = np.diff(path.V_daily)
        corr = np.corrcoef(dy, dv)[0, 1]
        corrs.append(corr)

    return float(np.mean(corrs))


def test_daily_dy_dv_correlation_tracks_rho_sign_and_monotonicity():
    corr_neg = _estimate_corr_for_rho(-0.8)
    corr_zero = _estimate_corr_for_rho(0.0)
    corr_pos = _estimate_corr_for_rho(0.8)

    assert corr_neg < corr_zero < corr_pos
    assert corr_neg < 0.0
    assert corr_pos > 0.0
