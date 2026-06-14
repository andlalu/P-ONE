import numpy as np

from Estimation.ISCGMM.heston_transition_cf import heston_p_transition_cf
from Estimation.ISCGMM.types import HestonEstimationParams


def _theta() -> HestonEstimationParams:
    return HestonEstimationParams(
        eta=5.0,
        kappa=7.0,
        vbar=0.0225,
        sigma_v=0.4,
        rho=-0.5,
        eta_v=5.0,
        r=0.02,
        q=0.0,
    )


def test_heston_transition_cf_zero_and_conjugacy():
    theta = _theta()
    x_prev = np.array([[0.01, 0.02], [-0.02, 0.04]], dtype=float)
    nodes = np.array([[0.0, 0.0], [0.4, -0.25], [-0.4, 0.25]], dtype=float)
    values = heston_p_transition_cf(nodes, x_prev, dt=5.0 / 252.0, theta=theta, rk_steps=24)
    assert values.shape == (3, 2)
    np.testing.assert_allclose(values[0], np.ones(2), atol=1e-12)
    np.testing.assert_allclose(values[2], np.conjugate(values[1]), rtol=1e-11, atol=1e-11)
    assert np.all(np.isfinite(values.real))
    assert np.all(np.isfinite(values.imag))


def test_heston_transition_cf_weekly_dt_stability():
    theta = _theta()
    x_prev = np.array([[0.0, 0.005], [0.0, 0.025], [0.0, 0.09]], dtype=float)
    nodes = np.array(
        [
            [0.0, 0.0],
            [0.25, 0.10],
            [1.0, -1.5],
            [-1.0, 1.5],
            [2.0, 0.5],
        ],
        dtype=float,
    )
    values = heston_p_transition_cf(nodes, x_prev, dt=5.0 / 252.0, theta=theta, rk_steps=32)
    assert values.shape == (nodes.shape[0], x_prev.shape[0])
    assert np.all(np.isfinite(values.real))
    assert np.all(np.isfinite(values.imag))
    assert np.max(np.abs(values)) <= 1.0 + 1e-12


def test_heston_transition_cf_matches_small_monte_carlo_check():
    theta = _theta()
    dt = 5.0 / 252.0
    v0 = 0.024
    node = np.array([[0.35, -0.20]], dtype=float)
    analytic = heston_p_transition_cf(node, np.array([[0.0, v0]], dtype=float), dt=dt, theta=theta, rk_steps=64)[0, 0]

    n_paths = 60_000
    n_steps = 40
    step = dt / n_steps
    rng = np.random.default_rng(117)
    variance = np.full(n_paths, v0, dtype=float)
    log_return = np.zeros(n_paths, dtype=float)
    sqrt_one_minus_rho_sq = np.sqrt(1.0 - theta.rho * theta.rho)

    for _ in range(n_steps):
        variance_pos = np.maximum(variance, 0.0)
        z1 = rng.standard_normal(n_paths)
        z2 = rng.standard_normal(n_paths)
        dw1 = np.sqrt(step) * z1
        dw2 = np.sqrt(step) * z2
        log_return += (theta.r - theta.q + (theta.eta - 0.5) * variance_pos) * step
        log_return += np.sqrt(variance_pos) * dw1
        variance += theta.kappa * (theta.vbar - variance_pos) * step
        variance += theta.sigma_v * np.sqrt(variance_pos) * (theta.rho * dw1 + sqrt_one_minus_rho_sq * dw2)
        variance = np.maximum(variance, 0.0)

    sample = np.exp(1j * (node[0, 0] * log_return + node[0, 1] * variance))
    empirical = np.mean(sample)
    standard_error = np.std(sample, ddof=1) / np.sqrt(n_paths)
    assert abs(empirical - analytic) < 6.0 * standard_error + 2e-5
