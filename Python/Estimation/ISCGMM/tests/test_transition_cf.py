import numpy as np
import pytest

from Estimation.ISCGMM.heston_transition_cf import heston_p_transition_cf, heston_p_transform_coefficients
from Models.Heston.parameters import HestonParameters


def _theta() -> HestonParameters:
    return HestonParameters(
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
    values = heston_p_transition_cf(nodes, x_prev, dt=5.0 / 252.0, theta=theta, rk_steps=24, method="analytic")
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
    for method in ("analytic", "rk4"):
        values = heston_p_transition_cf(nodes, x_prev, dt=5.0 / 252.0, theta=theta, rk_steps=32, method=method)
        assert values.shape == (nodes.shape[0], x_prev.shape[0])
        assert np.all(np.isfinite(values.real))
        assert np.all(np.isfinite(values.imag))
        assert np.max(np.abs(values)) <= 1.0 + 1e-12


def test_heston_transition_cf_analytic_matches_rk4_reference():
    theta = _theta()
    x_prev = np.array([[0.0, 0.005], [0.02, 0.025], [-0.03, 0.09]], dtype=float)
    nodes = np.array(
        [
            [0.0, 0.0],
            [0.25, 0.10],
            [1.0, -1.5],
            [-1.25, 0.8],
            [2.0, 0.5],
        ],
        dtype=float,
    )

    analytic = heston_p_transition_cf(nodes, x_prev, dt=5.0 / 252.0, theta=theta, method="analytic")
    rk4 = heston_p_transition_cf(nodes, x_prev, dt=5.0 / 252.0, theta=theta, rk_steps=512, method="rk4")
    np.testing.assert_allclose(analytic, rk4, rtol=5e-10, atol=5e-12)


def test_heston_transition_cf_first_moment_derivatives():
    theta = _theta()
    dt = 5.0 / 252.0
    v0 = np.array([0.008, 0.024, 0.07], dtype=float)
    x_prev = np.column_stack((np.zeros_like(v0), v0))
    eps = 1e-5

    phi_u_plus = heston_p_transition_cf(np.array([[eps, 0.0]]), x_prev, dt=dt, theta=theta, method="analytic")[0]
    phi_u_minus = heston_p_transition_cf(np.array([[-eps, 0.0]]), x_prev, dt=dt, theta=theta, method="analytic")[0]
    mean_return = ((phi_u_plus - phi_u_minus) / (2.0 * eps) / 1j).real

    phi_v_plus = heston_p_transition_cf(np.array([[0.0, eps]]), x_prev, dt=dt, theta=theta, method="analytic")[0]
    phi_v_minus = heston_p_transition_cf(np.array([[0.0, -eps]]), x_prev, dt=dt, theta=theta, method="analytic")[0]
    mean_variance = ((phi_v_plus - phi_v_minus) / (2.0 * eps) / 1j).real

    expected_variance = theta.vbar + (v0 - theta.vbar) * np.exp(-theta.kappa * dt)
    expected_integrated_variance = theta.vbar * dt + (v0 - theta.vbar) * (1.0 - np.exp(-theta.kappa * dt)) / theta.kappa
    expected_return = (theta.r - theta.q) * dt + (theta.eta - 0.5) * expected_integrated_variance

    np.testing.assert_allclose(mean_variance, expected_variance, rtol=2e-8, atol=2e-10)
    np.testing.assert_allclose(mean_return, expected_return, rtol=2e-8, atol=2e-10)


def test_heston_transition_coefficient_methods_match_at_zero_node():
    theta = _theta()
    nodes = np.array([[0.0, 0.0], [0.2, -0.3]], dtype=float)
    a_analytic, b_analytic = heston_p_transform_coefficients(nodes, dt=5.0 / 252.0, theta=theta, method="analytic")
    a_rk4, b_rk4 = heston_p_transform_coefficients(nodes, dt=5.0 / 252.0, theta=theta, method="rk4", rk_steps=512)
    np.testing.assert_allclose(a_analytic, a_rk4, rtol=5e-10, atol=5e-12)
    np.testing.assert_allclose(b_analytic, b_rk4, rtol=5e-10, atol=5e-12)


@pytest.mark.slow
def test_heston_transition_cf_matches_small_monte_carlo_check():
    theta = _theta()
    dt = 5.0 / 252.0
    v0 = 0.024
    node = np.array([[0.35, -0.20]], dtype=float)
    analytic = heston_p_transition_cf(node, np.array([[0.0, v0]], dtype=float), dt=dt, theta=theta, method="analytic")[0, 0]

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
