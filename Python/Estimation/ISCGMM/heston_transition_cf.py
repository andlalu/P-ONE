from __future__ import annotations

import numpy as np

from Estimation.ISCGMM.types import HestonEstimationParams


def heston_p_transform_coefficients(
    s_nodes: np.ndarray,
    *,
    dt: float,
    theta: HestonEstimationParams,
    rk_steps: int = 32,
) -> tuple[np.ndarray, np.ndarray]:
    """Affine coefficients for E[exp(i s1 r_T + i s2 V_T) | V_0].

    The estimation state is X_i = (r_i, V_i). The transition law only depends
    on V_i under the Heston P dynamics; r_i remains useful as an instrument.
    """

    theta.validate()
    nodes = np.asarray(s_nodes, dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("s_nodes must have shape (n_nodes, 2)")
    if dt < 0.0:
        raise ValueError("dt must be non-negative")
    if rk_steps <= 0:
        raise ValueError("rk_steps must be positive")

    u = nodes[:, 0]
    terminal_v = 1j * nodes[:, 1]
    a = np.zeros(nodes.shape[0], dtype=complex)
    b = terminal_v.astype(complex)
    if dt == 0.0:
        return a, b

    h = float(dt) / float(rk_steps)
    iu = 1j * u
    sigma = theta.sigma_v

    def rhs(a_val: np.ndarray, b_val: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        da = iu * (theta.r - theta.q) + theta.kappa * theta.vbar * b_val
        db = (
            iu * (theta.eta - 0.5)
            + (theta.rho * sigma * iu - theta.kappa) * b_val
            + 0.5 * sigma * sigma * b_val * b_val
            - 0.5 * u * u
        )
        return da, db

    for _ in range(rk_steps):
        k1a, k1b = rhs(a, b)
        k2a, k2b = rhs(a + 0.5 * h * k1a, b + 0.5 * h * k1b)
        k3a, k3b = rhs(a + 0.5 * h * k2a, b + 0.5 * h * k2b)
        k4a, k4b = rhs(a + h * k3a, b + h * k3b)
        a = a + (h / 6.0) * (k1a + 2.0 * k2a + 2.0 * k3a + k4a)
        b = b + (h / 6.0) * (k1b + 2.0 * k2b + 2.0 * k3b + k4b)

    if not np.all(np.isfinite(a.real)) or not np.all(np.isfinite(a.imag)):
        raise FloatingPointError("non-finite Heston transition A coefficients")
    if not np.all(np.isfinite(b.real)) or not np.all(np.isfinite(b.imag)):
        raise FloatingPointError("non-finite Heston transition B coefficients")
    return a, b


def heston_p_transition_cf(
    s_nodes: np.ndarray,
    x_prev: np.ndarray,
    *,
    dt: float,
    theta: HestonEstimationParams,
    rk_steps: int = 32,
) -> np.ndarray:
    """Conditional CF of next X=(return, variance) given current X.

    Returns an array with shape (n_nodes, n_observations).
    """

    x = np.asarray(x_prev, dtype=float)
    if x.ndim != 2 or x.shape[1] != 2:
        raise ValueError("x_prev must have shape (n_observations, 2)")
    if np.any(x[:, 1] < 0.0):
        raise ValueError("current variances must be non-negative")
    a, b = heston_p_transform_coefficients(
        s_nodes,
        dt=dt,
        theta=theta,
        rk_steps=rk_steps,
    )
    return np.exp(a[:, None] + b[:, None] * x[:, 1][None, :])
