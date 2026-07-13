from __future__ import annotations

import numpy as np


def solve_heston_pricing_riccati_branch_stable(
    *,
    u: np.ndarray,
    tau: np.ndarray,
    kappa: float,
    vbar: float,
    sigma_v: float,
    rho: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return risk-neutral Heston A/B using the established log branch.

    The generic terminal-loading solver below is appropriate for the physical
    transition transform. Pricing retains this branch-stable representation:
    substituting the generic complex-log primitive can introduce non-equivalent
    ``2*pi*i`` branch jumps in ``A`` on wide transform grids.
    """

    u_values = np.asarray(u, dtype=float).reshape(-1)
    maturities = np.asarray(tau, dtype=float).reshape(-1)
    iu = 1j * u_values[:, None]
    sigma_sq = sigma_v * sigma_v
    d = np.sqrt((rho * sigma_v * iu - kappa) ** 2 + sigma_sq * (iu + u_values[:, None] ** 2))
    g = (kappa - rho * sigma_v * iu - d) / (kappa - rho * sigma_v * iu + d)
    exp_dt = np.exp(-d * maturities[None, :])
    one_minus_g_exp = 1.0 - g * exp_dt
    one_minus_g = 1.0 - g
    a = (kappa * vbar / sigma_sq) * (
        (kappa - rho * sigma_v * iu - d) * maturities[None, :]
        - 2.0 * np.log(one_minus_g_exp / one_minus_g)
    )
    b = ((kappa - rho * sigma_v * iu - d) / sigma_sq) * ((1.0 - exp_dt) / one_minus_g_exp)
    return a, b


def solve_constant_riccati(
    *,
    a: np.ndarray,
    b: np.ndarray,
    c: float,
    b0: np.ndarray,
    tau: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve ``B'=a+bB+cB^2`` and its time integral.

    This terminal-loading form is shared by risk-neutral pricing and the
    physical return/variance transition transform. Pricing keeps its existing
    branch-stable closed form where appropriate.
    """

    a_arr, b_arr, b0_arr = np.broadcast_arrays(
        np.asarray(a, dtype=complex),
        np.asarray(b, dtype=complex),
        np.asarray(b0, dtype=complex),
    )
    a_arr = a_arr.reshape(-1)
    b_arr = b_arr.reshape(-1)
    b0_arr = b0_arr.reshape(-1)
    t = np.asarray(tau, dtype=float).reshape(-1)
    if np.any(t < 0.0):
        raise ValueError("tau must be non-negative")

    if abs(c) < 1e-15:
        b_col = b_arr[:, None]
        a_col = a_arr[:, None]
        b0_col = b0_arr[:, None]
        t_row = t[None, :]
        small_b = np.abs(b_col) < 1e-14
        exp_bt = np.exp(b_col * t_row)
        riccati_b = np.where(
            small_b,
            b0_col + a_col * t_row,
            b0_col * exp_bt + a_col * (exp_bt - 1.0) / b_col,
        )
        integral_b = np.where(
            small_b,
            b0_col * t_row + 0.5 * a_col * t_row * t_row,
            b0_col * (exp_bt - 1.0) / b_col
            + a_col * ((exp_bt - 1.0) / (b_col * b_col) - t_row / b_col),
        )
        return riccati_b, integral_b

    d = np.sqrt(b_arr * b_arr - 4.0 * c * a_arr)
    r1 = (-b_arr + d) / (2.0 * c)
    r2 = (-b_arr - d) / (2.0 * c)
    riccati_b = np.empty((a_arr.shape[0], t.shape[0]), dtype=complex)
    integral_b = np.empty_like(riccati_b)

    near_r1 = np.abs(b0_arr - r1) < 1e-13
    near_r2 = np.abs(b0_arr - r2) < 1e-13
    root_mask = near_r1 | near_r2
    if np.any(near_r1):
        riccati_b[near_r1] = r1[near_r1, None]
        integral_b[near_r1] = r1[near_r1, None] * t[None, :]
    if np.any(near_r2):
        riccati_b[near_r2] = r2[near_r2, None]
        integral_b[near_r2] = r2[near_r2, None] * t[None, :]

    regular = ~root_mask
    if np.any(regular):
        y0 = (b0_arr[regular] - r1[regular]) / (b0_arr[regular] - r2[regular])
        y_tau = y0[:, None] * np.exp(d[regular, None] * t[None, :])
        riccati_b[regular] = (r1[regular, None] - y_tau * r2[regular, None]) / (1.0 - y_tau)
        integral_b[regular] = r1[regular, None] * t[None, :]
        integral_b[regular] += ((r1[regular] - r2[regular]) / d[regular])[:, None] * (
            np.log(1.0 - y0)[:, None] - np.log(1.0 - y_tau)
        )
    return riccati_b, integral_b
