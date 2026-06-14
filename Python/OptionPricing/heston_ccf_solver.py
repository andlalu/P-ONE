from __future__ import annotations

import numpy as np

from OptionPricing.base import CcfSolver
from OptionPricing.types import CoefficientTensor, HestonPricingParamsQ


def solve_constant_riccati(
    *,
    a: np.ndarray,
    b: np.ndarray,
    c: float,
    b0: np.ndarray,
    tau: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve B'=a+bB+cB^2 and int_0^tau B(s) ds for complex arrays.

    This Heston-oriented helper is shared by pricing and IS-CGMM transition
    code. The current pricing solver keeps its established closed form, while
    the P-transition CCF uses this general terminal-loading version.
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
            b0_col * (exp_bt - 1.0) / b_col + a_col * ((exp_bt - 1.0) / (b_col * b_col) - t_row / b_col),
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


class HestonAnalyticCcfSolver(CcfSolver):
    def solve_coefficients(self, u_grid, maturity_grid, model_params: HestonPricingParamsQ, measure_params=None):
        model_params.validate()
        u = np.asarray(u_grid, dtype=float)
        t = np.asarray(maturity_grid, dtype=float)
        if u.ndim != 1 or t.ndim != 1:
            raise ValueError("u_grid and maturity_grid must be 1D arrays")
        if np.any(t < 0.0):
            raise ValueError("maturities must be non-negative")

        iu = 1j * u[:, None]
        kappa = model_params.kappa + model_params.v_risk_premium
        sigma = model_params.sigma_v
        rho = model_params.rho
        vbar = model_params.vbar

        if sigma == 0.0:
            raise NotImplementedError("deterministic variance limit is not implemented in the Heston CCF solver")

        d = np.sqrt((rho * sigma * iu - kappa) ** 2 + sigma * sigma * (iu + u[:, None] ** 2))
        g = (kappa - rho * sigma * iu - d) / (kappa - rho * sigma * iu + d)
        exp_dt = np.exp(-d * t[None, :])

        one_minus_g_exp = 1.0 - g * exp_dt
        one_minus_g = 1.0 - g
        c = (kappa * vbar / (sigma * sigma)) * ((kappa - rho * sigma * iu - d) * t[None, :] - 2.0 * np.log(one_minus_g_exp / one_minus_g))
        d_term = ((kappa - rho * sigma * iu - d) / (sigma * sigma)) * ((1.0 - exp_dt) / one_minus_g_exp)
        if not np.all(np.isfinite(c.real)) or not np.all(np.isfinite(c.imag)):
            raise FloatingPointError("non-finite affine A coefficient values")
        if not np.all(np.isfinite(d_term.real)) or not np.all(np.isfinite(d_term.imag)):
            raise FloatingPointError("non-finite affine B coefficient values")

        return CoefficientTensor(u_grid=u, maturities=t, cf_a=c, cf_b=d_term)
