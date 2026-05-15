from __future__ import annotations

import numpy as np

from OptionPricing.base import CcfSolver
from OptionPricing.types import CoefficientTensor, HestonPricingParamsQ


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
            a = np.exp(-0.5 * u[:, None] ** 2 * vbar * t[None, :])
            return CoefficientTensor(u_grid=u, maturities=t, cf_values=a.astype(complex))

        d = np.sqrt((rho * sigma * iu - kappa) ** 2 + sigma * sigma * (iu + u[:, None] ** 2))
        g = (kappa - rho * sigma * iu - d) / (kappa - rho * sigma * iu + d)
        exp_dt = np.exp(-d * t[None, :])

        one_minus_g_exp = 1.0 - g * exp_dt
        one_minus_g = 1.0 - g
        c = (kappa * vbar / (sigma * sigma)) * ((kappa - rho * sigma * iu - d) * t[None, :] - 2.0 * np.log(one_minus_g_exp / one_minus_g))
        d_term = ((kappa - rho * sigma * iu - d) / (sigma * sigma)) * ((1.0 - exp_dt) / one_minus_g_exp)
        cf_values = np.exp(c + d_term * vbar)

        if not np.all(np.isfinite(cf_values.real)) or not np.all(np.isfinite(cf_values.imag)):
            raise FloatingPointError("non-finite coefficient values")

        return CoefficientTensor(u_grid=u, maturities=t, cf_values=cf_values)
