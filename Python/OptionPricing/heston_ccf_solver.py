from __future__ import annotations

import numpy as np

from Models.Heston.affine_transform import solve_heston_pricing_riccati_branch_stable
from Models.Heston.parameters import HestonRiskNeutralParameters
from OptionPricing.base import CcfSolver
from OptionPricing.types import CoefficientTensor


class HestonAnalyticCcfSolver(CcfSolver):
    """Thin risk-neutral Heston wrapper preserving the established log branch."""

    def solve_coefficients(self, u_grid, maturity_grid, model_params: HestonRiskNeutralParameters, measure_params=None):
        model_params.validate()
        u = np.asarray(u_grid, dtype=float)
        t = np.asarray(maturity_grid, dtype=float)
        if u.ndim != 1 or t.ndim != 1:
            raise ValueError("u_grid and maturity_grid must be 1D arrays")
        if np.any(t < 0.0):
            raise ValueError("maturities must be non-negative")

        kappa = model_params.kappa
        sigma = model_params.sigma_v
        if sigma == 0.0:
            raise NotImplementedError("deterministic variance limit is not implemented in the Heston CCF solver")
        c, d_term = solve_heston_pricing_riccati_branch_stable(
            u=u,
            tau=t,
            kappa=kappa,
            vbar=model_params.vbar,
            sigma_v=sigma,
            rho=model_params.rho,
        )
        if not np.all(np.isfinite(c.real)) or not np.all(np.isfinite(c.imag)):
            raise FloatingPointError("non-finite affine A coefficient values")
        if not np.all(np.isfinite(d_term.real)) or not np.all(np.isfinite(d_term.imag)):
            raise FloatingPointError("non-finite affine B coefficient values")

        return CoefficientTensor(u_grid=u, maturities=t, cf_a=c, cf_b=d_term)
