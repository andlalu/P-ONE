from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from Models.Heston.affine_transform import solve_heston_pricing_riccati_branch_stable
from Models.Heston.parameters import HestonRiskNeutralParameters


@dataclass(frozen=True)
class CoefficientTensor:
    """Affine Heston coefficients on the requested grids."""

    u_grid: np.ndarray
    maturities: np.ndarray
    cf_a: np.ndarray
    cf_b: np.ndarray


class HestonCcf:
    r"""Heston's pricing CCF is affine in today's variance:

    ``phi(u, tau | V_t) = exp(A(u, tau) + B(u, tau) V_t)``.
    Here ``A`` is independent of variance and ``B`` carries its effect. The
    branch-stable form keeps wide transform grids on a consistent log branch.
    """

    def coefficients(
        self,
        u_grid: np.ndarray,
        maturity_grid: np.ndarray,
        model_params: HestonRiskNeutralParameters,
    ) -> CoefficientTensor:
        model_params.validate()
        u = np.asarray(u_grid, dtype=float)
        maturities = np.asarray(maturity_grid, dtype=float)
        if u.ndim != 1 or maturities.ndim != 1:
            raise ValueError("u_grid and maturity_grid must be 1D arrays")
        if np.any(maturities < 0.0):
            raise ValueError("maturities must be non-negative")

        sigma_v = model_params.sigma_v
        if sigma_v == 0.0:
            # This form divides by sigma_v**2; sigma_v=0 needs its own limit.
            raise NotImplementedError("deterministic variance limit is not implemented in the Heston CCF")
        affine_a, affine_b = solve_heston_pricing_riccati_branch_stable(
            u=u,
            tau=maturities,
            kappa=model_params.kappa,
            vbar=model_params.vbar,
            sigma_v=sigma_v,
            rho=model_params.rho,
        )
        if not np.all(np.isfinite(affine_a.real)) or not np.all(np.isfinite(affine_a.imag)):
            raise FloatingPointError("non-finite affine A coefficient values")
        if not np.all(np.isfinite(affine_b.real)) or not np.all(np.isfinite(affine_b.imag)):
            raise FloatingPointError("non-finite affine B coefficient values")

        return CoefficientTensor(
            u_grid=u,
            maturities=maturities,
            cf_a=affine_a,
            cf_b=affine_b,
        )
