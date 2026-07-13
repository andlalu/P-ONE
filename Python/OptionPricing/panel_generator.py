from __future__ import annotations

import numpy as np

from Models.Heston.parameters import HestonRiskNeutralParameters
from OptionPricing.base import OptionPriceCubeGenerator
from OptionPricing.cos_basis import FixedCosBasisConfig
from OptionPricing.types import (
    OptionPriceCube,
    OptionPriceCubeConfig,
    PricingStack,
)


class HestonOptionPriceCubeGenerator(OptionPriceCubeGenerator):
    def __init__(
        self,
        pricing_stack: PricingStack,
        panel_config: OptionPriceCubeConfig,
        solver_params: HestonRiskNeutralParameters,
        cos_basis: FixedCosBasisConfig,
    ):
        self.pricing_stack = pricing_stack
        self.panel_config = panel_config
        self.solver_params = solver_params
        self.cos_basis = cos_basis
        self._basis_cache: dict[float, object] = {}
        self.panel_config.validate()
        self.cos_basis.validate_requested_maturities(self.panel_config.maturities)

    def generate_panel(self, log_s_week, v_week):
        log_s = np.asarray(log_s_week, dtype=float)
        v_week = np.asarray(v_week, dtype=float)
        n_obs = log_s.shape[0]
        prices = np.empty((n_obs, self.panel_config.strikes.size, self.panel_config.maturities.size))
        for maturity_index, (maturity, rate) in enumerate(zip(self.panel_config.maturities, self.panel_config.rates)):
            maturity_value = float(maturity)
            basis = self._basis_cache.get(maturity_value)
            if basis is None:
                basis = self.pricing_stack.option_pricer.prepare_fixed_basis(
                    maturity=maturity_value,
                    effective_width=self.cos_basis.width_for_maturity(maturity_value),
                    n_cos=self.cos_basis.n_cos,
                    model_params=self.solver_params,
                )
                self._basis_cache[maturity_value] = basis
            prices[:, :, maturity_index] = self.pricing_stack.option_pricer.price_matrix_fixed_basis(
                log_s=log_s,
                variance=v_week,
                strike_grid=self.panel_config.strikes,
                rate=float(rate),
                dividend_yield=self.solver_params.q,
                basis=basis,
                option_type="call",
            )[:, :, 0]
        return OptionPriceCube(
            prices=prices,
            observation_index=np.arange(n_obs),
            strikes=self.panel_config.strikes,
            maturities=self.panel_config.maturities,
            metadata={"engine": "heston_cos", "cos_basis": self.cos_basis.to_metadata()},
        )
