from __future__ import annotations

import numpy as np

from OptionPricing.base import OptionPanelGenerator
from OptionPricing.types import (
    CosPricingConfig,
    HestonPricingParamsQ,
    OptionPanel,
    OptionPanelConfig,
    PricingStack,
)


class HestonOptionPanelGenerator(OptionPanelGenerator):
    def __init__(
        self,
        pricing_stack: PricingStack,
        panel_config: OptionPanelConfig,
        solver_params: HestonPricingParamsQ,
        pricing_config: CosPricingConfig,
    ):
        self.pricing_stack = pricing_stack
        self.panel_config = panel_config
        self.solver_params = solver_params
        self.pricing_config = pricing_config
        self.panel_config.validate()

    def generate_panel(self, log_s_week, v_week):
        log_s = np.asarray(log_s_week, dtype=float)
        _ = np.asarray(v_week, dtype=float)
        n_obs = log_s.shape[0]
        u_grid = np.linspace(0.0, 100.0, self.pricing_config.n_cos)
        coef = self.pricing_stack.ccf_solver.solve_coefficients(
            u_grid, self.panel_config.maturities, model_params=self.solver_params
        )
        prices = self.pricing_stack.option_pricer.price_matrix(
            state_matrix=log_s,
            strike_grid=self.panel_config.strikes,
            maturity_grid=self.panel_config.maturities,
            rate_grid=self.panel_config.rates,
            coefficients=coef,
            pricing_config=self.pricing_config,
        )
        return OptionPanel(
            prices=prices,
            observation_index=np.arange(n_obs),
            strikes=self.panel_config.strikes,
            maturities=self.panel_config.maturities,
            metadata={"engine": "heston_cos"},
        )
