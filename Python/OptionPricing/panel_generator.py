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
        self._coeff_cache: dict[tuple[int, bytes], object] = {}
        self.panel_config.validate()

    def _cache_key(self, truncation_width: float) -> tuple[int, float, bytes]:
        return (
            self.pricing_config.n_cos,
            truncation_width,
            self.panel_config.maturities.astype(float).tobytes(),
        )

    def generate_panel(self, log_s_week, v_week):
        log_s = np.asarray(log_s_week, dtype=float)
        v_week = np.asarray(v_week, dtype=float)
        n_obs = log_s.shape[0]
        truncation_width = self.pricing_stack.option_pricer.effective_truncation_width(
            v_week,
            self.panel_config.maturities,
            self.pricing_config,
        )
        key = self._cache_key(truncation_width)
        coef = self._coeff_cache.get(key)
        if coef is None:
            u_grid, _ = self.pricing_stack.option_pricer._get_static_terms(
                self.pricing_config.n_cos,
                truncation_width,
            )
            coef = self.pricing_stack.ccf_solver.solve_coefficients(
                u_grid, self.panel_config.maturities, model_params=self.solver_params
            )
            self._coeff_cache[key] = coef

        prices = self.pricing_stack.option_pricer.price_matrix(
            log_s=log_s,
            variance=v_week,
            strike_grid=self.panel_config.strikes,
            maturity_grid=self.panel_config.maturities,
            rate_grid=self.panel_config.rates,
            dividend_yield_grid=np.full_like(self.panel_config.rates, self.solver_params.q),
            coefficients=coef,
            pricing_config=self.pricing_config,
            option_type="call",
        )
        return OptionPanel(
            prices=prices,
            observation_index=np.arange(n_obs),
            strikes=self.panel_config.strikes,
            maturities=self.panel_config.maturities,
            metadata={"engine": "heston_cos"},
        )
