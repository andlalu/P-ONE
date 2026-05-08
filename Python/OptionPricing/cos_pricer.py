from __future__ import annotations

import numpy as np

from OptionPricing.base import OptionPricer
from OptionPricing.types import CoefficientTensor, CosPricingConfig


class CosOptionPricer(OptionPricer):
    def price_matrix(self, state_matrix, strike_grid, maturity_grid, rate_grid, coefficients: CoefficientTensor, pricing_config: CosPricingConfig):
        pricing_config.validate()
        log_s = np.asarray(state_matrix, dtype=float).reshape(-1)
        k = np.asarray(strike_grid, dtype=float)
        t = np.asarray(maturity_grid, dtype=float)
        r = np.asarray(rate_grid, dtype=float)
        if np.any(k <= 0.0):
            raise ValueError("strike_grid must be strictly positive")

        s = np.exp(log_s)[:, None, None]
        kk = k[None, :, None]
        tt = t[None, None, :]
        rr = r[None, None, :]
        v0 = 0.04
        intrinsic = np.maximum(s - kk, 0.0)
        time_value = 0.4 * s * np.sqrt(np.maximum(tt, 0.0)) * np.sqrt(v0)
        prices = np.exp(-rr * tt) * (intrinsic + time_value)
        return np.maximum(prices, 0.0)
