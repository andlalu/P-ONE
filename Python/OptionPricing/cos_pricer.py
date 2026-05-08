from __future__ import annotations

import numpy as np

from OptionPricing.base import OptionPricer
from OptionPricing.types import CoefficientTensor, CosPricingConfig


class CosOptionPricer(OptionPricer):
    """Fang-Oosterlee COS pricer for European calls using provided CF coefficients."""

    @staticmethod
    def _chi_psi_terms(k_grid: np.ndarray, a: float, b: float) -> np.ndarray:
        c = 0.0
        d = b
        omega = k_grid * np.pi / (b - a)

        exp_d = np.exp(d)
        exp_c = np.exp(c)
        cos_d = np.cos(omega * (d - a))
        cos_c = np.cos(omega * (c - a))
        sin_d = np.sin(omega * (d - a))
        sin_c = np.sin(omega * (c - a))

        chi = (cos_d * exp_d - cos_c * exp_c) + omega * (sin_d * exp_d - sin_c * exp_c)
        chi = chi / (1.0 + omega * omega)

        psi = np.empty_like(omega)
        psi[0] = d - c
        psi[1:] = (sin_d[1:] - sin_c[1:]) * (b - a) / (k_grid[1:] * np.pi)

        return 2.0 / (b - a) * (chi - psi)

    def price_matrix(self, state_matrix, strike_grid, maturity_grid, rate_grid, coefficients: CoefficientTensor, pricing_config: CosPricingConfig):
        pricing_config.validate()

        log_s = np.asarray(state_matrix, dtype=float).reshape(-1)
        strikes = np.asarray(strike_grid, dtype=float)
        maturities = np.asarray(maturity_grid, dtype=float)
        rates = np.asarray(rate_grid, dtype=float)

        if strikes.ndim != 1 or np.any(strikes <= 0.0):
            raise ValueError("strike_grid must be a strictly positive 1D array")
        if maturities.ndim != 1 or rates.ndim != 1 or maturities.shape[0] != rates.shape[0]:
            raise ValueError("maturity_grid and rate_grid must be 1D with equal length")

        cf_values = np.asarray(coefficients.cf_values)
        if cf_values.shape != (pricing_config.n_cos, maturities.shape[0]):
            raise ValueError("coefficient tensor shape must match (n_cos, n_maturities)")

        k_grid = np.arange(pricing_config.n_cos, dtype=float)
        a = -pricing_config.truncation_width
        b = pricing_config.truncation_width
        u = k_grid * np.pi / (b - a)

        v_call = self._chi_psi_terms(k_grid, a, b)
        weights = np.ones_like(k_grid)
        weights[0] = 0.5

        prices = np.empty((log_s.shape[0], strikes.shape[0], maturities.shape[0]), dtype=float)

        for t_idx, (t_val, r_val) in enumerate(zip(maturities, rates)):
            if t_val == 0.0:
                s0 = np.exp(log_s)[:, None]
                prices[:, :, t_idx] = np.maximum(s0 - strikes[None, :], 0.0)
                continue

            base_cf = cf_values[:, t_idx]
            discount = np.exp(-r_val * t_val)

            x = log_s[:, None] - np.log(strikes)[None, :]
            phase = np.exp(1j * np.multiply.outer(u, x.ravel()))
            payoff_transform = (weights * v_call * np.exp(-1j * u * a) * base_cf)[:, None]
            series = np.real(np.sum(payoff_transform * phase, axis=0)).reshape(x.shape)
            prices[:, :, t_idx] = discount * strikes[None, :] * series

        return np.maximum(prices, 0.0)
