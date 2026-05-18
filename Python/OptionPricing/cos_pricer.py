from __future__ import annotations

import math

import numpy as np

from OptionPricing.base import OptionPricer
from OptionPricing.types import CoefficientTensor, CosPricingConfig, HestonPricingParamsQ


class CosOptionPricer(OptionPricer):
    """Fang-Oosterlee COS pricer for European options using affine Heston CF coefficients."""

    _static_term_cache: dict[tuple[int, float], tuple[np.ndarray, np.ndarray]] = {}

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

    @classmethod
    def _get_static_terms(cls, n_cos: int, truncation_width: float) -> tuple[np.ndarray, np.ndarray]:
        key = (n_cos, truncation_width)
        cached = cls._static_term_cache.get(key)
        if cached is not None:
            return cached

        k_grid = np.arange(n_cos, dtype=float)
        a = -truncation_width
        b = truncation_width
        u = k_grid * np.pi / (b - a)
        v_call = cls._chi_psi_terms(k_grid, a, b)
        weights = np.ones_like(k_grid)
        weights[0] = 0.5
        common = weights * v_call * np.exp(-1j * u * a)
        cls._static_term_cache[key] = (u, common)
        return u, common

    @staticmethod
    def effective_truncation_width(variance, maturity_grid, pricing_config: CosPricingConfig) -> float:
        variance_arr = np.asarray(variance, dtype=float)
        maturities = np.asarray(maturity_grid, dtype=float)
        if variance_arr.size == 0 or maturities.size == 0:
            raise ValueError("variance and maturity_grid must be non-empty")
        if np.any(variance_arr < 0.0) or np.any(maturities < 0.0):
            raise ValueError("variance and maturity_grid must be non-negative")
        vol_time_scale = math.sqrt(max(float(np.max(variance_arr)), 1e-12) * max(float(np.max(maturities)), 1e-12))
        return max(0.5, pricing_config.truncation_width * vol_time_scale)

    def price_matrix(
        self,
        *,
        log_s,
        variance,
        strike_grid,
        maturity_grid,
        rate_grid,
        coefficients: CoefficientTensor,
        pricing_config: CosPricingConfig,
        dividend_yield_grid=None,
        option_type="call",
    ):
        pricing_config.validate()

        log_s = np.asarray(log_s, dtype=float).reshape(-1)
        variance = np.asarray(variance, dtype=float).reshape(-1)
        strikes = np.asarray(strike_grid, dtype=float)
        maturities = np.asarray(maturity_grid, dtype=float)
        rates = np.asarray(rate_grid, dtype=float)
        dividend_yields = np.zeros_like(rates) if dividend_yield_grid is None else np.asarray(dividend_yield_grid, dtype=float)
        option_types = np.asarray(option_type)

        if variance.shape != log_s.shape:
            raise ValueError("log_s and variance must have the same shape")
        if np.any(variance < 0.0):
            raise ValueError("variance must be non-negative")

        if maturities.ndim != 1 or rates.ndim != 1 or dividend_yields.ndim != 1:
            raise ValueError("maturity_grid, rate_grid, and dividend_yield_grid must be 1D")
        if rates.shape[0] != maturities.shape[0] or dividend_yields.shape[0] != maturities.shape[0]:
            raise ValueError("rate_grid and dividend_yield_grid must match maturity_grid length")

        if strikes.ndim == 1:
            if np.any(strikes <= 0.0):
                raise ValueError("strike_grid must contain strictly positive strikes")
            strike_tensor = np.broadcast_to(strikes[None, :, None], (log_s.shape[0], strikes.shape[0], maturities.shape[0]))
        elif strikes.ndim == 2:
            if strikes.shape != (log_s.shape[0], maturities.shape[0]):
                raise ValueError("2D strike_grid must have shape (n_observations, n_maturities)")
            if np.any(strikes <= 0.0):
                raise ValueError("strike_grid must contain strictly positive strikes")
            strike_tensor = strikes[:, None, :]
        elif strikes.ndim == 3:
            if strikes.shape[0] != log_s.shape[0] or strikes.shape[2] != maturities.shape[0]:
                raise ValueError("3D strike_grid must have shape (n_observations, n_strikes, n_maturities)")
            if np.any(strikes <= 0.0):
                raise ValueError("strike_grid must contain strictly positive strikes")
            strike_tensor = strikes
        else:
            raise ValueError("strike_grid must be 1D, 2D, or 3D")

        cf_a = np.asarray(coefficients.cf_a)
        cf_b = np.asarray(coefficients.cf_b)
        expected_shape = (pricing_config.n_cos, maturities.shape[0])
        truncation_width = self.effective_truncation_width(variance, maturities, pricing_config)
        expected_u, common = self._get_static_terms(pricing_config.n_cos, truncation_width)
        if cf_a.shape != expected_shape or cf_b.shape != expected_shape:
            raise ValueError("coefficient tensors must match (n_cos, n_maturities)")
        if not np.allclose(coefficients.u_grid, expected_u):
            raise ValueError("coefficient u_grid does not match pricing_config COS grid")
        u = expected_u

        n_obs = log_s.shape[0]
        n_strikes = strike_tensor.shape[1]
        prices_call = np.empty((n_obs, n_strikes, maturities.shape[0]), dtype=float)
        s0 = np.exp(log_s)

        for t_idx, (t_val, r_val, q_val) in enumerate(zip(maturities, rates, dividend_yields)):
            if t_val == 0.0:
                prices_call[:, :, t_idx] = np.maximum(s0[:, None] - strike_tensor[:, :, t_idx], 0.0)
                continue

            log_k = np.log(strike_tensor[:, :, t_idx])
            x = log_s[:, None] + (r_val - q_val) * t_val - log_k
            phase = np.exp(1j * u[:, None] * x.ravel()[None, :])
            affine = np.exp(cf_a[:, t_idx][:, None] + cf_b[:, t_idx][:, None] * variance[None, :])
            affine_by_strike = np.repeat(affine, n_strikes, axis=1)
            series_coeff = (common[:, None] * affine_by_strike)
            series = np.real(np.sum(series_coeff * phase, axis=0)).reshape(x.shape)
            prices_call[:, :, t_idx] = np.exp(-r_val * t_val) * strike_tensor[:, :, t_idx] * series

        prices_call = np.maximum(prices_call, 0.0)
        if option_types.ndim == 0:
            kind = str(option_types.item()).lower()
            if kind == "call":
                return prices_call
            if kind != "put":
                raise ValueError("option_type must be 'call' or 'put'")
            forward = s0[:, None, None] * np.exp((rates - dividend_yields)[None, None, :] * maturities[None, None, :])
            return np.maximum(prices_call - np.exp(-rates[None, None, :] * maturities[None, None, :]) * (forward - strike_tensor), 0.0)

        if option_types.shape != prices_call.shape:
            raise ValueError("option_type array must match output price shape")
        lower_types = np.char.lower(option_types.astype(str))
        if not np.all(np.isin(lower_types, ["call", "put"])):
            raise ValueError("option_type entries must be 'call' or 'put'")
        forward = s0[:, None, None] * np.exp((rates - dividend_yields)[None, None, :] * maturities[None, None, :])
        put_prices = np.maximum(prices_call - np.exp(-rates[None, None, :] * maturities[None, None, :]) * (forward - strike_tensor), 0.0)
        return np.where(lower_types == "call", prices_call, put_prices)

    def price_one(
        self,
        *,
        S: float,
        V: float,
        tau: float,
        K: float,
        option_type: str,
        model_params: HestonPricingParamsQ,
        coefficients: CoefficientTensor | None = None,
        pricing_config: CosPricingConfig | None = None,
    ) -> float:
        if S <= 0.0 or K <= 0.0:
            raise ValueError("S and K must be strictly positive")
        config = CosPricingConfig() if pricing_config is None else pricing_config
        coef = coefficients
        if coef is None:
            from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver

            width = self.effective_truncation_width(np.array([V], dtype=float), np.array([tau], dtype=float), config)
            u, _ = self._get_static_terms(config.n_cos, width)
            coef = HestonAnalyticCcfSolver().solve_coefficients(
                u,
                np.array([tau], dtype=float),
                model_params=model_params,
            )
        out = self.price_matrix(
            log_s=np.array([math.log(S)]),
            variance=np.array([V]),
            strike_grid=np.array([K]),
            maturity_grid=np.array([tau]),
            rate_grid=np.array([model_params.r]),
            dividend_yield_grid=np.array([model_params.q]),
            coefficients=coef,
            pricing_config=config,
            option_type=option_type,
        )
        return float(out[0, 0, 0])
