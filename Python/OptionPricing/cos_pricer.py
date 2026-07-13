from __future__ import annotations

import math

import numpy as np

from Models.Heston.parameters import HestonRiskNeutralParameters
from OptionPricing.base import OptionPricer
from OptionPricing.types import CoefficientTensor, PreparedFixedCosBasis, VarianceScaledCosConfig


class CosOptionPricer(OptionPricer):
    """Fang-Oosterlee COS pricing with an explicit effective-width convention."""

    _static_term_cache: dict[tuple[int, float], tuple[np.ndarray, np.ndarray]] = {}

    @staticmethod
    def _chi_psi_terms(k_grid: np.ndarray, a: float, b: float) -> np.ndarray:
        c = 0.0
        d = b
        omega = k_grid * np.pi / (b - a)
        exp_d = np.exp(d)
        exp_c = np.exp(c)
        chi = (
            np.cos(omega * (d - a)) * exp_d
            - np.cos(omega * (c - a)) * exp_c
            + omega * (np.sin(omega * (d - a)) * exp_d - np.sin(omega * (c - a)) * exp_c)
        ) / (1.0 + omega * omega)
        psi = np.empty_like(omega)
        psi[0] = d - c
        psi[1:] = (
            (np.sin(omega[1:] * (d - a)) - np.sin(omega[1:] * (c - a)))
            * (b - a)
            / (k_grid[1:] * np.pi)
        )
        return 2.0 / (b - a) * (chi - psi)

    @classmethod
    def _get_static_terms(cls, n_cos: int, effective_width: float) -> tuple[np.ndarray, np.ndarray]:
        if n_cos <= 1 or effective_width <= 0.0:
            raise ValueError("n_cos must be > 1 and effective_width must be positive")
        key = (int(n_cos), float(effective_width))
        cached = cls._static_term_cache.get(key)
        if cached is not None:
            return cached
        k_grid = np.arange(n_cos, dtype=float)
        a = -effective_width
        b = effective_width
        u = k_grid * np.pi / (b - a)
        payoff = cls._chi_psi_terms(k_grid, a, b)
        weights = np.ones_like(k_grid)
        weights[0] = 0.5
        terms = weights * payoff * np.exp(-1j * u * a)
        cls._static_term_cache[key] = (u, terms)
        return u, terms

    @staticmethod
    def variance_scaled_effective_width(
        variance,
        maturity_grid,
        config: VarianceScaledCosConfig,
    ) -> float:
        """Legacy/reference-only variance-scaled effective width."""

        config.validate()
        variance_arr = np.asarray(variance, dtype=float)
        maturities = np.asarray(maturity_grid, dtype=float)
        if variance_arr.size == 0 or maturities.size == 0:
            raise ValueError("variance and maturity_grid must be non-empty")
        if np.any(variance_arr < 0.0) or np.any(maturities < 0.0):
            raise ValueError("variance and maturity_grid must be non-negative")
        scale = math.sqrt(max(float(np.max(variance_arr)), 1e-12) * max(float(np.max(maturities)), 1e-12))
        return max(0.5, config.width_multiplier * scale)

    def prepare_fixed_basis(
        self,
        *,
        maturity: float,
        effective_width: float,
        n_cos: int,
        model_params: HestonRiskNeutralParameters,
    ) -> PreparedFixedCosBasis:
        """Prepare the grid, payoff expansion, and affine A/B once.

        The returned affine ``B`` remains available so a later implementation
        can evaluate price and initial-variance sensitivity together. That
        semi-analytical derivative is valid because this effective basis is
        fixed; it must not be used with a variance-dependent truncation rule.
        """

        from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver

        if maturity <= 0.0:
            raise ValueError("maturity must be strictly positive")
        u_grid, payoff_terms = self._get_static_terms(n_cos, effective_width)
        coefficients = HestonAnalyticCcfSolver().solve_coefficients(
            u_grid,
            np.array([maturity], dtype=float),
            model_params=model_params,
        )
        return PreparedFixedCosBasis(
            maturity=float(maturity),
            effective_width=float(effective_width),
            n_cos=int(n_cos),
            u_grid=u_grid,
            payoff_terms=payoff_terms,
            coefficients=coefficients,
        )

    def price_matrix_fixed_basis(
        self,
        *,
        log_s,
        variance,
        strike_grid,
        rate: float,
        dividend_yield: float,
        basis: PreparedFixedCosBasis,
        option_type="call",
    ):
        """Price one maturity using a prepared, variance-independent basis."""

        return self._price_matrix_with_effective_width(
            log_s=log_s,
            variance=variance,
            strike_grid=strike_grid,
            maturity_grid=np.array([basis.maturity], dtype=float),
            rate_grid=np.array([rate], dtype=float),
            dividend_yield_grid=np.array([dividend_yield], dtype=float),
            coefficients=basis.coefficients,
            n_cos=basis.n_cos,
            effective_width=basis.effective_width,
            option_type=option_type,
            prepared_payoff_terms=basis.payoff_terms,
        )

    def price_matrix_variance_scaled_reference(
        self,
        *,
        log_s,
        variance,
        strike_grid,
        maturity_grid,
        rate_grid,
        model_params: HestonRiskNeutralParameters,
        config: VarianceScaledCosConfig,
        dividend_yield_grid=None,
        option_type="call",
    ):
        """Explicit legacy/reference path whose grid scales with variance."""

        from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver

        width = self.variance_scaled_effective_width(variance, maturity_grid, config)
        u_grid, _ = self._get_static_terms(config.n_cos, width)
        coefficients = HestonAnalyticCcfSolver().solve_coefficients(
            u_grid,
            np.asarray(maturity_grid, dtype=float),
            model_params=model_params,
        )
        return self._price_matrix_with_effective_width(
            log_s=log_s,
            variance=variance,
            strike_grid=strike_grid,
            maturity_grid=maturity_grid,
            rate_grid=rate_grid,
            dividend_yield_grid=dividend_yield_grid,
            coefficients=coefficients,
            n_cos=config.n_cos,
            effective_width=width,
            option_type=option_type,
        )

    def price_matrix_with_explicit_effective_width(
        self,
        *,
        log_s,
        variance,
        strike_grid,
        maturity_grid,
        rate_grid,
        coefficients: CoefficientTensor,
        n_cos: int,
        effective_width: float,
        dividend_yield_grid=None,
        option_type="call",
    ):
        """Price a supplied coefficient tensor on one explicit fixed grid."""

        return self._price_matrix_with_effective_width(
            log_s=log_s,
            variance=variance,
            strike_grid=strike_grid,
            maturity_grid=maturity_grid,
            rate_grid=rate_grid,
            dividend_yield_grid=dividend_yield_grid,
            coefficients=coefficients,
            n_cos=n_cos,
            effective_width=effective_width,
            option_type=option_type,
        )

    def _price_matrix_with_effective_width(
        self,
        *,
        log_s,
        variance,
        strike_grid,
        maturity_grid,
        rate_grid,
        coefficients: CoefficientTensor,
        n_cos: int,
        effective_width: float,
        dividend_yield_grid=None,
        option_type="call",
        prepared_payoff_terms: np.ndarray | None = None,
    ):
        log_s = np.asarray(log_s, dtype=float).reshape(-1)
        variance = np.asarray(variance, dtype=float).reshape(-1)
        strikes = np.asarray(strike_grid, dtype=float)
        maturities = np.asarray(maturity_grid, dtype=float)
        rates = np.asarray(rate_grid, dtype=float)
        dividend_yields = np.zeros_like(rates) if dividend_yield_grid is None else np.asarray(dividend_yield_grid, dtype=float)
        option_types = np.asarray(option_type)
        if variance.shape != log_s.shape or np.any(variance < 0.0):
            raise ValueError("log_s and non-negative variance must have the same shape")
        if maturities.ndim != 1 or rates.shape != maturities.shape or dividend_yields.shape != maturities.shape:
            raise ValueError("maturity, rate, and dividend-yield grids must be aligned 1D arrays")

        if strikes.ndim == 1:
            if np.any(strikes <= 0.0):
                raise ValueError("strike_grid must contain strictly positive strikes")
            strike_tensor = np.broadcast_to(strikes[None, :, None], (log_s.size, strikes.size, maturities.size))
        elif strikes.ndim == 2:
            if strikes.shape != (log_s.size, maturities.size) or np.any(strikes <= 0.0):
                raise ValueError("2D strike_grid must have shape (n_observations, n_maturities) and be positive")
            strike_tensor = strikes[:, None, :]
        elif strikes.ndim == 3:
            if strikes.shape[0] != log_s.size or strikes.shape[2] != maturities.size or np.any(strikes <= 0.0):
                raise ValueError("3D strike_grid must align with observations/maturities and be positive")
            strike_tensor = strikes
        else:
            raise ValueError("strike_grid must be 1D, 2D, or 3D")

        expected_shape = (n_cos, maturities.size)
        u_grid, payoff_terms = self._get_static_terms(n_cos, effective_width)
        if prepared_payoff_terms is not None and not np.array_equal(prepared_payoff_terms, payoff_terms):
            raise ValueError("prepared fixed-basis payoff terms do not match the requested COS basis")
        if coefficients.cf_a.shape != expected_shape or coefficients.cf_b.shape != expected_shape:
            raise ValueError("coefficient tensors must match (n_cos, n_maturities)")
        if not np.allclose(coefficients.u_grid, u_grid, rtol=0.0, atol=1e-14):
            raise ValueError("coefficient u_grid does not match the explicit effective COS width")

        n_strikes = strike_tensor.shape[1]
        prices_call = np.empty((log_s.size, n_strikes, maturities.size), dtype=float)
        s0 = np.exp(log_s)
        for t_idx, (tau, rate, q) in enumerate(zip(maturities, rates, dividend_yields)):
            if tau == 0.0:
                prices_call[:, :, t_idx] = np.maximum(s0[:, None] - strike_tensor[:, :, t_idx], 0.0)
                continue
            x = log_s[:, None] + (rate - q) * tau - np.log(strike_tensor[:, :, t_idx])
            phase = np.exp(1j * u_grid[:, None] * x.ravel()[None, :])
            affine = np.exp(coefficients.cf_a[:, t_idx][:, None] + coefficients.cf_b[:, t_idx][:, None] * variance[None, :])
            series = np.real(np.sum(payoff_terms[:, None] * np.repeat(affine, n_strikes, axis=1) * phase, axis=0)).reshape(x.shape)
            prices_call[:, :, t_idx] = np.exp(-rate * tau) * strike_tensor[:, :, t_idx] * series

        prices_call = np.maximum(prices_call, 0.0)
        if option_types.ndim == 0:
            kind = str(option_types.item()).lower()
            if kind == "call":
                return prices_call
            if kind != "put":
                raise ValueError("option_type must be 'call' or 'put'")
            option_types = np.full(prices_call.shape, kind)
        if option_types.shape != prices_call.shape:
            raise ValueError("option_type array must match output price shape")
        lower_types = np.char.lower(option_types.astype(str))
        if not np.all(np.isin(lower_types, ["call", "put"])):
            raise ValueError("option_type entries must be 'call' or 'put'")
        forward = s0[:, None, None] * np.exp((rates - dividend_yields)[None, None, :] * maturities[None, None, :])
        puts = np.maximum(prices_call - np.exp(-rates[None, None, :] * maturities[None, None, :]) * (forward - strike_tensor), 0.0)
        return np.where(lower_types == "call", prices_call, puts)

    def price_one_variance_scaled_reference(
        self,
        *,
        S: float,
        V: float,
        tau: float,
        K: float,
        option_type: str,
        model_params: HestonRiskNeutralParameters,
        config: VarianceScaledCosConfig,
    ) -> float:
        if S <= 0.0 or K <= 0.0:
            raise ValueError("S and K must be strictly positive")
        output = self.price_matrix_variance_scaled_reference(
            log_s=np.array([math.log(S)]),
            variance=np.array([V]),
            strike_grid=np.array([K]),
            maturity_grid=np.array([tau]),
            rate_grid=np.array([model_params.r]),
            dividend_yield_grid=np.array([model_params.q]),
            model_params=model_params,
            config=config,
            option_type=option_type,
        )
        return float(output[0, 0, 0])
