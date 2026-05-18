from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from DGPSimulation.heston_simulator import HestonPathSimulator
from DGPSimulation.types import HestonParamsP, HestonSimConfig
from DGPSimulation.variance_drawers import AndersenQeVarianceDrawer
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver
from OptionPricing.types import CosPricingConfig, HestonPricingParamsQ


@dataclass(frozen=True)
class MonteCarloEstimate:
    mean: float
    standard_error: float
    n_paths: int


def cos_heston_price(
    *,
    spot: float,
    variance: float,
    tau: float,
    strike: float,
    option_type: str,
    params: HestonPricingParamsQ,
    config: CosPricingConfig,
) -> float:
    pricer = CosOptionPricer()
    width = pricer.effective_truncation_width(np.array([variance]), np.array([tau]), config)
    u_grid, _ = pricer._get_static_terms(config.n_cos, width)
    coefficients = HestonAnalyticCcfSolver().solve_coefficients(
        u_grid,
        np.array([tau], dtype=float),
        model_params=params,
    )
    prices = pricer.price_matrix(
        log_s=np.array([math.log(spot)]),
        variance=np.array([variance]),
        strike_grid=np.array([strike]),
        maturity_grid=np.array([tau]),
        rate_grid=np.array([params.r]),
        dividend_yield_grid=np.array([params.q]),
        coefficients=coefficients,
        pricing_config=config,
        option_type=option_type,
    )
    return float(prices[0, 0, 0])


def mc_heston_price(
    *,
    spot: float,
    variance: float,
    n_weeks: int,
    strike: float,
    option_type: str,
    params: HestonPricingParamsQ,
    n_paths: int,
    base_seed: int = 37_000,
) -> MonteCarloEstimate:
    if option_type not in {"call", "put"}:
        raise ValueError("option_type must be 'call' or 'put'")

    dgp_params = HestonParamsP(
        eta=0.0,
        kappa=params.kappa,
        vbar=params.vbar,
        sigma_v=params.sigma_v,
        rho=params.rho,
        r=params.r,
        q=params.q,
    )
    tau = n_weeks * 5.0 / 252.0
    discounted_payoffs = np.empty(n_paths, dtype=float)

    for path_idx in range(n_paths):
        sim_config = HestonSimConfig(
            delta=1.0 / 252.0,
            m_week=5,
            t_week=n_weeks,
            burnin_days=0,
            s0=spot,
            v0=variance,
            seed=base_seed + path_idx,
            return_daily=False,
        )
        path = HestonPathSimulator(
            params=dgp_params,
            config=sim_config,
            variance_drawer=AndersenQeVarianceDrawer(),
        ).simulate()
        terminal_spot = math.exp(float(path.logS_week[-1]))
        if option_type == "call":
            payoff = max(terminal_spot - strike, 0.0)
        else:
            payoff = max(strike - terminal_spot, 0.0)
        discounted_payoffs[path_idx] = math.exp(-params.r * tau) * payoff

    return MonteCarloEstimate(
        mean=float(np.mean(discounted_payoffs)),
        standard_error=float(np.std(discounted_payoffs, ddof=1) / math.sqrt(n_paths)),
        n_paths=n_paths,
    )


def forward_strike(*, spot: float, tau: float, log_moneyness: float, r: float, q: float) -> float:
    forward = spot * math.exp((r - q) * tau)
    return forward * math.exp(log_moneyness)
