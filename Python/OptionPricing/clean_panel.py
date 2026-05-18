from __future__ import annotations

import csv
import importlib.util
import math
import warnings
from pathlib import Path
from typing import Iterable

import numpy as np

from DGPSimulation.types import HestonParamsP, HestonPath
from ImpliedVolatility.black_iv import implied_vol_black76
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver
from OptionPricing.types import CosPricingConfig, HestonPricingParamsQ

PANEL_COLUMNS = [
    "run_id",
    "sample_id",
    "week_index",
    "t",
    "S",
    "logS",
    "V",
    "r",
    "q",
    "maturity_years",
    "expiry_time",
    "forward",
    "log_moneyness",
    "strike",
    "option_type",
    "is_otm",
    "pricing_method",
    "model_price",
    "model_iv",
    "iv_method",
]


def heston_q_from_p(params_p: HestonParamsP, eta_v: float) -> HestonPricingParamsQ:
    kappa_q = params_p.kappa - eta_v
    if kappa_q <= 0.0:
        raise ValueError("kappa - eta_v must be strictly positive")
    return HestonPricingParamsQ(
        kappa=kappa_q,
        vbar=params_p.kappa * params_p.vbar / kappa_q,
        sigma_v=params_p.sigma_v,
        rho=params_p.rho,
        r=params_p.r,
        q=params_p.q,
    )


def option_type_for_log_moneyness(log_moneyness: float, atm_option_type: str = "call") -> str:
    if log_moneyness < 0.0:
        return "put"
    if log_moneyness > 0.0:
        return "call"
    if atm_option_type not in {"call", "put"}:
        raise ValueError("atm_option_type must be 'call' or 'put'")
    return atm_option_type


def generate_clean_option_panel_rows(
    *,
    run_id: str,
    sample_id: int,
    path: HestonPath,
    params_p: HestonParamsP,
    eta_v: float,
    maturities_years: Iterable[float],
    log_moneyness: Iterable[float],
    atm_option_type: str,
    pricing_method: str,
    iv_method: str,
    pricer: CosOptionPricer,
    pricing_config: CosPricingConfig,
) -> list[dict[str, object]]:
    if pricing_method.upper() != "COS":
        raise NotImplementedError(f"unsupported pricing_method: {pricing_method}")

    params_q = heston_q_from_p(params_p, eta_v)
    maturities = [float(x) for x in maturities_years]
    moneyness_grid = [float(x) for x in log_moneyness]
    rows: list[dict[str, object]] = []
    log_s_week = np.asarray(path.logS_week, dtype=float)
    v_week = np.asarray(path.V_week, dtype=float)
    spots = np.exp(log_s_week)
    n_obs = log_s_week.shape[0]
    n_moneyness = len(moneyness_grid)
    n_maturities = len(maturities)

    truncation_width = pricer.effective_truncation_width(v_week, np.asarray(maturities, dtype=float), pricing_config)
    u_grid, _ = pricer._get_static_terms(pricing_config.n_cos, truncation_width)
    coefficients = HestonAnalyticCcfSolver().solve_coefficients(
        u_grid,
        np.asarray(maturities, dtype=float),
        model_params=params_q,
    )

    strike_tensor = np.empty((n_obs, n_moneyness, n_maturities), dtype=float)
    option_types = np.empty((n_obs, n_moneyness, n_maturities), dtype=object)
    forwards = np.empty((n_obs, n_maturities), dtype=float)

    for maturity_index, tau in enumerate(maturities):
        if tau <= 0.0:
            raise ValueError("maturities_years must be strictly positive")
        forwards[:, maturity_index] = spots * math.exp((params_p.r - params_p.q) * tau)
        for moneyness_index, lm in enumerate(moneyness_grid):
            strike_tensor[:, moneyness_index, maturity_index] = forwards[:, maturity_index] * math.exp(lm)
            option_types[:, moneyness_index, maturity_index] = option_type_for_log_moneyness(lm, atm_option_type)

    prices = pricer.price_matrix(
        log_s=log_s_week,
        variance=v_week,
        strike_grid=strike_tensor,
        maturity_grid=np.asarray(maturities, dtype=float),
        rate_grid=np.full(n_maturities, params_p.r, dtype=float),
        dividend_yield_grid=np.full(n_maturities, params_p.q, dtype=float),
        coefficients=coefficients,
        pricing_config=pricing_config,
        option_type=option_types,
    )

    for week_index, (t, log_s, v) in enumerate(zip(path.t_week, path.logS_week, path.V_week)):
        S = float(math.exp(float(log_s)))
        for maturity_index, tau in enumerate(maturities):
            discount = math.exp(-params_p.r * tau)
            forward = float(forwards[week_index, maturity_index])
            expiry_time = float(t) + tau
            for moneyness_index, lm in enumerate(moneyness_grid):
                strike = float(strike_tensor[week_index, moneyness_index, maturity_index])
                option_type = str(option_types[week_index, moneyness_index, maturity_index])
                model_price = float(prices[week_index, moneyness_index, maturity_index])
                model_iv = implied_vol_black76(
                    price=model_price,
                    forward=forward,
                    strike=strike,
                    tau=tau,
                    discount_factor=discount,
                    option_type=option_type,
                    on_bounds="raise",
                )
                rows.append(
                    {
                        "run_id": run_id,
                        "sample_id": sample_id,
                        "week_index": week_index,
                        "t": float(t),
                        "S": S,
                        "logS": float(log_s),
                        "V": float(v),
                        "r": params_p.r,
                        "q": params_p.q,
                        "maturity_years": tau,
                        "expiry_time": expiry_time,
                        "forward": forward,
                        "log_moneyness": lm,
                        "strike": strike,
                        "option_type": option_type,
                        "is_otm": True,
                        "pricing_method": pricing_method.upper(),
                        "model_price": model_price,
                        "model_iv": model_iv,
                        "iv_method": iv_method,
                    }
                )
    return rows


def parquet_available() -> bool:
    return importlib.util.find_spec("pandas") is not None and (
        importlib.util.find_spec("pyarrow") is not None or importlib.util.find_spec("fastparquet") is not None
    )


def write_panel(rows: list[dict[str, object]], target_without_suffix: Path) -> Path:
    target_without_suffix.parent.mkdir(parents=True, exist_ok=True)
    if parquet_available():
        import pandas as pd  # type: ignore[import-not-found]

        out = target_without_suffix.with_suffix(".parquet")
        pd.DataFrame(rows, columns=PANEL_COLUMNS).to_parquet(out, index=False)
        return out

    warnings.warn("Parquet dependencies unavailable; writing clean option panel as CSV.", RuntimeWarning)
    out = target_without_suffix.with_suffix(".csv")
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PANEL_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return out


def read_panel_csv(file_path: str | Path) -> list[dict[str, str]]:
    with Path(file_path).open(newline="") as fh:
        return list(csv.DictReader(fh))
