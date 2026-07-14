from __future__ import annotations

import csv
import importlib.util
import math
import os
from pathlib import Path
from typing import Iterable

import numpy as np

from DGPSimulation.types import HestonPath
from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_vega
from Models.Heston.parameters import HestonParameters, HestonPhysicalParameters, HestonRiskNeutralParameters
from OptionData.io import write_panel_metadata
from OptionPricing.cos_basis import FixedCosBasisConfig
from OptionPricing.cos_pricer import CosOptionPricer

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
    "model_vega",
    "iv_method",
]


def heston_q_from_p(params_p: HestonPhysicalParameters, eta_v: float) -> HestonRiskNeutralParameters:
    return HestonParameters(
        eta=params_p.eta,
        kappa=params_p.kappa,
        vbar=params_p.vbar,
        sigma_v=params_p.sigma_v,
        rho=params_p.rho,
        eta_v=eta_v,
        r=params_p.r,
        q=params_p.q,
    ).to_risk_neutral()


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
    params_p: HestonPhysicalParameters,
    eta_v: float,
    maturities_years: Iterable[float],
    log_moneyness: Iterable[float],
    atm_option_type: str,
    pricing_method: str,
    iv_method: str,
    pricer: CosOptionPricer,
    cos_basis: FixedCosBasisConfig,
) -> list[dict[str, object]]:
    if pricing_method.upper() != "COS":
        raise NotImplementedError(f"unsupported pricing_method: {pricing_method}")

    params_q = heston_q_from_p(params_p, eta_v)
    maturities = [float(x) for x in maturities_years]
    cos_basis.validate_requested_maturities(maturities)
    moneyness_grid = [float(x) for x in log_moneyness]
    rows: list[dict[str, object]] = []
    log_s_week = np.asarray(path.logS_week, dtype=float)
    v_week = np.asarray(path.V_week, dtype=float)
    spots = np.exp(log_s_week)
    n_obs = log_s_week.shape[0]
    n_moneyness = len(moneyness_grid)
    n_maturities = len(maturities)

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

    prices = np.empty_like(strike_tensor)
    for maturity_index, maturity in enumerate(maturities):
        basis = pricer.prepare_fixed_basis(
            maturity=maturity,
            effective_width=cos_basis.width_for_maturity(maturity),
            n_cos=cos_basis.generation_n_cos,
            model_params=params_q,
        )
        prices[:, :, maturity_index] = pricer.price_matrix_fixed_basis(
            log_s=log_s_week,
            variance=v_week,
            strike_grid=strike_tensor[:, :, maturity_index : maturity_index + 1],
            rate=params_p.r,
            dividend_yield=params_p.q,
            basis=basis,
            option_type=option_types[:, :, maturity_index : maturity_index + 1],
        )[:, :, 0]

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
                model_vega = black76_vega(
                    forward=forward,
                    strike=strike,
                    tau=tau,
                    vol=model_iv,
                    discount_factor=discount,
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
                        "model_vega": model_vega,
                        "iv_method": iv_method,
                    }
                )
    return rows


def parquet_available() -> bool:
    return importlib.util.find_spec("pandas") is not None and (
        importlib.util.find_spec("pyarrow") is not None or importlib.util.find_spec("fastparquet") is not None
    )


def write_panel(
    rows: list[dict[str, object]],
    target_without_suffix: Path,
    *,
    metadata: dict[str, object],
    panel_format: str,
) -> Path:
    target_without_suffix.parent.mkdir(parents=True, exist_ok=True)
    if panel_format == "parquet":
        if not parquet_available():
            raise RuntimeError("panel_format='parquet' requires pandas and pyarrow or fastparquet")
        import pandas as pd  # type: ignore[import-not-found]

        out = target_without_suffix.with_suffix(".parquet")
        temporary = out.with_name(out.stem + ".tmp.parquet")
        pd.DataFrame(rows, columns=PANEL_COLUMNS).to_parquet(temporary, index=False)
        if len(pd.read_parquet(temporary)) != len(rows):
            temporary.unlink(missing_ok=True)
            raise RuntimeError("atomic Parquet validation failed before publication")
        os.replace(temporary, out)
        write_panel_metadata(out, metadata)
        return out
    if panel_format != "csv":
        raise ValueError("panel_format must be 'parquet' or 'csv'")
    out = target_without_suffix.with_suffix(".csv")
    temporary = out.with_name(out.stem + ".tmp.csv")
    with temporary.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PANEL_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
        fh.flush()
        os.fsync(fh.fileno())
    with temporary.open(newline="") as fh:
        if sum(1 for _ in csv.DictReader(fh)) != len(rows):
            temporary.unlink(missing_ok=True)
            raise RuntimeError("atomic CSV validation failed before publication")
    os.replace(temporary, out)
    write_panel_metadata(out, metadata)
    return out


def read_panel_csv(file_path: str | Path) -> list[dict[str, str]]:
    with Path(file_path).open(newline="") as fh:
        return list(csv.DictReader(fh))
