from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import numpy as np

from ImpliedVolatility.black_iv import implied_vol_black76
from OptionPricing.types import OptionPanel

from Estimation.ISCGMM.types import EstimationPanel, PanelDate


def _read_records(file_path: str | Path) -> list[dict[str, Any]]:
    path = Path(file_path)
    if path.suffix == ".parquet":
        import pandas as pd  # type: ignore[import-not-found]

        return pd.read_parquet(path).to_dict("records")
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _choose_column(rows: list[dict[str, Any]], preferred: str | None, fallback: tuple[str, ...], label: str) -> str:
    if not rows:
        raise ValueError("panel file is empty")
    keys = set(rows[0].keys())
    if preferred is not None:
        if preferred not in keys:
            raise ValueError(f"{label} column {preferred!r} is missing")
        return preferred
    for name in fallback:
        if name in keys:
            return name
    raise ValueError(f"panel rows do not contain any supported {label} column: {fallback}")


def _optional_float(row: dict[str, Any], name: str) -> float | None:
    if name not in row:
        return None
    value = row[name]
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    return float(value)


def load_option_panel_data(
    file_path: str | Path,
    *,
    iv_column: str | None = None,
    price_column: str | None = None,
    max_dates: int | None = None,
) -> EstimationPanel:
    """Load an existing generated clean or observed panel into estimator slices."""

    rows = _read_records(file_path)
    if not rows:
        raise ValueError(f"{file_path} contains no rows")

    required = {"week_index", "t", "S", "logS", "maturity_years", "strike", "option_type", "r", "q"}
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"panel rows are missing required columns: {sorted(missing)}")

    iv_name = _choose_column(rows, iv_column, ("observed_iv", "model_iv", "clean_iv"), "IV")
    price_name = _choose_column(rows, price_column, ("observed_price", "model_price"), "price") if price_column is not None else None
    if price_name is None and any(name in rows[0] for name in ("observed_price", "model_price")):
        price_name = "observed_price" if "observed_price" in rows[0] else "model_price"

    sorted_rows = sorted(
        rows,
        key=lambda row: (
            int(row["week_index"]),
            float(row["maturity_years"]),
            float(row.get("log_moneyness", 0.0)),
            float(row["strike"]),
        ),
    )

    grouped: dict[int, list[dict[str, Any]]] = {}
    for row in sorted_rows:
        grouped.setdefault(int(row["week_index"]), []).append(row)

    slices: list[PanelDate] = []
    for week_index in sorted(grouped):
        group = grouped[week_index]
        first = group[0]
        true_v = _optional_float(first, "V")
        observed_price = None
        if price_name is not None:
            observed_price = np.array([float(row[price_name]) for row in group], dtype=float)
        clean_iv = None
        if "model_iv" in first:
            clean_iv = np.array([float(row["model_iv"]) for row in group], dtype=float)
        slices.append(
            PanelDate(
                date_index=week_index,
                time=float(first["t"]),
                log_spot=float(first["logS"]),
                spot=float(first["S"]),
                strikes=np.array([float(row["strike"]) for row in group], dtype=float),
                maturities=np.array([float(row["maturity_years"]) for row in group], dtype=float),
                option_types=np.array([str(row["option_type"]).lower() for row in group], dtype=str),
                observed_iv=np.array([float(row[iv_name]) for row in group], dtype=float),
                rates=np.array([float(row["r"]) for row in group], dtype=float),
                dividend_yields=np.array([float(row["q"]) for row in group], dtype=float),
                true_variance=true_v,
                observed_price=observed_price,
                clean_iv=clean_iv,
            )
        )

    panel = EstimationPanel(
        dates=tuple(slices),
        metadata={
            "source": str(file_path),
            "iv_column": iv_name,
            "price_column": price_name,
        },
    )
    return panel.truncate_dates(max_dates)


def panel_from_option_panel(
    option_panel: OptionPanel,
    *,
    log_spot: np.ndarray,
    rates: np.ndarray,
    dividend_yields: np.ndarray | None = None,
    option_type: str = "call",
    true_variance: np.ndarray | None = None,
) -> EstimationPanel:
    """Adapter for the older dense NPZ OptionPanel schema.

    The NPZ schema stores prices but not IVs, so IVs are reconstructed using
    Black-76. This adapter is intentionally thin and assumes one option type.
    """

    log_s = np.asarray(log_spot, dtype=float)
    prices = np.asarray(option_panel.prices, dtype=float)
    maturities = np.asarray(option_panel.maturities, dtype=float)
    strikes = np.asarray(option_panel.strikes, dtype=float)
    rate_grid = np.asarray(rates, dtype=float)
    q_grid = np.zeros_like(rate_grid) if dividend_yields is None else np.asarray(dividend_yields, dtype=float)
    if prices.shape != (log_s.shape[0], strikes.shape[0], maturities.shape[0]):
        raise ValueError("option_panel prices shape is inconsistent with log_spot, strikes, and maturities")
    if rate_grid.shape != maturities.shape or q_grid.shape != maturities.shape:
        raise ValueError("rates and dividend_yields must match maturities")

    slices: list[PanelDate] = []
    for date_idx, log_value in enumerate(log_s):
        spot = float(np.exp(log_value))
        flat_strikes: list[float] = []
        flat_maturities: list[float] = []
        flat_types: list[str] = []
        flat_iv: list[float] = []
        flat_rates: list[float] = []
        flat_q: list[float] = []
        flat_prices: list[float] = []
        for strike_idx, strike in enumerate(strikes):
            for maturity_idx, tau in enumerate(maturities):
                rate = float(rate_grid[maturity_idx])
                q_val = float(q_grid[maturity_idx])
                forward = spot * math.exp((rate - q_val) * float(tau))
                discount = math.exp(-rate * float(tau))
                price = float(prices[date_idx, strike_idx, maturity_idx])
                flat_iv.append(
                    implied_vol_black76(
                        price=price,
                        forward=forward,
                        strike=float(strike),
                        tau=float(tau),
                        discount_factor=discount,
                        option_type=option_type,
                        on_bounds="clip",
                    )
                )
                flat_strikes.append(float(strike))
                flat_maturities.append(float(tau))
                flat_types.append(option_type)
                flat_rates.append(rate)
                flat_q.append(q_val)
                flat_prices.append(price)
        slices.append(
            PanelDate(
                date_index=int(option_panel.observation_index[date_idx]),
                time=float(date_idx),
                log_spot=float(log_value),
                spot=spot,
                strikes=np.array(flat_strikes, dtype=float),
                maturities=np.array(flat_maturities, dtype=float),
                option_types=np.array(flat_types, dtype=str),
                observed_iv=np.array(flat_iv, dtype=float),
                rates=np.array(flat_rates, dtype=float),
                dividend_yields=np.array(flat_q, dtype=float),
                true_variance=None if true_variance is None else float(true_variance[date_idx]),
                observed_price=np.array(flat_prices, dtype=float),
            )
        )
    return EstimationPanel(dates=tuple(slices), metadata=dict(option_panel.metadata))
