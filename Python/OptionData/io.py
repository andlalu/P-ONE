from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np

from OptionData.panel import OptionPanel, OptionPanelDate


def panel_metadata_path(file_path: str | Path) -> Path:
    path = Path(file_path)
    return path.with_suffix(path.suffix + ".metadata.json")


def _read_records(file_path: str | Path) -> list[dict[str, Any]]:
    path = Path(file_path)
    if path.suffix == ".parquet":
        import pandas as pd  # type: ignore[import-not-found]

        return pd.read_parquet(path).to_dict("records")
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _read_metadata(file_path: str | Path) -> dict[str, Any]:
    sidecar = panel_metadata_path(file_path)
    if not sidecar.exists():
        return {}
    with sidecar.open() as fh:
        value = json.load(fh)
    if not isinstance(value, dict):
        raise ValueError(f"panel metadata sidecar {sidecar} must contain a JSON object")
    return value


def write_panel_metadata(file_path: str | Path, metadata: dict[str, Any]) -> Path:
    sidecar = panel_metadata_path(file_path)
    temporary_sidecar = sidecar.with_name(sidecar.name + ".tmp")
    with temporary_sidecar.open("w") as fh:
        json.dump(metadata, fh, indent=2, sort_keys=True)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(temporary_sidecar, sidecar)
    return sidecar


def _choose_column(rows: list[dict[str, Any]], preferred: str | None, fallback: tuple[str, ...], label: str) -> str:
    keys = set(rows[0])
    if preferred is not None:
        if preferred not in keys:
            raise ValueError(f"{label} column {preferred!r} is missing")
        return preferred
    for name in fallback:
        if name in keys:
            return name
    raise ValueError(f"panel rows do not contain any supported {label} column: {fallback}")


def _optional_float(row: dict[str, Any], name: str) -> float | None:
    value = row.get(name)
    if value is None or value == "":
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    return float(value)


def load_option_panel(
    file_path: str | Path,
    *,
    iv_column: str | None = None,
    price_column: str | None = None,
    max_dates: int | None = None,
) -> OptionPanel:
    """Load a generated CSV/Parquet panel into canonical date slices."""

    rows = _read_records(file_path)
    if not rows:
        raise ValueError(f"{file_path} contains no rows")
    required = {"week_index", "t", "S", "logS", "maturity_years", "strike", "option_type", "r", "q"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"panel rows are missing required columns: {sorted(missing)}")

    iv_name = _choose_column(rows, iv_column, ("observed_iv", "model_iv", "clean_iv"), "IV")
    price_name = None
    if price_column is not None:
        price_name = _choose_column(rows, price_column, (), "price")
    else:
        for candidate in ("observed_price", "model_price"):
            if candidate in rows[0]:
                price_name = candidate
                break

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

    dates: list[OptionPanelDate] = []
    for week_index in sorted(grouped):
        group = grouped[week_index]
        first = group[0]
        dates.append(
            OptionPanelDate(
                date_index=week_index,
                time=float(first["t"]),
                spot=float(first["S"]),
                log_spot=float(first["logS"]),
                strikes=np.array([float(row["strike"]) for row in group]),
                maturities=np.array([float(row["maturity_years"]) for row in group]),
                option_types=np.array([str(row["option_type"]).lower() for row in group]),
                observed_iv=np.array([float(row[iv_name]) for row in group]),
                rates=np.array([float(row["r"]) for row in group]),
                dividend_yields=np.array([float(row["q"]) for row in group]),
                observed_price=None if price_name is None else np.array([float(row[price_name]) for row in group]),
                clean_iv=None if "model_iv" not in first else np.array([float(row["model_iv"]) for row in group]),
                true_variance=_optional_float(first, "V"),
            )
        )

    metadata = _read_metadata(file_path)
    metadata.update({"source": str(file_path), "iv_column": iv_name, "price_column": price_name})
    return OptionPanel(dates=tuple(dates), metadata=metadata).truncate_dates(max_dates)
