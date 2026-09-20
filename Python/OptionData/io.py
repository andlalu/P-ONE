from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np

from OptionData.noise_common import NOISE_SCENARIOS
from OptionData.panel import OptionPanel, OptionPanelDate


def panel_metadata_path(file_path: str | Path) -> Path:
    path = Path(file_path)
    return path.with_suffix(path.suffix + ".metadata.json")


def parquet_available() -> bool:
    """Check whether Parquet support is installed."""

    return importlib.util.find_spec("pandas") is not None and (
        importlib.util.find_spec("pyarrow") is not None
        or importlib.util.find_spec("fastparquet") is not None
    )


def read_records(file_path: str | Path) -> list[dict[str, Any]]:
    """Read rows from CSV or Parquet."""

    path = Path(file_path)
    if path.suffix == ".parquet":
        import pandas as pd  # type: ignore[import-not-found]

        return pd.read_parquet(path).to_dict("records")
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def read_panel_metadata(file_path: str | Path) -> dict[str, Any]:
    """Read embedded Parquet metadata or a legacy panel sidecar."""

    path = Path(file_path)
    metadata: dict[str, Any] = {}
    if path.suffix == ".parquet" and path.exists():
        try:
            import pyarrow.parquet as pq  # type: ignore[import-not-found]

            raw_metadata = pq.read_schema(path).metadata or {}
            for raw_key, raw_value in raw_metadata.items():
                key = raw_key.decode("utf-8")
                if not key.startswith("p_one."):
                    continue
                name = key.removeprefix("p_one.")
                value = raw_value.decode("utf-8")
                if name in {"format_version", "sample_id"}:
                    metadata[name] = int(value)
                elif name in {"scenario_order", "cos_basis"}:
                    metadata[name] = json.loads(value)
                else:
                    metadata[name] = value
        except (ImportError, OSError, ValueError):
            pass
    sidecar = panel_metadata_path(path)
    if not sidecar.exists():
        return metadata
    with sidecar.open() as fh:
        value = json.load(fh)
    if not isinstance(value, dict):
        raise ValueError(f"panel metadata sidecar {sidecar} must contain a JSON object")
    metadata.update(value)
    return metadata


def write_panel_metadata(file_path: str | Path, metadata: dict[str, Any]) -> Path:
    sidecar = panel_metadata_path(file_path)
    temporary_sidecar = sidecar.with_name(sidecar.name + ".tmp")
    with temporary_sidecar.open("w") as fh:
        json.dump(metadata, fh, indent=2, sort_keys=True)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(temporary_sidecar, sidecar)
    return sidecar


def write_records(
    rows: list[dict[str, Any]],
    target_without_suffix: str | Path,
    *,
    panel_format: str,
    metadata: dict[str, Any] | None = None,
    columns: list[str] | None = None,
) -> Path:
    """Write rows and metadata without exposing a partial file."""

    target = Path(target_without_suffix)
    target.parent.mkdir(parents=True, exist_ok=True)
    if panel_format == "parquet":
        if not parquet_available():
            raise RuntimeError("panel_format='parquet' requires pandas and pyarrow or fastparquet")
        import pandas as pd  # type: ignore[import-not-found]

        output = target.with_suffix(".parquet")
        temporary = output.with_name(output.stem + ".tmp.parquet")
        pd.DataFrame(rows, columns=columns).to_parquet(temporary, index=False)
        if len(pd.read_parquet(temporary)) != len(rows):
            temporary.unlink(missing_ok=True)
            raise RuntimeError("atomic Parquet validation failed before publication")
        os.replace(temporary, output)
    elif panel_format == "csv":
        output = target.with_suffix(".csv")
        temporary = output.with_name(output.stem + ".tmp.csv")
        fieldnames = columns if columns is not None else (list(rows[0].keys()) if rows else [])
        with temporary.open("w", newline="") as file_handle:
            writer = csv.DictWriter(file_handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
            file_handle.flush()
            os.fsync(file_handle.fileno())
        with temporary.open(newline="") as file_handle:
            if sum(1 for _ in csv.DictReader(file_handle)) != len(rows):
                temporary.unlink(missing_ok=True)
                raise RuntimeError("atomic CSV validation failed before publication")
        os.replace(temporary, output)
    else:
        raise ValueError("panel_format must be 'parquet' or 'csv'")
    if metadata is not None:
        write_panel_metadata(output, metadata)
    return output


def write_panel(
    rows: list[dict[str, Any]],
    target_without_suffix: str | Path,
    *,
    metadata: dict[str, Any],
    panel_format: str,
) -> Path:
    """Write a clean panel in its standard column order."""

    from OptionData.clean_panel import PANEL_COLUMNS

    return write_records(
        rows,
        target_without_suffix,
        panel_format=panel_format,
        metadata=metadata,
        columns=PANEL_COLUMNS,
    )


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
    scenario: str | None = None,
    iv_column: str | None = None,
    price_column: str | None = None,
    max_dates: int | None = None,
) -> OptionPanel:
    """Load a generated CSV/Parquet panel into canonical date slices."""

    rows = read_records(file_path)
    if not rows:
        raise ValueError(f"{file_path} contains no rows")
    combined_scenarios = ("clean",) + NOISE_SCENARIOS
    has_scenario = "scenario" in rows[0]
    if scenario is not None:
        if scenario not in combined_scenarios:
            raise ValueError(f"unknown panel scenario: {scenario!r}")
        if not has_scenario:
            raise ValueError("scenario selection requires a combined panel with a scenario column")
        present = tuple(dict.fromkeys(str(row["scenario"]) for row in rows))
        if present != combined_scenarios:
            raise ValueError(
                "combined panel must contain exactly the deterministic scenario order "
                f"{combined_scenarios}, found {present}"
            )
        scenario_rank = {name: index for index, name in enumerate(combined_scenarios)}
        combined_order = [
            (
                scenario_rank[str(row["scenario"])],
                int(row["week_index"]),
                float(row["maturity_years"]),
                float(row["log_moneyness"]),
            )
            for row in rows
        ]
        if combined_order != sorted(combined_order):
            raise ValueError("combined panel rows are not in deterministic scenario and contract order")
        rows = [row for row in rows if str(row["scenario"]) == scenario]
        if not rows:
            raise ValueError(f"combined panel contains no rows for scenario {scenario!r}")
        sample_ids = {int(row["sample_id"]) for row in rows}
        if len(sample_ids) != 1:
            raise ValueError("one scenario panel must contain exactly one sample")
        order_keys = [
            (
                int(row["week_index"]),
                float(row["maturity_years"]),
                float(row["log_moneyness"]),
            )
            for row in rows
        ]
        if order_keys != sorted(order_keys):
            raise ValueError(f"scenario {scenario!r} rows are not in deterministic panel order")
        iv_column = "estimation_iv"
        price_column = "estimation_price"
    elif has_scenario and len({str(row["scenario"]) for row in rows}) > 1:
        raise ValueError("a scenario must be selected when loading a combined panel")
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

    metadata = read_panel_metadata(file_path)
    metadata.update(
        {
            "source": str(file_path),
            "scenario": scenario if scenario is not None else metadata.get("scenario"),
            "sample_id": (
                int(rows[0]["sample_id"])
                if scenario is not None
                else metadata.get("sample_id")
            ),
            "iv_column": iv_name,
            "price_column": price_name,
        }
    )
    return OptionPanel(dates=tuple(dates), metadata=metadata).truncate_dates(max_dates)
