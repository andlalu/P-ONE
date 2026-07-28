from __future__ import annotations

import csv
import importlib.util
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
    """Read a panel sidecar, or return an empty dictionary."""

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


def clean_panel_file(run_root: str | Path, sample_id: int) -> Path:
    """Find one generated clean panel."""

    root = Path(run_root)
    stem = f"sample_{sample_id:03d}"
    for suffix in (".parquet", ".csv"):
        candidate = root / "panels_clean" / f"{stem}{suffix}"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"missing clean panel for sample {sample_id:03d} under {root / 'panels_clean'}"
    )


def noisy_panel_file(
    run_root: str | Path,
    scenario: str,
    sample_id: int,
    panel_format: str,
) -> Path:
    """Build the file name for one noisy panel."""

    if panel_format not in {"parquet", "csv"}:
        raise ValueError("panel_format must be 'parquet' or 'csv'")
    return (
        Path(run_root)
        / "panels_observed"
        / scenario
        / f"sample_{sample_id:03d}.{panel_format}"
    )


def write_persistent_factor_records(
    factors: list[dict[str, Any]],
    run_root: str | Path,
    scenario: str,
    sample_id: int,
) -> Path:
    """Write the persistent factors for one noisy panel."""

    target = Path(run_root) / "noise_factors" / scenario / f"sample_{sample_id:03d}.csv"
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.stem + ".tmp.csv")
    with temporary.open("w", newline="") as file_handle:
        writer = csv.DictWriter(
            file_handle,
            fieldnames=["week_index", "factor_0", "factor_1", "factor_2"],
        )
        writer.writeheader()
        writer.writerows(factors)
        file_handle.flush()
        os.fsync(file_handle.fileno())
    os.replace(temporary, target)
    return target


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

    rows = read_records(file_path)
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

    metadata = read_panel_metadata(file_path)
    metadata.update({"source": str(file_path), "iv_column": iv_name, "price_column": price_name})
    return OptionPanel(dates=tuple(dates), metadata=metadata).truncate_dates(max_dates)
