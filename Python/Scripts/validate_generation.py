from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

from DGPSimulation.io import load_heston_path_npz
from OptionPricing.clean_panel import option_type_for_log_moneyness


def _panel_files(root: Path) -> list[Path]:
    panel_dir = root / "panels_clean"
    return sorted(panel_dir.glob("sample_*.parquet")) + sorted(panel_dir.glob("sample_*.csv"))


def _read_csv_panel(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _read_panel(path: Path) -> list[dict[str, object]]:
    if path.suffix == ".csv":
        return _read_csv_panel(path)
    import pandas as pd  # type: ignore[import-not-found]

    return pd.read_parquet(path).to_dict("records")


def validate_run(*, run_root: str | Path, expected_samples: int) -> None:
    root = Path(run_root)
    path_files = sorted((root / "paths").glob("sample_*.npz"))
    panel_files = _panel_files(root)
    if len(path_files) != expected_samples:
        raise AssertionError(f"expected {expected_samples} path files, found {len(path_files)}")
    if len(panel_files) != expected_samples:
        raise AssertionError(f"expected {expected_samples} panel files, found {len(panel_files)}")

    for path_file, panel_file in zip(path_files, panel_files):
        path, params, config = load_heston_path_npz(path_file)
        if len(path.t_week) != config.t_week + 1:
            raise AssertionError(f"{path_file} has wrong weekly length")
        if len(path.dlogS_week) != config.t_week:
            raise AssertionError(f"{path_file} has wrong weekly return length")
        if path.logS_daily is None or path.V_daily is None:
            raise AssertionError(f"{path_file} missing daily arrays")
        for name, values in {
            "t_week": path.t_week,
            "logS_week": path.logS_week,
            "V_week": path.V_week,
            "dlogS_week": path.dlogS_week,
            "logS_daily": path.logS_daily,
            "V_daily": path.V_daily,
        }.items():
            if not np.all(np.isfinite(values)):
                raise AssertionError(f"{path_file} has non-finite {name}")
        if np.any(path.V_week < 0.0) or np.any(path.V_daily < 0.0):
            raise AssertionError(f"{path_file} has negative variance")

        rows = _read_panel(panel_file)
        if not rows:
            raise AssertionError(f"{panel_file} is empty")
        maturities = {float(row["maturity_years"]) for row in rows}
        moneyness = {float(row["log_moneyness"]) for row in rows}
        expected_rows = (config.t_week + 1) * len(maturities) * len(moneyness)
        if len(rows) != expected_rows:
            raise AssertionError(f"{panel_file} expected {expected_rows} rows, found {len(rows)}")
        for row in rows:
            price = float(row["model_price"])
            iv = float(row["model_iv"])
            if not math.isfinite(price) or price < 0.0:
                raise AssertionError(f"{panel_file} has invalid price")
            if not math.isfinite(iv) or iv < 0.0:
                raise AssertionError(f"{panel_file} has invalid IV")
            option_type = str(row["option_type"])
            lm = float(row["log_moneyness"])
            atm_type = option_type if lm == 0.0 else "call"
            if lm != 0.0 and option_type != option_type_for_log_moneyness(lm, atm_type):
                raise AssertionError(f"{panel_file} violates OTM option type rule")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    args = parser.parse_args()
    validate_run(run_root=args.run_root, expected_samples=args.expected_samples)
    print(f"validated {args.expected_samples} samples under {args.run_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
