from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from OptionPricing.noisy_panel import NOISE_SCENARIOS, _price_bounds, read_table

REQUIRED_COLUMNS = {
    "noise_scenario",
    "raw_contaminated_iv",
    "raw_price_before_rounding",
    "price_after_rounding",
    "observed_price",
    "observed_iv",
    "was_price_capped",
    "cap_direction",
    "noise_draw",
    "noise_seed",
}


def _load_noise_config(run_root: Path) -> dict:
    with (run_root / "config" / "experiment_config.json").open() as file_handle:
        experiment = json.load(file_handle)
    noise = experiment.get("noise")
    if not isinstance(noise, dict):
        raise ValueError("experiment configuration does not contain noise settings")
    return noise


def _panel_files(panel_dir: Path) -> list[Path]:
    return sorted(panel_dir.glob("sample_*.parquet")) + sorted(panel_dir.glob("sample_*.csv"))


def _enabled_scenarios(config: dict) -> list[str]:
    scenarios = config.get("scenarios", {})
    return [scenario for scenario in NOISE_SCENARIOS if scenario in scenarios]


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).lower() in {"true", "1", "yes"}


def validate_noisy_panels(*, run_root: str | Path) -> None:
    root = Path(run_root)
    noise_config = _load_noise_config(root)
    sigma_min = float(noise_config["sigma_min"])
    price_epsilon = float(noise_config["price_epsilon"])
    clean_files = _panel_files(root / "panels_clean")
    if not clean_files:
        raise AssertionError(f"no clean panel files found under {root / 'panels_clean'}")

    for clean_file in clean_files:
        sample_name = clean_file.name
        clean_rows = read_table(clean_file)
        for scenario in _enabled_scenarios(noise_config):
            observed_file = root / "panels_observed" / scenario / sample_name
            if not observed_file.exists() and observed_file.suffix == ".parquet":
                observed_file = observed_file.with_suffix(".csv")
            if not observed_file.exists():
                raise AssertionError(f"missing observed panel {observed_file}")
            rows = read_table(observed_file)
            if len(rows) != len(clean_rows):
                raise AssertionError(f"{observed_file} row count differs from {clean_file}")
            missing = REQUIRED_COLUMNS - set(rows[0])
            if missing:
                raise AssertionError(f"{observed_file} missing columns: {sorted(missing)}")
            for row in rows:
                observed_price = float(row["observed_price"])
                observed_iv = float(row["observed_iv"])
                if not math.isfinite(observed_price):
                    raise AssertionError(f"{observed_file} has non-finite observed_price")
                if not math.isfinite(observed_iv) or observed_iv < sigma_min:
                    raise AssertionError(f"{observed_file} has invalid observed_iv")
                lower, upper = _price_bounds(
                    float(row["S"]),
                    float(row["strike"]),
                    float(row["maturity_years"]),
                    float(row["r"]),
                    float(row["q"]),
                    str(row["option_type"]),
                )
                if observed_price < lower + price_epsilon - 1e-11 or observed_price > upper - price_epsilon + 1e-11:
                    raise AssertionError(f"{observed_file} violates no-arbitrage bounds")
                cap_direction = str(row["cap_direction"])
                was_capped = _as_bool(row["was_price_capped"])
                if cap_direction not in {"none", "lower", "upper"}:
                    raise AssertionError(f"{observed_file} has invalid cap_direction")
                if was_capped != (cap_direction != "none"):
                    raise AssertionError(f"{observed_file} has inconsistent cap flag")
            if scenario == "persistent_factor":
                factor_file = root / "noise_factors" / scenario / sample_name.replace(".parquet", ".csv")
                factor_file = factor_file.with_suffix(".csv")
                if not factor_file.exists():
                    raise AssertionError(f"missing persistent factor file {factor_file}")
                with factor_file.open(newline="") as fh:
                    factor_rows = list(csv.DictReader(fh))
                weeks = {int(row["week_index"]) for row in clean_rows}
                if len(factor_rows) != len(weeks):
                    raise AssertionError(f"{factor_file} has wrong number of factor rows")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    args = parser.parse_args()
    validate_noisy_panels(run_root=args.run_root)
    print(f"validated noisy panels under {args.run_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
