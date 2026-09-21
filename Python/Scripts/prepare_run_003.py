"""Extend the frozen run_001 samples with design D for local run_003 generation.

Only samples 0-99 are handled here. The four original panel slices are copied
from the checksum-verified archive, never regenerated. Their run_id is the sole
changed original field; design D is newly generated with run_003 settings.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

import pandas as pd
from DGPSimulation.io import load_heston_path_npz
from OptionData.add_noise import generate_noisy_panel_rows
from OptionData.noise_common import scenario_seed
from Scripts.experiment_config import load_experiment_config
from Scripts.sample_run import (
    _combined_noisy_row,
    _git_state,
    _initial_record,
    _record_artifacts,
    _write_combined_panel,
    atomic_json,
    initialise_or_verify_run,
    sample_directory,
    sha256_file,
    utc_now,
    validate_sample_artifacts,
)

SOURCE_CONFIG_HASH = "b74031ebf443c40496c2b9b3b4088df11ce49d761cc154d78050268dd015873f"
SOURCE_GIT_SHA = "a9e17ee3880eddcd5cee31dfaf4854bdcf15d2e1"
ORIGINAL_SCENARIOS = ("clean", "low_iid", "spatial_corr", "persistent_factor")


def _read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _verify_source(source_root: Path, sample_id: int) -> tuple[Path, Path, dict]:
    source_run = _read_json(source_root / "run.json")
    if (
        source_run.get("run_id") != "run_001"
        or source_run.get("configuration_hash") != SOURCE_CONFIG_HASH
        or source_run.get("git_sha") != SOURCE_GIT_SHA
    ):
        raise ValueError("source is not the frozen 100-sample tick-rounded run_001")
    source_dir = sample_directory(source_root, sample_id)
    source_record = _read_json(source_dir / "record.json")
    if (
        source_record.get("run_id") != "run_001"
        or source_record.get("sample_id") != sample_id
        or source_record.get("configuration_hash") != SOURCE_CONFIG_HASH
    ):
        raise ValueError(f"sample {sample_id}: source record provenance mismatch")
    for name, filename in (("path", "path.npz"), ("panels", "panels.parquet")):
        file_path = source_dir / filename
        artifact = source_record["artifacts"][name]
        if sha256_file(file_path) != artifact["sha256"] or file_path.stat().st_size != artifact["bytes"]:
            raise ValueError(f"sample {sample_id}: source {name} checksum mismatch")
    return source_dir / "path.npz", source_dir / "panels.parquet", source_record


def _assert_original_slices(source: pd.DataFrame, output: pd.DataFrame) -> None:
    original = output[output["scenario"].isin(ORIGINAL_SCENARIOS)].reset_index(drop=True)
    if len(original) != len(source):
        raise AssertionError("original four-scenario row count changed")
    expected = source.copy()
    expected["run_id"] = "run_003"
    pd.testing.assert_frame_equal(
        original[list(source.columns)], expected, check_dtype=False, check_exact=True
    )


def prepare_sample(config, source_root: Path, sample_id: int) -> str:
    if config.run_id != "run_003" or config.noise is None or config.noise.tick_size != 0.01:
        raise ValueError("this import requires run_003 with the original 0.01 tick")
    if not 0 <= sample_id < 100:
        raise ValueError("only frozen run_001 samples 0-99 may be imported")
    source_path, source_panel, source_record = _verify_source(source_root, sample_id)
    target_dir = sample_directory(config.output_root, sample_id)
    target_dir.mkdir(parents=True, exist_ok=True)
    target_path = target_dir / "path.npz"
    target_panel = target_dir / "panels.parquet"
    record_path = target_dir / "record.json"
    if record_path.exists():
        old = _read_json(record_path)
        if old.get("status") == "generated" and old.get("configuration_hash") == config.experiment_config_hash:
            if old.get("source_run_001", {}).get("panels_sha256") != source_record["artifacts"]["panels"]["sha256"]:
                raise ValueError(f"sample {sample_id}: source provenance changed")
            if (
                sha256_file(target_path) == old["artifacts"]["path"]["sha256"]
                and sha256_file(target_panel) == old["artifacts"]["panels"]["sha256"]
            ):
                return "verified-existing"

    source_frame = pd.read_parquet(source_panel)
    if tuple(source_frame["scenario"].drop_duplicates()) != ORIGINAL_SCENARIOS:
        raise ValueError(f"sample {sample_id}: source scenario order mismatch")
    if len(source_frame) != 4 * 7890 or source_frame["sample_id"].nunique() != 1:
        raise ValueError(f"sample {sample_id}: source panel shape mismatch")
    source_frame = source_frame.copy()
    source_frame["run_id"] = config.run_id
    clean_rows = source_frame[source_frame["scenario"] == "clean"].to_dict("records")
    _, parameters, _ = load_heston_path_npz(source_path)
    seed = scenario_seed(config.noise.base_seed, sample_id, "variance_linked_factor")
    noisy_rows, factors = generate_noisy_panel_rows(
        clean_rows,
        scenario="variance_linked_factor",
        seed=seed,
        config=config.noise,
        params_p=parameters,
    )
    factor_by_week = {int(item["week_index"]): item for item in factors}
    new_rows = source_frame.to_dict("records")
    for row in new_rows:
        row["raw_noisy_price"] = row["raw_price_before_rounding"]
        row["variance_linked_factor"] = None
    new_rows.extend(
        _combined_noisy_row(row, "variance_linked_factor", factor_by_week[int(row["week_index"])])
        for row in noisy_rows
    )
    temporary_path = target_dir / "path.tmp.npz"
    shutil.copy2(source_path, temporary_path)
    if sha256_file(temporary_path) != source_record["artifacts"]["path"]["sha256"]:
        raise ValueError(f"sample {sample_id}: copied path checksum mismatch")
    os.replace(temporary_path, target_path)
    git_sha, _ = _git_state()
    _write_combined_panel(new_rows, target_panel, config=config, sample_id=sample_id, git_sha=git_sha)
    _assert_original_slices(source_frame, pd.read_parquet(target_panel))
    validation = validate_sample_artifacts(target_path, target_panel, config, sample_id)
    if not validation["passed"]:
        raise AssertionError(f"sample {sample_id}: run_003 validation failed")
    record = _initial_record(
        config, sample_id, generation_only=True, state_only=False, scenario=None, max_dates=None
    )
    record["validation"] = validation
    _record_artifacts(record, target_path, target_panel, validation)
    record["source_run_001"] = {
        "configuration_hash": SOURCE_CONFIG_HASH,
        "git_sha": SOURCE_GIT_SHA,
        "path_sha256": source_record["artifacts"]["path"]["sha256"],
        "panels_sha256": source_record["artifacts"]["panels"]["sha256"],
        "original_scenarios": list(ORIGINAL_SCENARIOS),
        "original_rows_preserved": 4 * 7890,
    }
    record["status"] = "generated"
    record["current_stage"] = "validated"
    record["completed_at_utc"] = utc_now()
    atomic_json(record_path, record)
    return "imported"


def _prepare_one(config_path: str, source_root: str, sample_id: int) -> tuple[int, str]:
    config = load_experiment_config(config_path)
    return sample_id, prepare_sample(config, Path(source_root), sample_id)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--sample-start", type=int, default=0)
    parser.add_argument("--sample-end", type=int, default=100)
    parser.add_argument("--sample-workers", type=int, default=1)
    args = parser.parse_args()
    if not 0 <= args.sample_start < args.sample_end <= 100:
        parser.error("source sample range must be within 0-100 (exclusive end)")
    if args.sample_workers < 1:
        parser.error("--sample-workers must be positive")
    config = load_experiment_config(args.config)
    initialise_or_verify_run(config.output_root, config)
    source_root = Path(args.source_root).expanduser().resolve()
    if args.sample_workers == 1:
        for sample_id in range(args.sample_start, args.sample_end):
            print(f"sample_{sample_id:03d} {prepare_sample(config, source_root, sample_id)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.sample_workers) as pool:
            futures = {
                pool.submit(_prepare_one, args.config, str(source_root), sample_id): sample_id
                for sample_id in range(args.sample_start, args.sample_end)
            }
            for future in as_completed(futures):
                sample_id, result = future.result()
                print(f"sample_{sample_id:03d} {result}", flush=True)


if __name__ == "__main__":
    main()
