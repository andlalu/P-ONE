from __future__ import annotations

import contextlib
import contextvars
import hashlib
import importlib.metadata
import json
import logging
import os
import platform
import shutil
import subprocess
import sys
import time
import traceback
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

THREAD_ENV_KEYS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def set_thread_env() -> None:
    for key in THREAD_ENV_KEYS:
        os.environ.setdefault(key, "1")


set_thread_env()

import numpy as np

from DGPSimulation.heston_simulator import HestonPathSimulator
from DGPSimulation.io import load_heston_path_npz, save_heston_path_npz
from DGPSimulation.variance_drawers import AndersenQeVarianceDrawer
from Estimation.ISCGMM.estimate import estimate_first_step
from Estimation.ISCGMM.implied_state import imply_heston_variance_path
from OptionData.add_noise import generate_noisy_panel_rows
from OptionData.clean_panel import generate_clean_option_panel_rows, option_type_for_log_moneyness
from OptionData.io import load_option_panel, parquet_available
from OptionData.noise_common import NOISE_SCENARIOS, price_bounds, scenario_seed
from OptionData.noise_variance_linked import variance_linked_gamma
from OptionPricing.cos_basis import cos_specification_metadata
from OptionPricing.cos_pricer import CosOptionPricer
from Scripts.experiment_config import ExperimentConfig

LOGGER = logging.getLogger(__name__)

FORMAT_VERSION = 2
SCENARIO_ORDER = ("clean",) + NOISE_SCENARIOS
NOISE_COLUMNS = (
    "noise_draw",
    "raw_noisy_iv",
    "raw_noisy_price",
    "observed_price",
    "observed_iv",
    "was_price_capped",
    "cap_direction",
    "noise_seed",
)
TICK_COLUMNS = ("raw_price_before_rounding", "price_after_rounding")
PERSISTENT_FACTOR_COLUMNS = (
    "persistent_factor_level",
    "persistent_factor_moneyness",
    "persistent_factor_maturity",
)
VARIANCE_LINKED_FACTOR_COLUMN = "variance_linked_factor"
COMBINED_PANEL_COLUMNS = (
    "run_id",
    "sample_id",
    "scenario",
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
    "estimation_price",
    "estimation_iv",
    *NOISE_COLUMNS,
    *PERSISTENT_FACTOR_COLUMNS,
    VARIANCE_LINKED_FACTOR_COLUMN,
)


def panel_columns(config: ExperimentConfig) -> tuple[str, ...]:
    """Keep run_002's schema unchanged; add tick diagnostics only when configured."""
    if config.noise is None or config.noise.tick_size is None:
        return COMBINED_PANEL_COLUMNS
    index = COMBINED_PANEL_COLUMNS.index("raw_noisy_price") + 1
    return COMBINED_PANEL_COLUMNS[:index] + TICK_COLUMNS + COMBINED_PANEL_COLUMNS[index:]


def panel_format_version(config: ExperimentConfig) -> int:
    return 3 if config.noise is not None and config.noise.tick_size is not None else FORMAT_VERSION

_LOG_STAGE = contextvars.ContextVar("sample_run_stage", default="initialising")
_LOG_SCENARIO = contextvars.ContextVar("sample_run_scenario", default="-")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def atomic_json(path: str | Path, payload: dict[str, Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f"{target.stem}.tmp{target.suffix}")
    with temporary.open("w", encoding="utf-8") as file_handle:
        json.dump(payload, file_handle, indent=2, sort_keys=True, allow_nan=False)
        file_handle.write("\n")
        file_handle.flush()
        os.fsync(file_handle.fileno())
    os.replace(temporary, target)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as file_handle:
        for block in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_directory(output_root: str | Path, sample_id: int) -> Path:
    if sample_id < 0:
        raise ValueError("sample_id must be non-negative")
    return Path(output_root).expanduser().resolve() / f"sample_{sample_id:03d}"


def _git_state() -> tuple[str, bool]:
    repository = Path(__file__).resolve().parents[2]
    try:
        sha = subprocess.run(
            ["git", "-C", str(repository), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(repository), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        return sha, bool(status.strip())
    except (OSError, subprocess.CalledProcessError):
        return "unavailable", True


def _package_version(distribution: str) -> str:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def _run_payload(config: ExperimentConfig) -> dict[str, Any]:
    git_sha, git_dirty = _git_state()
    return {
        "format_version": panel_format_version(config),
        "run_id": config.run_id,
        "git_sha": git_sha,
        "git_dirty": git_dirty,
        "configuration_hash": config.experiment_config_hash,
        "configuration": config.raw_config,
        "environment": {
            "python": platform.python_version(),
            "numpy": _package_version("numpy"),
            "scipy": _package_version("scipy"),
            "pandas": _package_version("pandas"),
            "pyarrow": _package_version("pyarrow"),
        },
        "platform": platform.platform(),
        "created_at_utc": utc_now(),
    }


def _verify_run_payload(
    payload: dict[str, Any],
    config: ExperimentConfig,
    *,
    allow_git_mismatch: bool,
) -> None:
    if payload.get("run_id") != config.run_id:
        raise ValueError("run.json run_id does not match the experiment configuration")
    if payload.get("configuration_hash") != config.experiment_config_hash:
        raise ValueError("run.json configuration hash does not match the experiment configuration")
    current_sha, _ = _git_state()
    if payload.get("git_sha") != current_sha and not allow_git_mismatch:
        raise ValueError(
            "run.json Git SHA does not match the current checkout; "
            "use --allow-dirty only for an intentional development override"
        )


def _read_json(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as file_handle:
        payload = json.load(file_handle)
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def initialise_or_verify_run(
    output_root: str | Path,
    config: ExperimentConfig,
    *,
    allow_git_mismatch: bool = False,
) -> dict[str, Any]:
    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    target = root / "run.json"
    if target.exists():
        payload = _read_json(target)
        _verify_run_payload(payload, config, allow_git_mismatch=allow_git_mismatch)
        return payload

    lock = root / ".run.json.lock"
    deadline = time.monotonic() + 30.0
    lock_descriptor: int | None = None
    while lock_descriptor is None:
        try:
            lock_descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            if target.exists():
                payload = _read_json(target)
                _verify_run_payload(payload, config, allow_git_mismatch=allow_git_mismatch)
                return payload
            if time.monotonic() >= deadline:
                raise TimeoutError(f"timed out waiting for concurrent run initialisation under {root}")
            time.sleep(0.05)

    try:
        os.close(lock_descriptor)
        if target.exists():
            payload = _read_json(target)
        else:
            payload = _run_payload(config)
            atomic_json(target, payload)
        _verify_run_payload(payload, config, allow_git_mismatch=allow_git_mismatch)
        return payload
    finally:
        lock.unlink(missing_ok=True)


def _initial_record(
    config: ExperimentConfig,
    sample_id: int,
    *,
    generation_only: bool,
    state_only: bool,
    scenario: str | None,
    max_dates: int | None,
) -> dict[str, Any]:
    now = utc_now()
    noise_seeds = {
        name: scenario_seed(config.noise.base_seed, sample_id, name)
        for name in NOISE_SCENARIOS
    } if config.noise is not None else {}
    git_sha, _ = _git_state()
    return {
        "format_version": panel_format_version(config),
        "run_id": config.run_id,
        "sample_id": sample_id,
        "configuration_hash": config.experiment_config_hash,
        "git_sha": git_sha,
        "status": "initialising",
        "current_stage": "initialising",
        "started_at_utc": now,
        "updated_at_utc": now,
        "completed_at_utc": None,
        "execution": {
            "generation_only": generation_only,
            "state_only": state_only,
            "scenario": scenario,
            "max_dates": max_dates,
        },
        "seeds": {
            "path": config.base_seed + sample_id,
            **noise_seeds,
        },
        "artifacts": {},
        "validation": {},
        "estimation": {},
        "diagnostics": {},
        "timings_seconds": {
            "simulation": 0.0,
            "panel_generation": 0.0,
            "validation": 0.0,
            "estimation_clean": 0.0,
            "estimation_low_iid": 0.0,
            "estimation_spatial_corr": 0.0,
            "estimation_persistent_factor": 0.0,
            "estimation_variance_linked_factor": 0.0,
            "total": 0.0,
        },
        "errors": [],
    }


def _publish_record(path: Path, record: dict[str, Any]) -> None:
    record["updated_at_utc"] = utc_now()
    atomic_json(path, record)


def _set_stage(
    record: dict[str, Any],
    record_path: Path,
    stage: str,
    *,
    scenario: str | None = None,
) -> None:
    record["status"] = stage
    record["current_stage"] = stage
    _LOG_STAGE.set(stage)
    _LOG_SCENARIO.set(scenario or "-")
    _publish_record(record_path, record)


class _ContextFilter(logging.Filter):
    def __init__(self, sample_id: int) -> None:
        super().__init__()
        self.sample_id = sample_id

    def filter(self, record: logging.LogRecord) -> bool:
        record.sample_id = self.sample_id
        record.stage = _LOG_STAGE.get()
        record.scenario = _LOG_SCENARIO.get()
        return True


@contextlib.contextmanager
def sample_logging(sample_dir: Path, sample_id: int, log_level: str) -> Iterator[None]:
    root_logger = logging.getLogger()
    level = getattr(logging, log_level)
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s sample=%(sample_id)03d "
        "stage=%(stage)s scenario=%(scenario)s %(name)s: %(message)s"
    )
    context_filter = _ContextFilter(sample_id)
    stream_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler(sample_dir / "sample.log", mode="a", encoding="utf-8")
    for handler in (stream_handler, file_handler):
        handler.setLevel(level)
        handler.setFormatter(formatter)
        handler.addFilter(context_filter)
        root_logger.addHandler(handler)
    old_level = root_logger.level
    root_logger.setLevel(min(old_level, level) if old_level else level)
    try:
        yield
    finally:
        for handler in (stream_handler, file_handler):
            root_logger.removeHandler(handler)
            handler.close()
        root_logger.setLevel(old_level)


def _save_path_atomically(
    target: Path,
    *,
    path: Any,
    parameters: Any,
    simulation: Any,
) -> None:
    temporary = target.with_name("path.tmp.npz")
    save_heston_path_npz(temporary, path=path, params=parameters, config=simulation)
    loaded_path, _, loaded_config = load_heston_path_npz(temporary)
    if (
        len(loaded_path.t_week) != loaded_config.t_week + 1
        or len(loaded_path.dlogS_week) != loaded_config.t_week
    ):
        temporary.unlink(missing_ok=True)
        raise RuntimeError("atomic path validation failed before publication")
    os.replace(temporary, target)


def simulate_and_write_path(
    config: ExperimentConfig,
    sample_id: int,
    target: Path,
) -> tuple[Any, Any]:
    seed = config.base_seed + sample_id
    simulation = replace(config.simulation, seed=seed)
    path = HestonPathSimulator(
        params=config.dgp,
        config=simulation,
        variance_drawer=AndersenQeVarianceDrawer(),
    ).simulate()
    _save_path_atomically(
        target,
        path=path,
        parameters=config.dgp,
        simulation=simulation,
    )
    return path, config.dgp


def _combined_clean_row(row: dict[str, Any]) -> dict[str, Any]:
    item = dict(row)
    item.update(
        {
            "scenario": "clean",
            "estimation_price": float(row["model_price"]),
            "estimation_iv": float(row["model_iv"]),
            "noise_draw": None,
            "raw_noisy_iv": None,
            "raw_noisy_price": None,
            "observed_price": None,
            "observed_iv": None,
            "was_price_capped": False,
            "cap_direction": "none",
            "noise_seed": None,
            "persistent_factor_level": None,
            "persistent_factor_moneyness": None,
            "persistent_factor_maturity": None,
            "variance_linked_factor": None,
        }
    )
    return item


def _combined_noisy_row(
    row: dict[str, Any],
    scenario: str,
    factor: dict[str, Any] | None,
) -> dict[str, Any]:
    item = dict(row)
    item.pop("noise_scenario", None)
    item.update(
        {
            "scenario": scenario,
            "estimation_price": float(row["observed_price"]),
            "estimation_iv": float(row["observed_iv"]),
            "persistent_factor_level": None if factor is None else float(factor["factor_0"]),
            "persistent_factor_moneyness": None if factor is None else float(factor["factor_1"]),
            "persistent_factor_maturity": None if factor is None else float(factor["factor_2"]),
            "variance_linked_factor": (
                None
                if factor is None or factor.get("variance_linked_factor") is None
                else float(factor["variance_linked_factor"])
            ),
        }
    )
    return item


def _write_combined_panel(
    rows: list[dict[str, Any]],
    target: Path,
    *,
    config: ExperimentConfig,
    sample_id: int,
    git_sha: str,
) -> None:
    if not parquet_available():
        raise RuntimeError("combined production panel output requires pandas and pyarrow")
    import pandas as pd  # type: ignore[import-not-found]
    import pyarrow as pa  # type: ignore[import-not-found]
    import pyarrow.parquet as pq  # type: ignore[import-not-found]

    columns = panel_columns(config)
    frame = pd.DataFrame(rows, columns=columns)
    frame["noise_seed"] = pd.array(frame["noise_seed"], dtype="Int64")
    table = pa.Table.from_pandas(frame, preserve_index=False)
    metadata = dict(table.schema.metadata or {})
    metadata.update(
        {
            b"p_one.run_id": config.run_id.encode("utf-8"),
            b"p_one.sample_id": str(sample_id).encode("ascii"),
            b"p_one.configuration_hash": config.experiment_config_hash.encode("ascii"),
            b"p_one.git_sha": git_sha.encode("ascii"),
            b"p_one.format_version": str(panel_format_version(config)).encode("ascii"),
            b"p_one.scenario_order": json.dumps(SCENARIO_ORDER).encode("utf-8"),
            b"p_one.cos_basis": json.dumps(
                cos_specification_metadata(config.cos_basis),
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8"),
        }
    )
    table = table.replace_schema_metadata(metadata)
    temporary = target.with_name("panels.tmp.parquet")
    pq.write_table(table, temporary)
    published = pq.read_table(temporary)
    if published.num_rows != len(rows):
        temporary.unlink(missing_ok=True)
        raise RuntimeError("atomic Parquet row-count validation failed before publication")
    if published.column_names != list(columns):
        temporary.unlink(missing_ok=True)
        raise RuntimeError("atomic Parquet schema validation failed before publication")
    os.replace(temporary, target)


def build_combined_panel(
    config: ExperimentConfig,
    sample_id: int,
    path: Any,
    parameters: Any,
    target: Path,
) -> dict[str, Any]:
    if config.noise is None or config.noise.scenario_names() != NOISE_SCENARIOS:
        raise ValueError("the sample runner requires all four production noise scenarios")
    clean_rows = generate_clean_option_panel_rows(
        run_id=config.run_id,
        sample_id=sample_id,
        path=path,
        params_p=parameters,
        eta_v=config.eta_v,
        maturities_years=config.cos_basis.maturities,
        log_moneyness=config.log_moneyness,
        atm_option_type=config.atm_option_type,
        pricing_method="COS",
        iv_method="lets_be_rational",
        pricer=CosOptionPricer(),
        cos_basis=config.cos_basis,
    )
    combined = [_combined_clean_row(row) for row in clean_rows]
    summaries: dict[str, Any] = {
        "clean": {"rows": len(clean_rows), "seed": config.base_seed + sample_id}
    }
    for scenario in NOISE_SCENARIOS:
        seed = scenario_seed(config.noise.base_seed, sample_id, scenario)
        noisy_rows, factors = generate_noisy_panel_rows(
            clean_rows,
            scenario=scenario,
            seed=seed,
            config=config.noise,
            params_p=parameters,
        )
        factor_by_week = {int(item["week_index"]): item for item in factors}
        combined.extend(
            _combined_noisy_row(
                row,
                scenario,
                (
                    factor_by_week.get(int(row["week_index"]))
                    if scenario in {"persistent_factor", "variance_linked_factor"}
                    else None
                ),
            )
            for row in noisy_rows
        )
        lower = sum(str(row["cap_direction"]) == "lower" for row in noisy_rows)
        upper = sum(str(row["cap_direction"]) == "upper" for row in noisy_rows)
        summaries[scenario] = {
            "rows": len(noisy_rows),
            "seed": seed,
            "n_capped_lower": lower,
            "n_capped_upper": upper,
            "n_capped_total": lower + upper,
        }
    git_sha, _ = _git_state()
    _write_combined_panel(
        combined,
        target,
        config=config,
        sample_id=sample_id,
        git_sha=git_sha,
    )
    return summaries


def validate_path_artifact(path_file: str | Path, config: ExperimentConfig, sample_id: int) -> dict[str, Any]:
    path, _, stored_config = load_heston_path_npz(path_file)
    arrays = [path.t_week, path.logS_week, path.V_week, path.dlogS_week]
    daily_present = path.logS_daily is not None and path.V_daily is not None
    if path.logS_daily is not None:
        arrays.append(path.logS_daily)
    if path.V_daily is not None:
        arrays.append(path.V_daily)
    checks = {
        "weekly_path_length_valid": len(path.t_week) == stored_config.t_week + 1,
        "weekly_return_length_valid": len(path.dlogS_week) == stored_config.t_week,
        "daily_arrays_present_when_requested": (not stored_config.return_daily) or daily_present,
        "all_path_values_finite": all(np.all(np.isfinite(values)) for values in arrays),
        "weekly_variance_nonnegative": bool(np.all(np.asarray(path.V_week) >= 0.0)),
        "daily_variance_nonnegative": (
            True if path.V_daily is None else bool(np.all(np.asarray(path.V_daily) >= 0.0))
        ),
        "path_seed_matches": stored_config.seed == config.base_seed + sample_id,
    }
    variance = np.asarray(path.V_week, dtype=float)
    spot = np.exp(np.asarray(path.logS_week, dtype=float))
    return {
        **checks,
        "spot_min": float(np.min(spot)),
        "spot_max": float(np.max(spot)),
        "variance_min": float(np.min(variance)),
        "variance_max": float(np.max(variance)),
        "variance_mean": float(np.mean(variance)),
        "variance_std": float(np.std(variance)),
        "passed": all(checks.values()),
    }


def _finite(frame: Any, column: str) -> bool:
    return bool(np.all(np.isfinite(frame[column].astype(float).to_numpy())))


def _panel_price_bounds_valid(frame: Any, config: ExperimentConfig, scenario: str) -> bool:
    assert config.noise is not None
    for row in frame.to_dict("records"):
        lower, upper = price_bounds(
            float(row["S"]),
            float(row["strike"]),
            float(row["maturity_years"]),
            float(row["r"]),
            float(row["q"]),
            str(row["option_type"]),
        )
        price = float(row["estimation_price"])
        if scenario == "clean":
            if price < lower - 1e-10 or price > upper + 1e-10:
                return False
        elif (
            price < lower + config.noise.price_epsilon - 1e-11
            or price > upper - config.noise.price_epsilon + 1e-11
        ):
            return False
    return True


def _scenario_validation(frame: Any, config: ExperimentConfig, sample_id: int, scenario: str) -> dict[str, Any]:
    expected_rows = (
        (config.simulation.t_week + 1)
        * len(config.cos_basis.maturities)
        * len(config.log_moneyness)
    )
    keys = list(
        zip(
            frame["week_index"].astype(int),
            frame["maturity_years"].astype(float),
            frame["log_moneyness"].astype(float),
        )
    )
    row_order_valid = keys == sorted(keys)
    option_rule = all(
        str(row.option_type)
        == option_type_for_log_moneyness(float(row.log_moneyness), config.atm_option_type)
        for row in frame.itertuples(index=False)
    )
    if scenario == "clean":
        seed_valid = bool(frame["noise_seed"].isna().all())
    else:
        assert config.noise is not None
        expected_seed = scenario_seed(config.noise.base_seed, sample_id, scenario)
        seed_valid = bool(
            frame["noise_seed"].notna().all()
            and np.all(frame["noise_seed"].astype(int).to_numpy() == expected_seed)
        )
    checks: dict[str, Any] = {
        "expected_row_count": expected_rows,
        "actual_row_count": int(len(frame)),
        "row_order_valid": row_order_valid,
        "finite_model_prices": _finite(frame, "model_price"),
        "finite_model_ivs": _finite(frame, "model_iv"),
        "finite_estimation_prices": _finite(frame, "estimation_price"),
        "finite_estimation_ivs": _finite(frame, "estimation_iv"),
        "nonnegative_ivs": bool(
            np.all(frame["model_iv"].astype(float).to_numpy() >= 0.0)
            and np.all(frame["estimation_iv"].astype(float).to_numpy() >= 0.0)
        ),
        "otm_option_rule_valid": option_rule,
        "price_bounds_valid": _panel_price_bounds_valid(frame, config, scenario),
        "scenario_seed_valid": seed_valid,
    }
    if scenario == "clean":
        nullable = [name for name in NOISE_COLUMNS if name not in {"was_price_capped", "cap_direction"}]
        checks["clean_noise_fields_nullable_valid"] = bool(
            all(frame[name].isna().all() for name in nullable)
            and not frame["was_price_capped"].astype(bool).any()
            and (frame["cap_direction"].astype(str) == "none").all()
        )
        if config.noise is not None and config.noise.tick_size is not None:
            checks["clean_tick_fields_nullable_valid"] = bool(
                all(frame[name].isna().all() for name in TICK_COLUMNS)
            )
    else:
        required = [
            "noise_draw",
            "raw_noisy_iv",
            "raw_noisy_price",
            "observed_price",
            "observed_iv",
            "was_price_capped",
            "cap_direction",
            "noise_seed",
        ]
        checks["required_noise_fields_present"] = bool(
            all(name in frame.columns and frame[name].notna().all() for name in required)
        )
        capped = frame["was_price_capped"].astype(bool).to_numpy()
        directions = frame["cap_direction"].astype(str).to_numpy()
        checks["cap_flag_consistency"] = bool(np.all(capped == (directions != "none")))
        checks["n_capped_lower"] = int(np.count_nonzero(directions == "lower"))
        checks["n_capped_upper"] = int(np.count_nonzero(directions == "upper"))
        checks["n_capped_total"] = int(np.count_nonzero(capped))
        if config.noise is not None and config.noise.tick_size is not None:
            raw = frame["raw_price_before_rounding"].astype(float).to_numpy()
            rounded = frame["price_after_rounding"].astype(float).to_numpy()
            alias = frame["raw_noisy_price"].astype(float).to_numpy()
            expected = config.noise.tick_size * np.rint(raw / config.noise.tick_size)
            checks["tick_rounding_valid"] = bool(
                np.isfinite(raw).all()
                and np.isfinite(rounded).all()
                and np.array_equal(raw, alias)
                and np.allclose(rounded, expected, rtol=0.0, atol=1e-12)
            )
    if scenario in {"persistent_factor", "variance_linked_factor"}:
        factor_frame = frame[["week_index", *PERSISTENT_FACTOR_COLUMNS]]
        checks["factor_rows_cover_all_dates"] = (
            factor_frame["week_index"].nunique() == config.simulation.t_week + 1
        )
        checks["factor_values_constant_within_date"] = bool(
            (
                factor_frame.groupby("week_index", sort=True)[list(PERSISTENT_FACTOR_COLUMNS)]
                .nunique(dropna=False)
                .to_numpy()
                == 1
            ).all()
        )
        checks["factor_values_finite"] = bool(
            np.all(
                np.isfinite(
                    factor_frame[list(PERSISTENT_FACTOR_COLUMNS)].astype(float).to_numpy()
                )
            )
        )
    else:
        checks["persistent_factor_values_nullable_valid"] = bool(
            all(frame[name].isna().all() for name in PERSISTENT_FACTOR_COLUMNS)
        )
    if scenario == "variance_linked_factor":
        variance_factor = frame[["week_index", VARIANCE_LINKED_FACTOR_COLUMN]]
        checks["variance_factor_rows_cover_all_dates"] = (
            variance_factor["week_index"].nunique() == config.simulation.t_week + 1
        )
        checks["variance_factor_values_constant_within_date"] = bool(
            (
                variance_factor.groupby("week_index", sort=True)[
                    VARIANCE_LINKED_FACTOR_COLUMN
                ]
                .nunique(dropna=False)
                .to_numpy()
                == 1
            ).all()
        )
        checks["variance_factor_values_finite"] = _finite(
            variance_factor,
            VARIANCE_LINKED_FACTOR_COLUMN,
        )
        assert config.noise is not None
        factor_config = config.noise.scenarios[scenario]
        gamma_v = variance_linked_gamma(
            frame["log_moneyness"].astype(float).to_numpy(),
            frame["maturity_years"].astype(float).to_numpy(),
            factor_config,
        )
        variance_state_std = np.sqrt(
            config.dgp.vbar
            * config.dgp.sigma_v
            * config.dgp.sigma_v
            / (2.0 * config.dgp.kappa)
        )
        expected = gamma_v * (
            frame["V"].astype(float).to_numpy() - config.dgp.vbar
        ) / variance_state_std
        actual = frame[VARIANCE_LINKED_FACTOR_COLUMN].astype(float).to_numpy()
        checks["variance_factor_matches_latent_variance"] = bool(
            np.allclose(actual, expected, rtol=1e-12, atol=1e-15)
        )
    else:
        checks["variance_linked_factor_nullable_valid"] = bool(
            frame[VARIANCE_LINKED_FACTOR_COLUMN].isna().all()
        )
    boolean_checks = [
        value
        for name, value in checks.items()
        if isinstance(value, (bool, np.bool_))
    ]
    checks["passed"] = (
        int(checks["actual_row_count"]) == int(checks["expected_row_count"])
        and all(bool(value) for value in boolean_checks)
    )
    return checks


def validate_panel_artifact(
    panel_file: str | Path,
    config: ExperimentConfig,
    sample_id: int,
) -> dict[str, Any]:
    import pandas as pd  # type: ignore[import-not-found]

    frame = pd.read_parquet(panel_file)
    missing = set(panel_columns(config)) - set(frame.columns)
    if missing:
        raise AssertionError(f"combined panel is missing columns: {sorted(missing)}")
    sample_ids = set(frame["sample_id"].astype(int))
    scenario_values = tuple(frame["scenario"].drop_duplicates().astype(str))
    scenario_rank = {name: index for index, name in enumerate(SCENARIO_ORDER)}
    combined_order = [
        (
            scenario_rank.get(str(row.scenario), len(SCENARIO_ORDER)),
            int(row.week_index),
            float(row.maturity_years),
            float(row.log_moneyness),
        )
        for row in frame.itertuples(index=False)
    ]
    scenarios_valid = (
        scenario_values == SCENARIO_ORDER
        and combined_order == sorted(combined_order)
    )
    sample_valid = sample_ids == {sample_id}
    results = {
        scenario: _scenario_validation(
            frame[frame["scenario"].astype(str) == scenario].reset_index(drop=True),
            config,
            sample_id,
            scenario,
        )
        for scenario in SCENARIO_ORDER
    }
    return {
        "sample_id_valid": sample_valid,
        "scenario_order_valid": scenarios_valid,
        "scenarios": results,
        "passed": sample_valid and scenarios_valid and all(item["passed"] for item in results.values()),
    }


def validate_sample_artifacts(
    path_file: str | Path,
    panel_file: str | Path,
    config: ExperimentConfig,
    sample_id: int,
) -> dict[str, Any]:
    path_result = validate_path_artifact(path_file, config, sample_id)
    panel_result = validate_panel_artifact(panel_file, config, sample_id)
    return {
        "path": path_result,
        "panels": panel_result,
        "passed": bool(path_result["passed"] and panel_result["passed"]),
    }


def _verify_recorded_artifact(record: dict[str, Any], name: str, path: Path) -> bool:
    artifact = record.get("artifacts", {}).get(name)
    if artifact is None:
        return False
    if not path.exists():
        return False
    expected_hash = artifact.get("sha256")
    if not expected_hash:
        return False
    actual_hash = sha256_file(path)
    if actual_hash != expected_hash:
        raise ValueError(f"{name} artifact checksum does not match record.json")
    if int(artifact.get("bytes", -1)) != path.stat().st_size:
        raise ValueError(f"{name} artifact size does not match record.json")
    return True


def _record_artifacts(
    record: dict[str, Any],
    path_file: Path,
    panel_file: Path,
    validation: dict[str, Any],
) -> None:
    rows_by_scenario = {
        scenario: int(result["actual_row_count"])
        for scenario, result in validation["panels"]["scenarios"].items()
    }
    record["artifacts"] = {
        "path": {
            "file": "path.npz",
            "sha256": sha256_file(path_file),
            "bytes": path_file.stat().st_size,
        },
        "panels": {
            "file": "panels.parquet",
            "sha256": sha256_file(panel_file),
            "bytes": panel_file.stat().st_size,
            "rows": sum(rows_by_scenario.values()),
            "rows_by_scenario": rows_by_scenario,
        },
    }


def _record_matches_execution(
    record: dict[str, Any],
    config: ExperimentConfig,
    sample_id: int,
) -> None:
    if record.get("sample_id") != sample_id:
        raise ValueError("record.json sample_id does not match the requested sample")
    if record.get("run_id") != config.run_id:
        raise ValueError("record.json run_id does not match the experiment configuration")
    if record.get("configuration_hash") != config.experiment_config_hash:
        raise ValueError("record.json configuration hash does not match the experiment configuration")


def _estimate_is_complete(result: Any) -> bool:
    if not isinstance(result, dict) or not isinstance(result.get("success"), bool):
        return False

    required = (
        "final_criterion",
        "estimated_parameters",
        "free_parameters",
        "function_evaluations",
        "penalty_evaluations",
        "final_diagnostics",
        "total_runtime_seconds",
    )
    if any(name not in result for name in required):
        return False
    try:
        if not np.isfinite(float(result["final_criterion"])):
            return False
        free_parameters = np.asarray(result["free_parameters"], dtype=float)
        if free_parameters.shape != (6,) or not np.all(np.isfinite(free_parameters)):
            return False
        required_parameters = {
            "eta",
            "kappa",
            "vbar",
            "sigma_v",
            "rho",
            "eta_v",
            "r",
            "q",
        }
        if set(result["estimated_parameters"]) != required_parameters:
            return False
        if not all(
            np.isfinite(float(value))
            for value in result["estimated_parameters"].values()
        ):
            return False
        if int(result["function_evaluations"]) < 0:
            return False
        if int(result["penalty_evaluations"]) < 0:
            return False
        if (
            not np.isfinite(float(result["total_runtime_seconds"]))
            or float(result["total_runtime_seconds"]) < 0.0
        ):
            return False
    except (AttributeError, TypeError, ValueError):
        return False
    if not isinstance(result["final_diagnostics"], dict) or not result["final_diagnostics"]:
        return False

    passes = result.get("powell_passes")
    return (
        isinstance(passes, list)
        and len(passes) == 2
        and all(isinstance(item, dict) for item in passes)
        and [item.get("pass_name") for item in passes] == ["coarse", "refinement"]
    )


def estimate_scenario(
    panel_file: str | Path,
    scenario: str,
    config: ExperimentConfig,
    *,
    max_dates: int | None,
) -> dict[str, Any]:
    panel = load_option_panel(
        panel_file,
        scenario=scenario,
        max_dates=max_dates,
    )
    return estimate_first_step(
        panel,
        criterion_config=config.criterion_config,
        powell_config=config.powell_config,
    ).to_dict()


def _state_only_diagnostic(
    panel_file: str | Path,
    scenario: str,
    config: ExperimentConfig,
    *,
    max_dates: int | None,
) -> dict[str, Any]:
    panel = load_option_panel(panel_file, scenario=scenario, max_dates=max_dates)
    result = imply_heston_variance_path(
        config.powell_config.base_start,
        panel,
        config.criterion_config.implied_state,
    )
    true_variance = panel.true_variance
    state_rmse = None
    if true_variance is not None:
        state_rmse = float(np.sqrt(np.mean((result.variance - true_variance) ** 2)))
    return {
        "base_start": asdict(config.powell_config.base_start),
        "state_rmse": state_rmse,
        "implied_state": result.to_dict(),
    }


def _prepare_sample_directory(
    output_root: str | Path,
    sample_id: int,
    *,
    resume: bool,
    overwrite: bool,
) -> Path:
    if resume and overwrite:
        raise ValueError("--resume and --overwrite cannot be combined")
    directory = sample_directory(output_root, sample_id)
    if overwrite and directory.exists():
        shutil.rmtree(directory)
    elif directory.exists() and any(directory.iterdir()) and not resume:
        raise FileExistsError(
            f"{directory} already contains outputs; use --resume or --overwrite"
        )
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def run_sample(
    *,
    config: ExperimentConfig,
    sample_id: int,
    output_root: str | Path,
    resume: bool = False,
    overwrite: bool = False,
    generation_only: bool = False,
    state_only: bool = False,
    scenario: str | None = None,
    max_dates: int | None = None,
    log_level: str = "INFO",
    allow_dirty: bool = False,
) -> dict[str, Any]:
    if sample_id >= config.n_samples:
        raise ValueError(f"sample_id must be below configured n_samples={config.n_samples}")
    if generation_only and state_only:
        raise ValueError("--generation-only and --state-only cannot be combined")
    if scenario is not None and scenario not in SCENARIO_ORDER:
        raise ValueError(f"unknown scenario: {scenario}")
    if max_dates is not None and max_dates <= 0:
        raise ValueError("--max-dates must be positive")
    if config.panel_format != "parquet":
        raise ValueError("the sample-centred production runner requires panel_format='parquet'")

    set_thread_env()
    sample_dir = _prepare_sample_directory(
        output_root,
        sample_id,
        resume=resume,
        overwrite=overwrite,
    )
    record_path = sample_dir / "record.json"
    path_file = sample_dir / "path.npz"
    panel_file = sample_dir / "panels.parquet"
    started = time.perf_counter()
    record: dict[str, Any] | None = None

    with sample_logging(sample_dir, sample_id, log_level):
        try:
            initialise_or_verify_run(
                output_root,
                config,
                allow_git_mismatch=allow_dirty,
            )
            if resume and record_path.exists():
                record = _read_json(record_path)
                _record_matches_execution(record, config, sample_id)
                record["execution"] = {
                    "generation_only": generation_only,
                    "state_only": state_only,
                    "scenario": scenario,
                    "max_dates": max_dates,
                }
            else:
                record = _initial_record(
                    config,
                    sample_id,
                    generation_only=generation_only,
                    state_only=state_only,
                    scenario=scenario,
                    max_dates=max_dates,
                )
                _publish_record(record_path, record)

            if resume and record.get("status") == "complete":
                path_valid = _verify_recorded_artifact(record, "path", path_file)
                panels_valid = _verify_recorded_artifact(record, "panels", panel_file)
                if path_valid and panels_valid:
                    LOGGER.info("complete sample verified; no computation required")
                    return record

            _set_stage(record, record_path, "generating")
            path_reused = False
            if resume and path_file.exists():
                _verify_recorded_artifact(record, "path", path_file)
                if validate_path_artifact(path_file, config, sample_id)["passed"]:
                    path, parameters, _ = load_heston_path_npz(path_file)
                    path_reused = True
                    LOGGER.info("path reused")
            if not path_reused:
                stage_started = time.perf_counter()
                path, parameters = simulate_and_write_path(config, sample_id, path_file)
                record["timings_seconds"]["simulation"] = time.perf_counter() - stage_started
                record["artifacts"] = {}
                record["validation"] = {}
                record["estimation"] = {}
                LOGGER.info(
                    "path completed elapsed=%.3fs",
                    record["timings_seconds"]["simulation"],
                )

            panels_reused = False
            if resume and panel_file.exists():
                _verify_recorded_artifact(record, "panels", panel_file)
                if validate_panel_artifact(panel_file, config, sample_id)["passed"]:
                    panels_reused = True
                    LOGGER.info("combined panels reused")
            if not panels_reused:
                stage_started = time.perf_counter()
                summaries = build_combined_panel(
                    config,
                    sample_id,
                    path,
                    parameters,
                    panel_file,
                )
                record["timings_seconds"]["panel_generation"] = (
                    time.perf_counter() - stage_started
                )
                record["artifacts"].pop("panels", None)
                record["validation"] = {}
                record["estimation"] = {}
                LOGGER.info(
                    "combined panels completed elapsed=%.3fs rows=%d summaries=%s",
                    record["timings_seconds"]["panel_generation"],
                    sum(item["rows"] for item in summaries.values()),
                    summaries,
                )
            _publish_record(record_path, record)

            _set_stage(record, record_path, "validating")
            stage_started = time.perf_counter()
            validation = validate_sample_artifacts(path_file, panel_file, config, sample_id)
            record["timings_seconds"]["validation"] = time.perf_counter() - stage_started
            record["validation"] = validation
            if validation["passed"]:
                _record_artifacts(record, path_file, panel_file, validation)
            _publish_record(record_path, record)
            if not validation["passed"]:
                raise AssertionError("sample artifact validation failed")
            LOGGER.info(
                "path and panels passed rows=%d elapsed=%.3fs",
                record["artifacts"]["panels"]["rows"],
                record["timings_seconds"]["validation"],
            )

            if generation_only:
                LOGGER.info("generation-only run stopped after successful validation")
                record["timings_seconds"]["total"] += time.perf_counter() - started
                record["status"] = "generated"
                record["current_stage"] = "validated"
                _LOG_STAGE.set("validated")
                _publish_record(record_path, record)
                return record

            requested_scenarios = (scenario,) if scenario is not None else SCENARIO_ORDER
            _set_stage(record, record_path, "estimating")
            if state_only:
                diagnostic_scenario = scenario or "clean"
                _LOG_SCENARIO.set(diagnostic_scenario)
                record["diagnostics"].setdefault("state_only", {})[diagnostic_scenario] = (
                    _state_only_diagnostic(
                        panel_file,
                        diagnostic_scenario,
                        config,
                        max_dates=max_dates,
                    )
                )
                record["timings_seconds"]["total"] += time.perf_counter() - started
                _publish_record(record_path, record)
                LOGGER.info("state-only diagnostic completed")
                return record

            for name in requested_scenarios:
                _LOG_SCENARIO.set(name)
                if (
                    resume
                    and name in record["estimation"]
                    and _estimate_is_complete(record["estimation"][name])
                ):
                    LOGGER.info("completed estimate reused")
                    continue
                record["estimation"].pop(name, None)
                estimate_started = time.perf_counter()
                result = estimate_scenario(
                    panel_file,
                    name,
                    config,
                    max_dates=max_dates,
                )
                record["estimation"][name] = result
                record["timings_seconds"][f"estimation_{name}"] = (
                    time.perf_counter() - estimate_started
                )
                _publish_record(record_path, record)
                LOGGER.info(
                    "estimate completed success=%s elapsed=%.3fs",
                    result["success"],
                    record["timings_seconds"][f"estimation_{name}"],
                )
                if not result["success"]:
                    LOGGER.warning(
                        "bounded Powell estimate retained with success=false status=%s message=%s",
                        result.get("status"),
                        result.get("message"),
                    )

            unrestricted = scenario is None and max_dates is None
            all_estimated = all(name in record["estimation"] for name in SCENARIO_ORDER)
            record["timings_seconds"]["total"] += time.perf_counter() - started
            if unrestricted and all_estimated:
                record["status"] = "complete"
                record["current_stage"] = "complete"
                record["completed_at_utc"] = utc_now()
                _LOG_STAGE.set("complete")
                _LOG_SCENARIO.set("-")
                LOGGER.info("sample completed")
            _publish_record(record_path, record)
            return record
        except Exception as exception:
            _LOG_SCENARIO.set(_LOG_SCENARIO.get())
            LOGGER.exception("sample run failed")
            if record is not None:
                failing_stage = str(record.get("current_stage", _LOG_STAGE.get()))
                record["status"] = "failed"
                record["current_stage"] = failing_stage
                record["completed_at_utc"] = None
                record["timings_seconds"]["total"] += time.perf_counter() - started
                record["errors"].append(
                    {
                        "stage": failing_stage,
                        "scenario": None if _LOG_SCENARIO.get() == "-" else _LOG_SCENARIO.get(),
                        "exception_type": type(exception).__name__,
                        "message": str(exception),
                        "traceback": traceback.format_exc(),
                        "occurred_at_utc": utc_now(),
                    }
                )
                _publish_record(record_path, record)
            raise
