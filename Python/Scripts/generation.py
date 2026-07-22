from __future__ import annotations

import csv
import json
import logging
import multiprocessing as mp
import os
import platform
import subprocess
import time
import traceback
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

from DGPSimulation.heston_simulator import HestonPathSimulator
from DGPSimulation.io import load_heston_path_npz, save_heston_path_npz
from DGPSimulation.types import HestonSimConfig
from DGPSimulation.variance_drawers import AndersenQeVarianceDrawer
from OptionData.io import panel_metadata_path
from OptionPricing.clean_panel import generate_clean_option_panel_rows, parquet_available, write_panel
from OptionPricing.cos_basis import cos_specification_metadata, validate_panel_cos_compatibility
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.noisy_panel import (
    NoisyPanelResult,
    generate_noisy_panel_file,
    read_table,
    write_noisy_manifest,
)
from Scripts.experiment_config import ExperimentConfig, load_experiment_config

LOGGER = logging.getLogger(__name__)
THREAD_ENV_KEYS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")


@dataclass(frozen=True)
class SampleGenerationResult:
    run_id: str
    sample_id: int
    seed: int
    path_file: str
    panel_file: str
    completion_file: str
    status: str
    elapsed_seconds: float
    noisy_results: tuple[NoisyPanelResult, ...] = ()
    error: str = ""


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w") as file_handle:
        json.dump(payload, file_handle, indent=2, sort_keys=True, allow_nan=False)
        file_handle.write("\n")
        file_handle.flush()
        os.fsync(file_handle.fileno())
    os.replace(temporary, path)


def load_config(config_path: str | Path) -> ExperimentConfig:
    return load_experiment_config(config_path)


def validate_panel_dependencies(panel_format: str) -> None:
    if panel_format == "parquet" and not parquet_available():
        raise RuntimeError("production Parquet output requires pandas and pyarrow or fastparquet")
    if panel_format not in {"parquet", "csv"}:
        raise ValueError("panel_format must be 'parquet' or 'csv'")


def set_thread_env() -> None:
    for key in THREAD_ENV_KEYS:
        os.environ.setdefault(key, "1")


def _git_sha() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def with_overrides(
    config: ExperimentConfig,
    *,
    n_samples: int | None = None,
    workers: int | None = None,
    output_root: str | None = None,
    skip_existing: bool | None = None,
    paths_only: bool | None = None,
    panels_only: bool | None = None,
) -> ExperimentConfig:
    """Apply direct execution overrides without changing the experiment hash."""

    updates: dict[str, Any] = {}
    for name, value in {
        "n_samples": n_samples,
        "workers": workers,
        "skip_existing": skip_existing,
        "paths_only": paths_only,
        "panels_only": panels_only,
    }.items():
        if value is not None:
            updates[name] = value
    if output_root is not None:
        updates["output_root"] = str(Path(output_root).expanduser().resolve())
    result = replace(config, **updates)
    if result.n_samples <= 0:
        raise ValueError("n_samples must be positive")
    if result.workers is not None and result.workers <= 0:
        raise ValueError("workers must be positive")
    if result.paths_only and result.panels_only:
        raise ValueError("paths_only and panels_only cannot both be true")
    return result


def default_workers(n_tasks: int, requested: int | None) -> int:
    if requested is not None:
        return max(1, min(int(requested), n_tasks))
    return min(max((os.cpu_count() or 2) - 1, 1), n_tasks)


def sample_artifact_paths(config: ExperimentConfig, sample_id: int) -> tuple[Path, Path, Path]:
    root = Path(config.output_root)
    stem = f"sample_{sample_id:03d}"
    return (
        root / "paths" / f"{stem}.npz",
        root / "panels_clean" / f"{stem}.{config.panel_format}",
        root / "status" / f"{stem}.complete.json",
    )


def _read_json(path: Path) -> dict[str, Any]:
    with path.open() as file_handle:
        payload = json.load(file_handle)
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def sample_is_complete(config: ExperimentConfig, sample_id: int) -> bool:
    path_file, panel_file, completion_file = sample_artifact_paths(config, sample_id)
    required: list[Path] = []
    if not config.panels_only:
        required.append(path_file)
    if not config.paths_only:
        required.extend([panel_file, panel_metadata_path(panel_file)])
        if config.noise is not None:
            for scenario in config.noise.scenario_names():
                noisy = Path(config.output_root) / "panels_observed" / scenario / f"sample_{sample_id:03d}.{config.panel_format}"
                required.extend([noisy, panel_metadata_path(noisy)])
                if scenario == "persistent_factor":
                    required.append(
                        Path(config.output_root) / "noise_factors" / scenario / f"sample_{sample_id:03d}.csv"
                    )
    required.append(completion_file)
    if not all(path.exists() and path.stat().st_size > 0 for path in required):
        return False
    try:
        completion = _read_json(completion_file)
        if completion.get("status") != "complete":
            return False
        if completion.get("experiment_config_hash") != config.experiment_config_hash:
            return False
        if not config.paths_only:
            metadata = _read_json(panel_metadata_path(panel_file))
            if metadata.get("experiment_config_hash") != config.experiment_config_hash:
                return False
            validate_panel_cos_compatibility(config.cos_basis, metadata["cos_basis"])
            expected_rows = (
                (config.simulation.t_week + 1)
                * len(config.cos_basis.maturities)
                * len(config.log_moneyness)
            )
            if (
                int(completion.get("clean_panel_rows", -1)) != expected_rows
                or len(read_table(panel_file)) != expected_rows
            ):
                return False
        if not config.panels_only:
            path, _, stored_config = load_heston_path_npz(path_file)
            if len(path.t_week) != stored_config.t_week + 1 or len(path.dlogS_week) != stored_config.t_week:
                return False
    except (KeyError, OSError, ValueError, json.JSONDecodeError):
        return False
    return True


def _save_path_atomically(path_file: Path, path: Any, params: Any, simulation: HestonSimConfig) -> None:
    temporary = path_file.with_name(path_file.stem + ".tmp.npz")
    save_heston_path_npz(temporary, path=path, params=params, config=simulation)
    loaded_path, _, loaded_config = load_heston_path_npz(temporary)
    if len(loaded_path.t_week) != loaded_config.t_week + 1 or len(loaded_path.dlogS_week) != loaded_config.t_week:
        temporary.unlink(missing_ok=True)
        raise RuntimeError("atomic path validation failed before publication")
    os.replace(temporary, path_file)


def generate_one_sample(sample_id: int, config: ExperimentConfig) -> SampleGenerationResult:
    started = time.perf_counter()
    seed = config.base_seed + sample_id
    path_file, panel_file, completion_file = sample_artifact_paths(config, sample_id)
    noisy_results: list[NoisyPanelResult] = []
    try:
        if config.skip_existing and sample_is_complete(config, sample_id):
            return SampleGenerationResult(
                config.run_id,
                sample_id,
                seed,
                str(path_file),
                str(panel_file if not config.paths_only else ""),
                str(completion_file),
                "skipped",
                time.perf_counter() - started,
            )
        simulation = replace(config.simulation, seed=seed)
        if config.panels_only:
            path, parameters, _ = load_heston_path_npz(path_file)
        else:
            path = HestonPathSimulator(
                params=config.dgp,
                config=simulation,
                variance_drawer=AndersenQeVarianceDrawer(),
            ).simulate()
            parameters = config.dgp
            _save_path_atomically(path_file, path, parameters, simulation)

        clean_rows = 0
        if not config.paths_only:
            rows = generate_clean_option_panel_rows(
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
            clean_rows = len(rows)
            panel_file = write_panel(
                rows,
                panel_file.with_suffix(""),
                panel_format=config.panel_format,
                metadata={
                    "run_id": config.run_id,
                    "sample_id": sample_id,
                    "seed": seed,
                    "scenario": "clean",
                    "panel_format": config.panel_format,
                    "experiment_config_hash": config.experiment_config_hash,
                    "git_sha": _git_sha(),
                    "cos_basis": cos_specification_metadata(config.cos_basis),
                },
            )
            if config.noise is not None:
                for scenario in config.noise.scenario_names():
                    result = generate_noisy_panel_file(
                        clean_panel_path=panel_file,
                        run_root=config.output_root,
                        sample_id=sample_id,
                        scenario=scenario,
                        config=config.noise,
                        panel_format=config.panel_format,
                        skip_existing=False,
                    )
                    noisy_results.append(result)
                    if result.status == "error":
                        raise RuntimeError(f"{scenario}: {result.error_message}")
        _atomic_json(
            completion_file,
            {
                "sample_id": sample_id,
                "seed": seed,
                "experiment_config_hash": config.experiment_config_hash,
                "git_sha": _git_sha(),
                "clean_panel_rows": clean_rows,
                "generated_scenarios": [result.noise_scenario for result in noisy_results],
                "status": "complete",
            },
        )
        return SampleGenerationResult(
            config.run_id,
            sample_id,
            seed,
            str(path_file),
            str(panel_file if not config.paths_only else ""),
            str(completion_file),
            "ok",
            time.perf_counter() - started,
            tuple(noisy_results),
        )
    except Exception as exception:
        completion_file.unlink(missing_ok=True)
        return SampleGenerationResult(
            config.run_id,
            sample_id,
            seed,
            str(path_file),
            str(panel_file),
            str(completion_file),
            "error",
            time.perf_counter() - started,
            tuple(noisy_results),
            f"{type(exception).__name__}: {exception}\n{traceback.format_exc()}",
        )


def _generate_one_star(arguments: tuple[int, ExperimentConfig]) -> SampleGenerationResult:
    return generate_one_sample(*arguments)


def run_samples(
    *, config: ExperimentConfig, sample_start: int = 0, sample_end: int | None = None
) -> list[SampleGenerationResult]:
    if config.paths_only and config.panels_only:
        raise ValueError("paths_only and panels_only cannot both be true")
    validate_panel_dependencies(config.panel_format)
    set_thread_env()
    end = config.n_samples if sample_end is None else min(sample_end, config.n_samples)
    if sample_start < 0 or end < sample_start:
        raise ValueError("sample range must satisfy 0 <= sample_start <= sample_end")
    sample_ids = list(range(sample_start, end))
    if not sample_ids:
        return []
    workers = default_workers(len(sample_ids), config.workers)
    LOGGER.info("generating samples [%d, %d) with %d worker(s)", sample_start, end, workers)
    if workers == 1:
        return [generate_one_sample(sample_id, config) for sample_id in sample_ids]
    context = mp.get_context("spawn" if "spawn" in mp.get_all_start_methods() else mp.get_all_start_methods()[0])
    with context.Pool(processes=workers) as pool:
        return list(pool.imap_unordered(_generate_one_star, [(sample_id, config) for sample_id in sample_ids]))


def write_run_metadata(config: ExperimentConfig, results: list[SampleGenerationResult]) -> None:
    config_directory = Path(config.output_root) / "config"
    requested_workers = config.workers
    resolved_workers = default_workers(max(len(results), 1), requested_workers)
    sample_range = (
        [min(result.sample_id for result in results), max(result.sample_id for result in results) + 1]
        if results
        else []
    )
    _atomic_json(config_directory / "experiment_config.json", config.raw_config)
    _atomic_json(
        config_directory / "run_metadata.json",
        {
            "git_sha": _git_sha(),
            "experiment_config_hash": config.experiment_config_hash,
            "requested_workers": requested_workers,
            "resolved_workers": resolved_workers,
            "sample_range": sample_range,
            "thread_environment": {key: os.environ.get(key) for key in THREAD_ENV_KEYS},
            "panel_format": config.panel_format,
            "python_version": platform.python_version(),
        },
    )
    manifest_path = config_directory / "manifest_generation.csv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = manifest_path.with_name(manifest_path.stem + ".tmp.csv")
    fieldnames = [
        "run_id", "sample_id", "seed", "path_file", "panel_file", "completion_file",
        "status", "elapsed_seconds", "error",
    ]
    with temporary.open("w", newline="") as file_handle:
        writer = csv.DictWriter(file_handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in sorted(results, key=lambda item: item.sample_id):
            payload = asdict(result)
            payload.pop("noisy_results")
            writer.writerow(payload)
        file_handle.flush()
        os.fsync(file_handle.fileno())
    os.replace(temporary, manifest_path)
    noisy_results = [noisy for result in results for noisy in result.noisy_results]
    if config.noise is not None:
        noisy_manifest = config_directory / "manifest_noisy_panels.csv"
        if noisy_results or not noisy_manifest.exists():
            write_noisy_manifest(config.output_root, noisy_results)


def write_noisy_metadata(config: ExperimentConfig, sample_results: list[SampleGenerationResult]) -> None:
    """Write the in-memory noisy manifest without rescanning panels."""

    if config.noise is not None:
        noisy_results = [noisy for result in sample_results for noisy in result.noisy_results]
        noisy_manifest = Path(config.output_root) / "config" / "manifest_noisy_panels.csv"
        if noisy_results or not noisy_manifest.exists():
            write_noisy_manifest(config.output_root, noisy_results)
