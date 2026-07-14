from __future__ import annotations

import csv
import hashlib
import json
import logging
import multiprocessing as mp
import os
import platform
import subprocess
import sys
import time
import traceback
from dataclasses import asdict, dataclass, replace
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any

import numpy as np

from DGPSimulation.heston_simulator import HestonPathSimulator
from DGPSimulation.io import load_heston_path_npz, save_heston_path_npz
from DGPSimulation.types import HestonSimConfig
from DGPSimulation.variance_drawers import AndersenQeVarianceDrawer
from Models.Heston.parameters import HestonPhysicalParameters
from OptionData.io import panel_metadata_path
from OptionPricing.clean_panel import generate_clean_option_panel_rows, parquet_available, write_panel
from OptionPricing.cos_basis import FixedCosBasisConfig
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.noisy_panel import (
    NoiseGenerationConfig,
    NoisyPanelResult,
    generate_noisy_panel_file,
    parse_noise_config,
    write_noise_config,
    write_noisy_manifest,
)

LOGGER = logging.getLogger(__name__)
THREAD_ENV_KEYS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")
SUPPORTED_OPTION_RULE = "otm_by_forward_moneyness"
SUPPORTED_PRICING_METHOD = "COS"
SUPPORTED_IV_METHOD = "lets_be_rational"


@dataclass(frozen=True)
class PanelConfig:
    maturities_years: tuple[float, ...]
    log_moneyness: tuple[float, ...]
    option_type_rule: str
    atm_option_type: str
    pricing_method: str
    iv_method: str

    def validate(self) -> None:
        if self.option_type_rule != SUPPORTED_OPTION_RULE:
            raise ValueError(f"option_type_rule must be {SUPPORTED_OPTION_RULE!r}")
        if self.pricing_method != SUPPORTED_PRICING_METHOD:
            raise ValueError(f"pricing_method must be {SUPPORTED_PRICING_METHOD!r}")
        if self.iv_method != SUPPORTED_IV_METHOD:
            raise ValueError(f"iv_method must be {SUPPORTED_IV_METHOD!r}")
        if self.atm_option_type not in {"call", "put"}:
            raise ValueError("atm_option_type must be 'call' or 'put'")
        if not self.maturities_years or any(value <= 0.0 for value in self.maturities_years):
            raise ValueError("panel maturities must be non-empty and positive")
        if not self.log_moneyness or not all(math_is_finite(value) for value in self.log_moneyness):
            raise ValueError("panel log-moneyness grid must be non-empty and finite")


@dataclass(frozen=True)
class GenerationConfig:
    run_id: str
    n_samples: int
    base_seed: int
    output_root: str
    dgp: HestonPhysicalParameters
    simulation: HestonSimConfig
    eta_v: float
    panel: PanelConfig
    panel_format: str
    cos_basis: FixedCosBasisConfig
    cos_basis_path: str
    daily_array_semantics: str = (
        "When return_daily is true, saved daily arrays include the initial state and burn-in observations, "
        "followed by the production observations."
    )
    source_config_path: str = ""
    configuration_hash: str = ""
    workers: int | None = None
    skip_existing: bool = False
    paths_only: bool = False
    panels_only: bool = False
    noise: NoiseGenerationConfig | None = None

    def validate(self) -> None:
        if not self.run_id or self.n_samples <= 0:
            raise ValueError("run_id must be non-empty and n_samples must be positive")
        if self.panel_format not in {"parquet", "csv"}:
            raise ValueError("panel_format must be 'parquet' or 'csv'")
        if self.paths_only and self.panels_only:
            raise ValueError("paths_only and panels_only cannot both be true")
        self.dgp.validate()
        self.simulation.validate()
        self.panel.validate()
        self.cos_basis.validate()
        self.cos_basis.validate_requested_maturities(self.panel.maturities_years)


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


def math_is_finite(value: float) -> bool:
    return bool(np.isfinite(float(value)))


def _stable_json_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w") as file_handle:
        json.dump(payload, file_handle, indent=2, sort_keys=True)
        file_handle.write("\n")
        file_handle.flush()
        os.fsync(file_handle.fileno())
    os.replace(temporary, path)


def _resolve_relative(base_directory: Path, configured_path: str) -> Path:
    path = Path(configured_path).expanduser()
    return (base_directory / path).resolve() if not path.is_absolute() else path.resolve()


def load_config(config_path: str | Path) -> GenerationConfig:
    source = Path(config_path).expanduser().resolve()
    with source.open() as file_handle:
        raw = json.load(file_handle)
    if not isinstance(raw, dict):
        raise ValueError("generation configuration must contain one JSON object")
    panel_raw = raw["panel"]
    basis_path = _resolve_relative(source.parent, str(raw["cos_basis_path"]))
    output_root = _resolve_relative(source.parent, str(raw["output_root"]))
    noise_raw = raw.get("noise")
    config = GenerationConfig(
        run_id=str(raw["run_id"]),
        n_samples=int(raw["n_samples"]),
        base_seed=int(raw["base_seed"]),
        output_root=str(output_root),
        dgp=HestonPhysicalParameters(**raw["dgp"]),
        simulation=HestonSimConfig(**raw["simulation"]),
        eta_v=float(raw["q_measure"]["eta_v"]),
        panel=PanelConfig(
            maturities_years=tuple(float(value) for value in panel_raw["maturities_years"]),
            log_moneyness=tuple(float(value) for value in panel_raw["log_moneyness"]),
            option_type_rule=str(panel_raw["option_type_rule"]),
            atm_option_type=str(panel_raw.get("atm_option_type", "call")),
            pricing_method=str(panel_raw["pricing_method"]),
            iv_method=str(panel_raw["iv_method"]),
        ),
        panel_format=str(raw["panel_format"]),
        cos_basis=FixedCosBasisConfig.load(basis_path),
        cos_basis_path=str(basis_path),
        daily_array_semantics=str(
            raw.get(
                "daily_array_semantics",
                "When return_daily is true, saved daily arrays include the initial state and burn-in "
                "observations, followed by the production observations.",
            )
        ),
        source_config_path=str(source),
        configuration_hash=_stable_json_hash(raw),
        workers=raw.get("parallel", {}).get("workers"),
        noise=None if noise_raw is None else parse_noise_config(noise_raw),
    )
    config.validate()
    validate_panel_dependencies(config.panel_format)
    return config


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


def _package_version(name: str) -> str:
    try:
        return importlib_metadata.version(name)
    except importlib_metadata.PackageNotFoundError:
        return "unavailable"


def resolved_config_dict(
    config: GenerationConfig, *, requested_workers: int | None = None, resolved_workers: int | None = None
) -> dict[str, Any]:
    payload = asdict(config)
    payload["cos_basis_hash"] = config.cos_basis.stable_hash()
    payload["execution"] = {
        "git_sha": _git_sha(),
        "configuration_hash": config.configuration_hash,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": _package_version("scipy"),
        "pandas_version": _package_version("pandas"),
        "pyarrow_version": _package_version("pyarrow"),
        "requested_worker_count": requested_workers,
        "resolved_worker_count": resolved_workers,
        "thread_environment": {key: os.environ.get(key) for key in THREAD_ENV_KEYS},
        "platform": platform.platform(),
        "panel_format": config.panel_format,
    }
    return payload


def with_overrides(
    config: GenerationConfig,
    *,
    n_samples: int | None = None,
    workers: int | None = None,
    output_root: str | None = None,
    skip_existing: bool | None = None,
    paths_only: bool | None = None,
    panels_only: bool | None = None,
) -> GenerationConfig:
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
    result.validate()
    return result


def default_workers(n_samples: int, requested: int | None) -> int:
    if requested is not None:
        return max(1, min(int(requested), n_samples))
    return min(max((os.cpu_count() or 2) - 1, 1), n_samples)


def sample_artifact_paths(config: GenerationConfig, sample_id: int) -> tuple[Path, Path, Path]:
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


def sample_is_complete(config: GenerationConfig, sample_id: int) -> bool:
    path_file, panel_file, completion_file = sample_artifact_paths(config, sample_id)
    required: list[Path] = []
    if not config.panels_only:
        required.append(path_file)
    if not config.paths_only:
        required.extend([panel_file, panel_metadata_path(panel_file)])
        if config.noise is not None:
            for scenario in config.noise.enabled_scenarios():
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
        if completion.get("configuration_hash") != config.configuration_hash:
            return False
        if completion.get("cos_basis_hash") != config.cos_basis.stable_hash():
            return False
        if not config.paths_only:
            metadata = _read_json(panel_metadata_path(panel_file))
            config.cos_basis.assert_panel_compatible(metadata["cos_basis"])
            expected_rows = (
                (config.simulation.t_week + 1)
                * len(config.panel.maturities_years)
                * len(config.panel.log_moneyness)
            )
            if int(completion.get("clean_panel_rows", -1)) != expected_rows:
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


def generate_one_sample(sample_id: int, config: GenerationConfig) -> SampleGenerationResult:
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
                maturities_years=config.panel.maturities_years,
                log_moneyness=config.panel.log_moneyness,
                atm_option_type=config.panel.atm_option_type,
                pricing_method=config.panel.pricing_method,
                iv_method=config.panel.iv_method,
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
                    "configuration_hash": config.configuration_hash,
                    "git_sha": _git_sha(),
                    "panel_format": config.panel_format,
                    "cos_basis": config.cos_basis.generation_metadata(
                        resolved_basis_path=config.cos_basis_path
                    ),
                },
            )
            if config.noise is not None:
                for scenario in config.noise.enabled_scenarios():
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
                "run_id": config.run_id,
                "sample_id": sample_id,
                "seed": seed,
                "configuration_hash": config.configuration_hash,
                "cos_basis_hash": config.cos_basis.stable_hash(),
                "git_sha": _git_sha(),
                "panel_format": config.panel_format,
                "clean_panel_rows": clean_rows,
                "noise_scenarios": [] if config.noise is None else list(config.noise.enabled_scenarios()),
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


def _generate_one_star(arguments: tuple[int, GenerationConfig]) -> SampleGenerationResult:
    return generate_one_sample(*arguments)


def run_samples(
    *, config: GenerationConfig, sample_start: int = 0, sample_end: int | None = None
) -> list[SampleGenerationResult]:
    config.validate()
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


def write_run_metadata(config: GenerationConfig, results: list[SampleGenerationResult]) -> None:
    config_directory = Path(config.output_root) / "config"
    requested_workers = config.workers
    resolved_workers = default_workers(max(len(results), 1), requested_workers)
    resolved_payload = resolved_config_dict(
        config, requested_workers=requested_workers, resolved_workers=resolved_workers
    )
    if results:
        resolved_payload["execution"]["sample_range"] = [
            min(result.sample_id for result in results),
            max(result.sample_id for result in results) + 1,
        ]
    else:
        resolved_payload["execution"]["sample_range"] = []
    _atomic_json(config_directory / "run_config_resolved.json", resolved_payload)
    if config.noise is not None:
        write_noise_config(config.output_root, config.noise)
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


def write_noisy_metadata(config: GenerationConfig, sample_results: list[SampleGenerationResult]) -> None:
    """Write the in-memory noisy manifest without rescanning or regenerating panels."""

    if config.noise is not None:
        noisy_results = [noisy for result in sample_results for noisy in result.noisy_results]
        noisy_manifest = Path(config.output_root) / "config" / "manifest_noisy_panels.csv"
        if noisy_results or not noisy_manifest.exists():
            write_noisy_manifest(config.output_root, noisy_results)
