from __future__ import annotations

import csv
import json
import multiprocessing as mp
import os
import time
import traceback
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

from DGPSimulation.heston_simulator import HestonPathSimulator
from DGPSimulation.io import load_heston_path_npz, save_heston_path_npz
from DGPSimulation.types import HestonParamsP, HestonSimConfig
from DGPSimulation.variance_drawers import AndersenQeVarianceDrawer
from OptionPricing.clean_panel import generate_clean_option_panel_rows, parquet_available, write_panel
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.types import CosPricingConfig

THREAD_ENV_KEYS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")


@dataclass(frozen=True)
class PanelConfig:
    maturities_years: tuple[float, ...]
    log_moneyness: tuple[float, ...]
    option_type_rule: str
    atm_option_type: str
    pricing_method: str
    iv_method: str
    cos_n_terms: int
    cos_truncation_L: float


@dataclass(frozen=True)
class GenerationConfig:
    run_id: str
    n_samples: int
    base_seed: int
    output_root: str
    dgp: HestonParamsP
    simulation: HestonSimConfig
    eta_v: float
    panel: PanelConfig
    workers: int | None = None
    skip_existing: bool = False
    paths_only: bool = False
    panels_only: bool = False


@dataclass(frozen=True)
class SampleGenerationResult:
    run_id: str
    sample_id: int
    seed: int
    path_file: str
    panel_file: str
    status: str
    elapsed_seconds: float
    error: str = ""


def set_thread_env() -> None:
    for key in THREAD_ENV_KEYS:
        os.environ.setdefault(key, "1")


def load_config(config_path: str | Path) -> GenerationConfig:
    with Path(config_path).open() as fh:
        raw = json.load(fh)
    sim_raw = dict(raw["simulation"])
    sim_raw["return_daily"] = True
    panel_raw = raw["panel"]
    cos_raw = panel_raw.get("cos", {})
    return GenerationConfig(
        run_id=raw["run_id"],
        n_samples=int(raw["n_samples"]),
        base_seed=int(raw["base_seed"]),
        output_root=raw["output_root"],
        dgp=HestonParamsP(**raw["dgp"]),
        simulation=HestonSimConfig(**sim_raw),
        eta_v=float(raw["q_measure"]["eta_v"]),
        panel=PanelConfig(
            maturities_years=tuple(float(x) for x in panel_raw["maturities_years"]),
            log_moneyness=tuple(float(x) for x in panel_raw["log_moneyness"]),
            option_type_rule=panel_raw["option_type_rule"],
            atm_option_type=panel_raw.get("atm_option_type", "call"),
            pricing_method=panel_raw["pricing_method"],
            iv_method=panel_raw["iv_method"],
            cos_n_terms=int(cos_raw.get("n_terms", 256)),
            cos_truncation_L=float(cos_raw.get("truncation_L", 10.0)),
        ),
        workers=raw.get("parallel", {}).get("workers"),
    )


def resolved_config_dict(config: GenerationConfig) -> dict[str, Any]:
    payload = asdict(config)
    payload["panel_format"] = "parquet" if parquet_available() else "csv"
    payload["thread_env"] = {key: os.environ.get(key, "1") for key in THREAD_ENV_KEYS}
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
    if n_samples is not None:
        updates["n_samples"] = n_samples
    if workers is not None:
        updates["workers"] = workers
    if output_root is not None:
        updates["output_root"] = output_root
    if skip_existing is not None:
        updates["skip_existing"] = skip_existing
    if paths_only is not None:
        updates["paths_only"] = paths_only
    if panels_only is not None:
        updates["panels_only"] = panels_only
    return replace(config, **updates)


def sample_artifact_paths(config: GenerationConfig, sample_id: int) -> tuple[Path, Path]:
    root = Path(config.output_root)
    stem = f"sample_{sample_id:03d}"
    path_file = root / "paths" / f"{stem}.npz"
    panel_base = root / "panels_clean" / stem
    panel_file = panel_base.with_suffix(".parquet" if parquet_available() else ".csv")
    return path_file, panel_file


def generate_one_sample(sample_id: int, config: GenerationConfig) -> SampleGenerationResult:
    started = time.perf_counter()
    seed = config.base_seed + sample_id
    path_file, panel_file = sample_artifact_paths(config, sample_id)

    try:
        if config.skip_existing:
            path_ok = config.panels_only or path_file.exists()
            panel_ok = config.paths_only or panel_file.exists()
            if path_ok and panel_ok:
                return SampleGenerationResult(
                    run_id=config.run_id,
                    sample_id=sample_id,
                    seed=seed,
                    path_file=str(path_file),
                    panel_file=str(panel_file if panel_file.exists() else ""),
                    status="skipped",
                    elapsed_seconds=time.perf_counter() - started,
                )

        sim_config = replace(config.simulation, seed=seed, return_daily=True)

        if config.panels_only:
            path, params, _ = load_heston_path_npz(path_file)
        else:
            simulator = HestonPathSimulator(
                params=config.dgp,
                config=sim_config,
                variance_drawer=AndersenQeVarianceDrawer(),
            )
            path = simulator.simulate()
            params = config.dgp
            save_heston_path_npz(path_file, path=path, params=params, config=sim_config)

        if not config.paths_only:
            pricer = CosOptionPricer()
            pricing_config = CosPricingConfig(
                n_cos=config.panel.cos_n_terms,
                truncation_width=config.panel.cos_truncation_L,
            )
            rows = generate_clean_option_panel_rows(
                run_id=config.run_id,
                sample_id=sample_id,
                path=path,
                params_p=params,
                eta_v=config.eta_v,
                maturities_years=config.panel.maturities_years,
                log_moneyness=config.panel.log_moneyness,
                atm_option_type=config.panel.atm_option_type,
                pricing_method=config.panel.pricing_method,
                iv_method=config.panel.iv_method,
                pricer=pricer,
                pricing_config=pricing_config,
            )
            panel_file = write_panel(rows, Path(config.output_root) / "panels_clean" / f"sample_{sample_id:03d}")

        return SampleGenerationResult(
            run_id=config.run_id,
            sample_id=sample_id,
            seed=seed,
            path_file=str(path_file),
            panel_file=str(panel_file if not config.paths_only else ""),
            status="ok",
            elapsed_seconds=time.perf_counter() - started,
        )
    except Exception as exc:
        return SampleGenerationResult(
            run_id=config.run_id,
            sample_id=sample_id,
            seed=seed,
            path_file=str(path_file),
            panel_file=str(panel_file),
            status="error",
            elapsed_seconds=time.perf_counter() - started,
            error=f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
        )


def _generate_one_star(args: tuple[int, GenerationConfig]) -> SampleGenerationResult:
    return generate_one_sample(*args)


def default_workers(n_samples: int, requested: int | None) -> int:
    if requested is not None:
        return max(1, min(requested, n_samples))
    return min(max((os.cpu_count() or 2) - 1, 1), n_samples)


def run_samples(
    *,
    config: GenerationConfig,
    sample_start: int = 0,
    sample_end: int | None = None,
) -> list[SampleGenerationResult]:
    set_thread_env()
    end = config.n_samples if sample_end is None else min(sample_end, config.n_samples)
    sample_ids = list(range(sample_start, end))
    if not sample_ids:
        return []

    workers = default_workers(len(sample_ids), config.workers)
    if workers == 1:
        return [generate_one_sample(sample_id, config) for sample_id in sample_ids]

    start_method = "spawn" if "spawn" in mp.get_all_start_methods() else mp.get_all_start_methods()[0]
    ctx = mp.get_context(start_method)
    with ctx.Pool(processes=workers) as pool:
        return list(pool.imap_unordered(_generate_one_star, [(sample_id, config) for sample_id in sample_ids]))


def write_run_metadata(config: GenerationConfig, results: list[SampleGenerationResult]) -> None:
    config_dir = Path(config.output_root) / "config"
    config_dir.mkdir(parents=True, exist_ok=True)
    with (config_dir / "run_config_resolved.json").open("w") as fh:
        json.dump(resolved_config_dict(config), fh, indent=2)

    manifest_path = config_dir / "manifest_generation.csv"
    fieldnames = [
        "run_id",
        "sample_id",
        "seed",
        "path_file",
        "panel_file",
        "status",
        "elapsed_seconds",
        "error",
    ]
    with manifest_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for result in sorted(results, key=lambda item: item.sample_id):
            writer.writerow(asdict(result))
