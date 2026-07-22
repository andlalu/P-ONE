from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from DGPSimulation.types import HestonSimConfig
from Estimation.ISCGMM.config import (
    CcfQuadratureConfig,
    CgmmConfig,
    ImpliedStateConfig,
    OptimizerConfig,
    PowellStageConfig,
)
from Models.Heston.parameters import HestonParameters, HestonPhysicalParameters
from OptionPricing.cos_basis import FixedCosBasisConfig
from OptionPricing.noisy_panel import (
    NoiseSettings,
    compute_q_diag_from_stationary_std,
    compute_stationary_std_from_q_diag,
    validate_noise_settings,
)


@dataclass(frozen=True)
class ExperimentConfig:
    """Resolved experiment inputs plus independent numerical settings."""

    source_path: str
    raw_config: dict[str, Any]
    experiment_config_hash: str
    run_id: str
    n_samples: int
    base_seed: int
    output_root: str
    panel_format: str
    workers: int | None
    dgp: HestonPhysicalParameters
    simulation: HestonSimConfig
    eta_v: float
    log_moneyness: tuple[float, ...]
    atm_option_type: str
    cos_basis: FixedCosBasisConfig
    noise: NoiseSettings | None
    criterion_config: CgmmConfig
    optimizer_config: OptimizerConfig
    skip_existing: bool = False
    paths_only: bool = False
    panels_only: bool = False


def canonical_experiment_json(raw_config: dict[str, Any]) -> str:
    return json.dumps(
        raw_config,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def experiment_config_hash(raw_config: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_experiment_json(raw_config).encode("utf-8")).hexdigest()


def _resolve_relative(base_directory: Path, configured_path: str) -> Path:
    path = Path(configured_path).expanduser()
    return (base_directory / path).resolve() if not path.is_absolute() else path.resolve()


def _normalise_noise_settings(raw_noise: dict[str, Any] | None) -> NoiseSettings | None:
    if raw_noise is None:
        return None
    raw_scenarios = raw_noise.get("scenarios", {})
    if not isinstance(raw_scenarios, dict):
        raise ValueError("noise.scenarios must be a JSON object")
    scenarios: dict[str, dict[str, Any]] = {}
    for name, raw_scenario in raw_scenarios.items():
        if not isinstance(raw_scenario, dict):
            raise ValueError(f"noise scenario {name} must be a JSON object")
        scenario = {
            "alpha_0": float(raw_scenario["alpha_0"]),
            "alpha_m": float(raw_scenario["alpha_m"]),
            "alpha_tau": float(raw_scenario["alpha_tau"]),
            "tau_min": float(raw_scenario["tau_min"]),
        }
        if name == "spatial_corr":
            scenario.update(
                ell_m=float(raw_scenario["ell_m"]),
                ell_tau=float(raw_scenario["ell_tau"]),
                correlation_jitter=float(raw_scenario["correlation_jitter"]),
                max_correlation_jitter=float(raw_scenario["max_correlation_jitter"]),
            )
        elif name == "persistent_factor":
            a_diag = np.asarray(raw_scenario["a_diag"], dtype=float)
            stationary_std_raw = raw_scenario.get("stationary_factor_std")
            q_diag_raw = raw_scenario.get("q_diag")
            if stationary_std_raw is None and q_diag_raw is None:
                raise ValueError("persistent_factor requires stationary_factor_std or q_diag")
            stationary_std = (
                compute_stationary_std_from_q_diag(a_diag, q_diag_raw)
                if stationary_std_raw is None
                else np.asarray(stationary_std_raw, dtype=float)
            )
            q_diag = (
                compute_q_diag_from_stationary_std(a_diag, stationary_std)
                if q_diag_raw is None
                else np.asarray(q_diag_raw, dtype=float)
            )
            scenario.update(
                a_diag=a_diag.tolist(),
                stationary_factor_std=stationary_std.tolist(),
                q_diag=q_diag.tolist(),
                residual_policy=str(raw_scenario.get("residual_policy", "match_total_marginal_scale")),
                factor_initialization=str(raw_scenario.get("factor_initialization", "zero")),
            )
            if raw_scenario.get("residual_scale_multiplier") is not None:
                scenario["residual_scale_multiplier"] = float(raw_scenario["residual_scale_multiplier"])
        scenarios[str(name)] = scenario
    settings = NoiseSettings(
        base_seed=int(raw_noise["base_seed"]),
        sigma_min=float(raw_noise["sigma_min"]),
        price_epsilon=float(raw_noise["price_epsilon"]),
        tick_size=float(raw_noise["tick_size"]),
        scenarios=scenarios,
    )
    validate_noise_settings(settings)
    return settings


def _validate_experiment(config: ExperimentConfig) -> None:
    if not config.run_id or config.n_samples <= 0:
        raise ValueError("run.run_id must be non-empty and run.n_samples must be positive")
    if config.panel_format not in {"parquet", "csv"}:
        raise ValueError("run.panel_format must be 'parquet' or 'csv'")
    if config.workers is not None and config.workers <= 0:
        raise ValueError("run.workers must be positive or null")
    config.dgp.validate()
    config.simulation.validate()
    if not math.isfinite(config.eta_v):
        raise ValueError("q_measure.eta_v must be finite")
    if config.atm_option_type not in {"call", "put"}:
        raise ValueError("panel.atm_option_type must be 'call' or 'put'")
    if not config.log_moneyness or not all(math.isfinite(value) for value in config.log_moneyness):
        raise ValueError("panel.log_moneyness must be non-empty and finite")
    config.cos_basis.validate()
    config.cos_basis.validate_requested_maturities(config.cos_basis.maturities)
    config.criterion_config.implied_state.validate()
    config.criterion_config.quadrature.validate()
    config.criterion_config.validate()
    config.optimizer_config.base_start.validate()
    config.optimizer_config.stage1.validate()
    config.optimizer_config.stage2.validate()
    config.optimizer_config.validate()
    if config.noise is not None:
        validate_noise_settings(config.noise)


def load_experiment_config(config_path: str | Path) -> ExperimentConfig:
    source = Path(config_path).expanduser().resolve()
    with source.open() as file_handle:
        raw = json.load(file_handle)
    if not isinstance(raw, dict):
        raise ValueError("experiment configuration must contain one JSON object")
    for section in ("run", "dgp", "simulation", "q_measure", "panel", "cos", "estimation"):
        if section not in raw:
            raise ValueError(f"experiment configuration is missing required section {section!r}")

    run = raw["run"]
    panel = raw["panel"]
    cos = raw["cos"]
    estimation = raw["estimation"]
    state = estimation["implied_state"]
    quadrature = estimation["quadrature"]
    cgmm = estimation["cgmm"]
    optimizer = estimation["optimizer"]
    stages = optimizer["stages"]
    if len(stages) != 2:
        raise ValueError("estimation.optimizer.stages must contain exactly two Powell stages")

    cos_basis = FixedCosBasisConfig(
        maturities=tuple(float(value) for value in cos["maturities_years"]),
        effective_widths=tuple(float(value) for value in cos["effective_widths"]),
        generation_n_cos=int(cos["generation_n_cos"]),
        estimation_n_cos=int(cos["estimation_n_cos"]),
        maturity_tolerance=float(cos["maturity_tolerance"]),
    )
    implied_state = ImpliedStateConfig(
        cos_basis=cos_basis,
        v_min=float(state["v_min"]),
        v_max=float(state["v_max"]),
        tol=float(state["tol"]),
        max_iter=int(state["max_iter"]),
        boundary_tol=float(state["boundary_tol"]),
        warm_start_window=None if state.get("warm_start_window") is None else float(state["warm_start_window"]),
        state_solver=str(state["state_solver"]),
        fallback_solver=state.get("fallback_solver"),
        finite_difference_relative_step=float(state["finite_difference_relative_step"]),
        finite_difference_absolute_step=float(state["finite_difference_absolute_step"]),
        minimum_black_vega=float(state["minimum_black_vega"]),
    )
    criterion = CgmmConfig(
        implied_state=implied_state,
        quadrature=CcfQuadratureConfig(
            dimension=int(quadrature["dimension"]),
            order=int(quadrature["order"]),
            scale=float(quadrature["scale"]),
        ),
        instrument_precision=tuple(float(value) for value in cgmm["instrument_precision"]),
        transition_rk_steps=int(cgmm["transition_rk_steps"]),
        transition_cf_method=str(cgmm["transition_cf_method"]),
        dt=None if cgmm.get("dt") is None else float(cgmm["dt"]),
        spacing_tolerance=float(cgmm["spacing_tolerance"]),
    )
    powell = OptimizerConfig(
        base_start=HestonParameters(**optimizer["base_start"]),
        natural_bounds=tuple(tuple(float(value) for value in pair) for pair in optimizer["natural_bounds"]),
        candidate_relative_perturbations=tuple(
            tuple(float(value) for value in perturbation)
            for perturbation in optimizer["candidate_relative_perturbations"]
        ),
        stage1=PowellStageConfig(
            max_evaluations=int(stages[0]["max_evaluations"]),
            xtol=float(stages[0]["xtol"]),
            ftol=float(stages[0]["ftol"]),
        ),
        stage2=PowellStageConfig(
            max_evaluations=int(stages[1]["max_evaluations"]),
            xtol=float(stages[1]["xtol"]),
            ftol=float(stages[1]["ftol"]),
        ),
        penalty_value=float(optimizer["penalty_value"]),
        progress_every=int(optimizer["progress_every"]),
    )
    config = ExperimentConfig(
        source_path=str(source),
        raw_config=raw,
        experiment_config_hash=experiment_config_hash(raw),
        run_id=str(run["run_id"]),
        n_samples=int(run["n_samples"]),
        base_seed=int(run["base_seed"]),
        output_root=str(_resolve_relative(source.parent, str(run["output_root"]))),
        panel_format=str(run["panel_format"]),
        workers=None if run.get("workers") is None else int(run["workers"]),
        dgp=HestonPhysicalParameters(**raw["dgp"]),
        simulation=HestonSimConfig(**raw["simulation"]),
        eta_v=float(raw["q_measure"]["eta_v"]),
        log_moneyness=tuple(float(value) for value in panel["log_moneyness"]),
        atm_option_type=str(panel["atm_option_type"]),
        cos_basis=cos_basis,
        noise=_normalise_noise_settings(raw.get("noise")),
        criterion_config=criterion,
        optimizer_config=powell,
    )
    _validate_experiment(config)
    return config
