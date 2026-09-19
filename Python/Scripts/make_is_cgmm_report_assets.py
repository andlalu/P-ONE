from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from Scripts.make_noise_visualisations import (
    LINE_STYLES,
    MATURITY_LABELS,
    PdfFigure,
    _draw_axes,
    _draw_colorbar,
    _draw_heatmap,
    _draw_legend_horizontal,
    _fmt_decimal,
    _maturity_label,
)


SCENARIOS = ("clean", "low_iid", "spatial_corr", "persistent_factor")
NOISY_SCENARIOS = SCENARIOS[1:]
SCENARIO_LABELS = {
    "clean": "Clean",
    "low_iid": "Design A: low i.i.d.",
    "spatial_corr": "Design B: spatial correlation",
    "persistent_factor": "Design C: persistent factor",
}
SHORT_SCENARIO_LABELS = {
    "clean": "Clean",
    "low_iid": "Design A",
    "spatial_corr": "Design B",
    "persistent_factor": "Design C",
}
PARAMETERS = ("eta", "kappa", "vbar", "sigma_v", "rho", "kappa_q")
PARAMETER_LABELS = {
    "eta": r"$\eta$",
    "kappa": r"$\kappa$",
    "vbar": r"$\bar V$",
    "sigma_v": r"$\sigma_v$",
    "rho": r"$\rho$",
    "kappa_q": r"$\kappa_{\mathbb Q}$",
}
TRUE_PARAMETERS = {
    "eta": 5.0,
    "kappa": 7.0,
    "vbar": 0.0225,
    "sigma_v": 0.4,
    "rho": -0.5,
    "kappa_q": 2.0,
}
NATURAL_BOUNDS = {
    "eta": (2.0, 8.0),
    "kappa": (4.0, 9.0),
    "vbar": (0.012, 0.04),
    "sigma_v": (0.2, 0.65),
    "rho": (-0.8, -0.2),
    "kappa_q": (1.5, 4.0),
}
BOUND_RELATIVE_TOLERANCE = 1.0e-4
EXPECTED_CONFIG_HASH = "b74031ebf443c40496c2b9b3b4088df11ce49d761cc154d78050268dd015873f"
EXPECTED_GIT_SHA = "a9e17ee3880eddcd5cee31dfaf4854bdcf15d2e1"
EXPECTED_RUN_ID = "run_001"
REPRESENTATIVE_CONTRACT = (0.0, 0.25)
HORIZON_YEARS = 1.0 / 12.0


@dataclass(frozen=True)
class ValidatedReplication:
    sample_id: int
    record: dict[str, Any]
    path_file: Path
    panel_file: Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_finite_array(value: Any, *, shape: tuple[int, ...], label: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{label} must be finite with shape {shape}")
    return array


def _validate_recorded_artifact(sample_dir: Path, artifact: dict[str, Any], label: str) -> Path:
    path = sample_dir / str(artifact.get("file", ""))
    if not path.is_file():
        raise ValueError(f"missing {label} artifact for {sample_dir.name}")
    if path.stat().st_size != int(artifact.get("bytes", -1)):
        raise ValueError(f"{label} byte count differs from record for {sample_dir.name}")
    if sha256_file(path) != artifact.get("sha256"):
        raise ValueError(f"{label} checksum differs from record for {sample_dir.name}")
    return path


def _validate_estimate(result: dict[str, Any], *, sample_id: int, scenario: str) -> None:
    required_parameters = {"eta", "kappa", "vbar", "sigma_v", "rho", "eta_v", "r", "q"}
    parameters = result.get("estimated_parameters")
    if not isinstance(parameters, dict) or set(parameters) != required_parameters:
        raise ValueError(f"sample {sample_id:03d} {scenario} has an incomplete parameter vector")
    if not all(math.isfinite(float(value)) for value in parameters.values()):
        raise ValueError(f"sample {sample_id:03d} {scenario} has non-finite parameters")
    if not math.isfinite(float(result.get("final_criterion", math.nan))):
        raise ValueError(f"sample {sample_id:03d} {scenario} has a non-finite criterion")
    diagnostics = result.get("final_diagnostics")
    if not isinstance(diagnostics, dict):
        raise ValueError(f"sample {sample_id:03d} {scenario} lacks final diagnostics")
    if (
        int(diagnostics.get("n_dates", -1)) != 526
        or int(diagnostics.get("n_transitions", -1)) != 524
        or int(diagnostics.get("n_contracts", -1)) != 7890
    ):
        raise ValueError(f"sample {sample_id:03d} {scenario} has unexpected result dimensions")
    implied = diagnostics.get("implied_state")
    if not isinstance(implied, dict):
        raise ValueError(f"sample {sample_id:03d} {scenario} lacks an implied-state result")
    _require_finite_array(
        implied.get("variance"),
        shape=(526,),
        label=f"sample {sample_id:03d} {scenario} implied variance",
    )
    _require_finite_array(
        implied.get("objective"),
        shape=(526,),
        label=f"sample {sample_id:03d} {scenario} implied-state objective",
    )


def _validate_panel(panel_file: Path, *, sample_id: int) -> None:
    columns = [
        "sample_id",
        "scenario",
        "week_index",
        "model_iv",
        "estimation_iv",
        "observed_iv",
        "raw_noisy_iv",
        "cap_direction",
        "log_moneyness",
        "maturity_years",
    ]
    frame = pd.read_parquet(panel_file, columns=columns)
    if len(frame) != 31_560 or set(frame["scenario"].astype(str)) != set(SCENARIOS):
        raise ValueError(f"sample {sample_id:03d} has an incomplete combined panel")
    if set(frame["sample_id"].astype(int)) != {sample_id}:
        raise ValueError(f"sample {sample_id:03d} panel has a mismatched sample identifier")
    finite_columns = ["model_iv", "estimation_iv", "log_moneyness", "maturity_years"]
    if not np.isfinite(frame[finite_columns].to_numpy(dtype=float)).all():
        raise ValueError(f"sample {sample_id:03d} panel has non-finite required fields")
    clean_observed = frame.loc[frame["scenario"] == "clean", "observed_iv"]
    noisy_observed = frame.loc[frame["scenario"] != "clean", "observed_iv"]
    if not clean_observed.isna().all() or not np.isfinite(noisy_observed.to_numpy(dtype=float)).all():
        raise ValueError(f"sample {sample_id:03d} panel has unexpected observed-IV availability")
    clean_raw = frame.loc[frame["scenario"] == "clean", "raw_noisy_iv"]
    noisy_raw = frame.loc[frame["scenario"] != "clean", "raw_noisy_iv"]
    noisy_directions = set(frame.loc[frame["scenario"] != "clean", "cap_direction"].astype(str))
    if (
        not clean_raw.isna().all()
        or not np.isfinite(noisy_raw.to_numpy(dtype=float)).all()
        or not noisy_directions.issubset({"none", "lower", "upper"})
    ):
        raise ValueError(f"sample {sample_id:03d} panel has incomplete price-mechanics fields")
    sizes = frame.groupby(["scenario", "week_index"], sort=False).size()
    for scenario in SCENARIOS:
        local = sizes.loc[scenario]
        if len(local) != 526 or not local.eq(15).all():
            raise ValueError(f"sample {sample_id:03d} {scenario} does not contain 526 panels of 15 contracts")


def validate_snapshot(source_root: Path) -> tuple[list[ValidatedReplication], dict[str, Any]]:
    root = source_root.expanduser().resolve()
    run_file = root / "run.json"
    production_file = root / "PRODUCTION_SUCCESS"
    aggregate_file = root / "all400_summary.json"
    for path in (run_file, production_file, aggregate_file):
        if not path.is_file():
            raise ValueError(f"missing required snapshot file: {path.name}")
    run = json.loads(run_file.read_text())
    production = json.loads(production_file.read_text())
    if (
        run.get("run_id") != EXPECTED_RUN_ID
        or run.get("git_sha") != EXPECTED_GIT_SHA
        or run.get("configuration_hash") != EXPECTED_CONFIG_HASH
    ):
        raise ValueError("snapshot run identifiers do not match the production experiment")
    if production.get("validated_samples") != 100 or production.get("status") != "success":
        raise ValueError("production marker does not certify 100 validated samples")
    if production.get("scenario_results") != {scenario: 100 for scenario in SCENARIOS}:
        raise ValueError("production marker does not certify 100 results in every scenario")

    record_files = sorted(root.glob("sample_*/record.json"))
    sample_ids = [int(path.parent.name.removeprefix("sample_")) for path in record_files]
    if sample_ids != list(range(100)):
        raise ValueError(f"snapshot sample identifiers are {sample_ids}, expected exactly 000--099")
    replications: list[ValidatedReplication] = []
    for sample_id, record_file in enumerate(record_files):
        record = json.loads(record_file.read_text())
        if (
            record.get("sample_id") != sample_id
            or record.get("run_id") != EXPECTED_RUN_ID
            or record.get("git_sha") != EXPECTED_GIT_SHA
            or record.get("configuration_hash") != EXPECTED_CONFIG_HASH
            or record.get("status") != "complete"
            or record.get("current_stage") != "complete"
            or record.get("errors") != []
        ):
            raise ValueError(f"sample {sample_id:03d} record is incomplete or belongs to another run")
        artifacts = record.get("artifacts", {})
        path_file = _validate_recorded_artifact(record_file.parent, artifacts.get("path", {}), "path")
        panel_file = _validate_recorded_artifact(record_file.parent, artifacts.get("panels", {}), "panel")
        with np.load(path_file) as path:
            _require_finite_array(path["V_week"], shape=(526,), label=f"sample {sample_id:03d} true variance")
            _require_finite_array(path["dlogS_week"], shape=(525,), label=f"sample {sample_id:03d} weekly returns")
        _validate_panel(panel_file, sample_id=sample_id)
        estimates = record.get("estimation")
        if not isinstance(estimates, dict) or set(estimates) != set(SCENARIOS):
            raise ValueError(f"sample {sample_id:03d} does not contain all four scenarios")
        for scenario in SCENARIOS:
            _validate_estimate(estimates[scenario], sample_id=sample_id, scenario=scenario)
        replications.append(ValidatedReplication(sample_id, record, path_file, panel_file))

    identifiers = {
        "production_run_id": EXPECTED_RUN_ID,
        "git_sha": EXPECTED_GIT_SHA,
        "configuration_hash": EXPECTED_CONFIG_HASH,
        "production_marker_sha256": sha256_file(production_file),
        "run_manifest_sha256": sha256_file(run_file),
        "aggregate_summary_sha256": sha256_file(aggregate_file),
        "completed_at_utc": production.get("completed_at_utc"),
    }
    return replications, identifiers


def average_variance(v: np.ndarray | float, *, kappa: float, vbar: float, h: float) -> np.ndarray:
    if kappa <= 0.0 or h < 0.0:
        raise ValueError("kappa must be positive and h must be non-negative")
    values = np.asarray(v, dtype=float)
    if h == 0.0:
        return values.copy()
    loading = -math.expm1(-kappa * h) / (kappa * h)
    return vbar + (values - vbar) * loading


def risk_premium_paths(
    v: np.ndarray,
    *,
    eta: float,
    kappa: float,
    vbar: float,
    kappa_q: float,
    h: float = HORIZON_YEARS,
) -> tuple[np.ndarray, np.ndarray]:
    if kappa_q <= 0.0:
        raise ValueError("kappa_q must be positive")
    a_p = average_variance(v, kappa=kappa, vbar=vbar, h=h)
    vbar_q = kappa * vbar / kappa_q
    a_q = average_variance(v, kappa=kappa_q, vbar=vbar_q, h=h)
    return eta * a_p, a_p - a_q


def _correlation(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def paired_path_metrics(estimate: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    estimated = np.asarray(estimate, dtype=float)
    true = np.asarray(truth, dtype=float)
    if estimated.shape != true.shape or estimated.ndim != 1 or not np.isfinite(estimated).all() or not np.isfinite(true).all():
        raise ValueError("paired paths must be same-length finite vectors")
    error = estimated - true
    return {
        "path_average_truth": float(np.mean(true)),
        "path_average_estimate": float(np.mean(estimated)),
        "mean_error": float(np.mean(error)),
        "mae": float(np.mean(np.abs(error))),
        "rmse": float(np.sqrt(np.mean(error * error))),
        "correlation": _correlation(estimated, true),
    }


def summarise(values: Iterable[float]) -> dict[str, float]:
    array = np.asarray(list(values), dtype=float)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return {key: math.nan for key in ("n", "mean", "median", "sd", "q10", "q90")}
    return {
        "n": int(finite.size),
        "mean": float(np.mean(finite)),
        "median": float(np.median(finite)),
        "sd": float(np.std(finite, ddof=1)) if finite.size > 1 else 0.0,
        "q10": float(np.quantile(finite, 0.10)),
        "q90": float(np.quantile(finite, 0.90)),
    }


def _parameter_vector(result: dict[str, Any]) -> dict[str, float]:
    raw = result["estimated_parameters"]
    return {
        "eta": float(raw["eta"]),
        "kappa": float(raw["kappa"]),
        "vbar": float(raw["vbar"]),
        "sigma_v": float(raw["sigma_v"]),
        "rho": float(raw["rho"]),
        "kappa_q": float(raw["kappa"]) - float(raw["eta_v"]),
    }


def build_structural_outputs(
    replications: list[ValidatedReplication],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    estimate_rows: list[dict[str, Any]] = []
    for replication in replications:
        for scenario in SCENARIOS:
            estimate_rows.append(
                {
                    "sample_id": replication.sample_id,
                    "scenario": scenario,
                    **_parameter_vector(replication.record["estimation"][scenario]),
                }
            )
    estimates = pd.DataFrame(estimate_rows)

    summary_rows: list[dict[str, Any]] = []
    contact_rows: list[dict[str, Any]] = []
    for scenario in SCENARIOS:
        local = estimates.loc[estimates["scenario"] == scenario]
        for parameter in PARAMETERS:
            values = local[parameter].to_numpy(dtype=float)
            truth = TRUE_PARAMETERS[parameter]
            errors = values - truth
            rmse = float(np.sqrt(np.mean(errors * errors)))
            squared_errors = errors * errors
            bias_mcse = float(np.std(errors, ddof=1) / math.sqrt(values.size))
            rmse_mcse = (
                float(np.std(squared_errors, ddof=1) / math.sqrt(values.size) / (2.0 * rmse))
                if rmse > 0.0
                else 0.0
            )
            statistics = {
                "true_value": truth,
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "bias": float(np.mean(errors)),
                "empirical_sd": float(np.std(values, ddof=1)),
                "rmse": rmse,
                "q10": float(np.quantile(values, 0.10)),
                "q90": float(np.quantile(values, 0.90)),
            }
            for statistic, value in statistics.items():
                summary_rows.append(
                    {
                        "scenario": scenario,
                        "statistic": statistic,
                        "parameter": parameter,
                        "value": value,
                        "n": int(values.size),
                        "mcse": bias_mcse if statistic == "bias" else rmse_mcse if statistic == "rmse" else math.nan,
                    }
                )
            lower, upper = NATURAL_BOUNDS[parameter]
            # Match the production audit: measure the relative contact band against
            # the width of each pre-specified natural parameter interval.
            contact_tolerance = BOUND_RELATIVE_TOLERANCE * (upper - lower)
            contact_rows.append(
                {
                    "scenario": scenario,
                    "parameter": parameter,
                    "lower_bound": lower,
                    "upper_bound": upper,
                    "lower_contact_percent": 100.0 * float(np.mean(values <= lower + contact_tolerance)),
                    "upper_contact_percent": 100.0 * float(np.mean(values >= upper - contact_tolerance)),
                    "relative_tolerance": BOUND_RELATIVE_TOLERANCE,
                    "n": int(values.size),
                }
            )
    return estimates, pd.DataFrame(summary_rows), pd.DataFrame(contact_rows)


def build_state_and_risk_outputs(
    replications: list[ValidatedReplication],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    state_rows: list[dict[str, Any]] = []
    risk_rows: list[dict[str, Any]] = []
    true_parameter_vector = TRUE_PARAMETERS
    for replication in replications:
        with np.load(replication.path_file) as path:
            true_variance = np.asarray(path["V_week"], dtype=float)
        true_erp, true_vrp = risk_premium_paths(
            true_variance,
            eta=true_parameter_vector["eta"],
            kappa=true_parameter_vector["kappa"],
            vbar=true_parameter_vector["vbar"],
            kappa_q=true_parameter_vector["kappa_q"],
        )
        for scenario in SCENARIOS:
            result = replication.record["estimation"][scenario]
            implied_variance = np.asarray(
                result["final_diagnostics"]["implied_state"]["variance"], dtype=float
            )
            state_error = implied_variance - true_variance
            state_rows.append(
                {
                    "sample_id": replication.sample_id,
                    "scenario": scenario,
                    "mean_error": float(np.mean(state_error)),
                    "mae": float(np.mean(np.abs(state_error))),
                    "rmse": float(np.sqrt(np.mean(state_error * state_error))),
                    "median_absolute_error": float(np.median(np.abs(state_error))),
                    "q95_absolute_error": float(np.quantile(np.abs(state_error), 0.95)),
                    "correlation": _correlation(implied_variance, true_variance),
                }
            )
            parameters = _parameter_vector(result)
            estimated_erp, estimated_vrp = risk_premium_paths(
                implied_variance,
                eta=parameters["eta"],
                kappa=parameters["kappa"],
                vbar=parameters["vbar"],
                kappa_q=parameters["kappa_q"],
            )
            for quantity, estimate, truth in (
                ("erp", estimated_erp, true_erp),
                ("vrp", estimated_vrp, true_vrp),
            ):
                risk_rows.append(
                    {
                        "sample_id": replication.sample_id,
                        "scenario": scenario,
                        "quantity": quantity,
                        **paired_path_metrics(estimate, truth),
                    }
                )

    state_per_replication = pd.DataFrame(state_rows)
    state_summary_rows: list[dict[str, Any]] = []
    state_metrics = (
        "mean_error",
        "mae",
        "rmse",
        "median_absolute_error",
        "q95_absolute_error",
        "correlation",
    )
    for scenario in SCENARIOS:
        local = state_per_replication.loc[state_per_replication["scenario"] == scenario]
        for metric in state_metrics:
            state_summary_rows.append({"scenario": scenario, "metric": metric, **summarise(local[metric])})

    risk_per_replication = pd.DataFrame(risk_rows)
    risk_summary_rows: list[dict[str, Any]] = []
    risk_metrics = (
        "path_average_truth",
        "path_average_estimate",
        "mean_error",
        "mae",
        "rmse",
        "correlation",
    )
    for scenario in SCENARIOS:
        for quantity in ("erp", "vrp"):
            local = risk_per_replication.loc[
                (risk_per_replication["scenario"] == scenario)
                & (risk_per_replication["quantity"] == quantity)
            ]
            for metric in risk_metrics:
                risk_summary_rows.append(
                    {
                        "scenario": scenario,
                        "quantity": quantity,
                        "metric": metric,
                        **summarise(local[metric]),
                    }
                )
    return (
        state_per_replication,
        pd.DataFrame(state_summary_rows),
        risk_per_replication,
        pd.DataFrame(risk_summary_rows),
    )


def _lag_one_correlation(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return _correlation(array[:-1], array[1:])


def _build_noise_mechanics_frame(
    raw_cell_errors: dict[tuple[str, float, float], list[np.ndarray]],
    lower_cell_adjustments: dict[tuple[str, float, float], list[np.ndarray]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    contract_cells = sorted({(maturity, moneyness) for _, maturity, moneyness in raw_cell_errors})
    for scenario in (*NOISY_SCENARIOS, "pooled_noisy"):
        scenarios = NOISY_SCENARIOS if scenario == "pooled_noisy" else (scenario,)
        for maturity, moneyness in contract_cells:
            raw_values = np.concatenate(
                [
                    array
                    for local_scenario in scenarios
                    for array in raw_cell_errors[(local_scenario, maturity, moneyness)]
                ]
            ) * 10_000.0
            lower_adjusted = np.concatenate(
                [
                    array
                    for local_scenario in scenarios
                    for array in lower_cell_adjustments[(local_scenario, maturity, moneyness)]
                ]
            )
            rows.append(
                {
                    "scenario": scenario,
                    "maturity_years": maturity,
                    "log_moneyness": moneyness,
                    "q95_absolute_raw_error_bp": float(np.quantile(np.abs(raw_values), 0.95)),
                    "lower_bound_adjustment_percent": 100.0 * float(np.mean(lower_adjusted)),
                    "observations": int(raw_values.size),
                }
            )
    return pd.DataFrame(rows)


def build_noise_outputs(
    replications: list[ValidatedReplication],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    errors: dict[str, list[np.ndarray]] = {scenario: [] for scenario in SCENARIOS}
    cell_errors: dict[tuple[str, float, float], list[np.ndarray]] = {}
    cell_levels: dict[tuple[str, float, float], list[np.ndarray]] = {}
    raw_cell_errors: dict[tuple[str, float, float], list[np.ndarray]] = {}
    lower_cell_adjustments: dict[tuple[str, float, float], list[np.ndarray]] = {}
    dependence_rows: list[dict[str, Any]] = []
    for replication in replications:
        columns = [
            "scenario",
            "week_index",
            "maturity_years",
            "log_moneyness",
            "model_iv",
            "estimation_iv",
            "observed_iv",
            "raw_noisy_iv",
            "cap_direction",
        ]
        frame = pd.read_parquet(replication.panel_file, columns=columns)
        clean = frame.loc[frame["scenario"] == "clean"].copy()
        clean_error = np.zeros(len(clean), dtype=float)
        errors["clean"].append(clean_error)
        for (tau, moneyness), group in clean.groupby(["maturity_years", "log_moneyness"], sort=True):
            cell_levels.setdefault(("clean", float(tau), float(moneyness)), []).append(
                group["model_iv"].to_numpy(dtype=float)
            )
        for scenario in NOISY_SCENARIOS:
            local = frame.loc[frame["scenario"] == scenario].copy()
            local["signed_error"] = local["observed_iv"].astype(float) - local["model_iv"].astype(float)
            errors[scenario].append(local["signed_error"].to_numpy(dtype=float))
            for (tau, moneyness), group in local.groupby(["maturity_years", "log_moneyness"], sort=True):
                key = (scenario, float(tau), float(moneyness))
                cell_levels.setdefault(key, []).append(group["observed_iv"].to_numpy(dtype=float))
                cell_errors.setdefault(key, []).append(group["signed_error"].to_numpy(dtype=float))
                raw_cell_errors.setdefault(key, []).append(
                    group["raw_noisy_iv"].to_numpy(dtype=float) - group["model_iv"].to_numpy(dtype=float)
                )
                lower_cell_adjustments.setdefault(key, []).append(
                    group["cap_direction"].astype(str).eq("lower").to_numpy(dtype=bool)
                )

            pivot = local.pivot_table(
                index="week_index",
                columns=["maturity_years", "log_moneyness"],
                values="signed_error",
            ).sort_index(axis=1)
            correlation = pivot.corr().to_numpy(dtype=float)
            upper = correlation[np.triu_indices_from(correlation, k=1)]
            moneyness, maturity = REPRESENTATIVE_CONTRACT
            representative = local.loc[
                np.isclose(local["log_moneyness"].astype(float), moneyness)
                & np.isclose(local["maturity_years"].astype(float), maturity)
            ].sort_values("week_index")
            dependence_rows.append(
                {
                    "sample_id": replication.sample_id,
                    "scenario": scenario,
                    "mean_off_diagonal_cross_contract_correlation": float(np.mean(upper)),
                    "representative_contract_moneyness": moneyness,
                    "representative_contract_maturity_years": maturity,
                    "representative_contract_lag1_correlation": _lag_one_correlation(
                        representative["signed_error"].to_numpy(dtype=float)
                    ),
                }
            )

    noise_summary_rows: list[dict[str, Any]] = []
    for scenario in SCENARIOS:
        values = np.concatenate(errors[scenario]) * 10_000.0
        absolute = np.abs(values)
        local_dependence = [row for row in dependence_rows if row["scenario"] == scenario]
        noise_summary_rows.append(
            {
                "scenario": scenario,
                "mean_error_bp": float(np.mean(values)),
                "mae_bp": float(np.mean(absolute)),
                "rmse_bp": float(np.sqrt(np.mean(values * values))),
                "median_absolute_error_bp": float(np.median(absolute)),
                "q95_absolute_error_bp": float(np.quantile(absolute, 0.95)),
                "maximum_absolute_error_bp": float(np.max(absolute)),
                "observations": int(values.size),
                "mean_cross_contract_correlation": (
                    float(
                        np.mean(
                            [
                                row["mean_off_diagonal_cross_contract_correlation"]
                                for row in local_dependence
                            ]
                        )
                    )
                    if local_dependence
                    else math.nan
                ),
                "median_representative_contract_lag1_correlation": (
                    float(
                        np.median(
                            [
                                row["representative_contract_lag1_correlation"]
                                for row in local_dependence
                            ]
                        )
                    )
                    if local_dependence
                    else math.nan
                ),
            }
        )

    level_rows: list[dict[str, Any]] = []
    for key, arrays in sorted(cell_levels.items()):
        scenario, maturity, moneyness = key
        values = np.concatenate(arrays)
        level_rows.append(
            {
                "scenario": scenario,
                "maturity_years": maturity,
                "log_moneyness": moneyness,
                "mean_iv": float(np.mean(values)),
                "observations": int(values.size),
            }
        )
    error_surface_rows: list[dict[str, Any]] = []
    for key, arrays in sorted(cell_errors.items()):
        scenario, maturity, moneyness = key
        values = np.concatenate(arrays) * 10_000.0
        error_surface_rows.append(
            {
                "scenario": scenario,
                "maturity_years": maturity,
                "log_moneyness": moneyness,
                "mean_signed_error_bp": float(np.mean(values)),
                "q95_absolute_error_bp": float(np.quantile(np.abs(values), 0.95)),
                "observations": int(values.size),
            }
        )
    mechanics = _build_noise_mechanics_frame(raw_cell_errors, lower_cell_adjustments)
    return (
        pd.DataFrame(noise_summary_rows),
        pd.DataFrame(level_rows),
        pd.DataFrame(error_surface_rows),
        pd.DataFrame(dependence_rows),
        mechanics,
    )


def _surface_matrix(frame: pd.DataFrame, value_column: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    maturities = np.array(sorted(frame["maturity_years"].unique()), dtype=float)
    moneyness = np.array(sorted(frame["log_moneyness"].unique()), dtype=float)
    matrix = np.empty((maturities.size, moneyness.size), dtype=float)
    for row, maturity in enumerate(maturities):
        for column, log_moneyness in enumerate(moneyness):
            value = frame.loc[
                np.isclose(frame["maturity_years"], maturity)
                & np.isclose(frame["log_moneyness"], log_moneyness),
                value_column,
            ]
            if len(value) != 1:
                raise ValueError("surface data do not contain one value per contract cell")
            matrix[row, column] = float(value.iloc[0])
    return maturities, moneyness, matrix


def make_noise_surface_figure(levels: pd.DataFrame, output: Path) -> None:
    all_values = levels["mean_iv"].to_numpy(dtype=float)
    y_min = float(np.min(all_values)) - 0.002
    y_max = float(np.max(all_values)) + 0.002
    pdf = PdfFigure(7.2, 5.15)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, "Mean implied-volatility surfaces across 100 samples", size=10.0, align="center", bold=True)
    left, right, bottom, top, gap_x, gap_y = 64.0, 18.0, 62.0, 58.0, 32.0, 40.0
    panel_w = (pdf.width - left - right - gap_x) / 2.0
    panel_h = (pdf.height - bottom - top - gap_y) / 2.0
    rects = [
        (left, bottom + panel_h + gap_y, panel_w, panel_h),
        (left + panel_w + gap_x, bottom + panel_h + gap_y, panel_w, panel_h),
        (left, bottom, panel_w, panel_h),
        (left + panel_w + gap_x, bottom, panel_w, panel_h),
    ]
    for index, (scenario, rect) in enumerate(zip(SCENARIOS, rects)):
        local = levels.loc[levels["scenario"] == scenario]
        maturities, moneyness, matrix = _surface_matrix(local, "mean_iv")
        sx, sy = _draw_axes(
            pdf,
            rect,
            xmin=float(moneyness.min()),
            xmax=float(moneyness.max()),
            ymin=y_min,
            ymax=y_max,
            title=SHORT_SCENARIO_LABELS[scenario],
            xlabel="log-forward moneyness" if index >= 2 else "",
            ylabel="mean implied vol." if index % 2 == 0 else "",
            x_ticks=[float(value) for value in moneyness],
        )
        for maturity_index, maturity in enumerate(maturities):
            gray, dash = LINE_STYLES[maturity_index]
            points = [
                (sx(float(log_moneyness)), sy(float(matrix[maturity_index, column])))
                for column, log_moneyness in enumerate(moneyness)
            ]
            pdf.polyline(points, gray=gray, width=1.0, dash=dash)
    _draw_legend_horizontal(
        pdf,
        143.0,
        pdf.height - 38.0,
        [_maturity_label(value) for value in sorted(levels["maturity_years"].unique())],
    )
    pdf.save(output)


def make_noise_error_surface_figure(surfaces: pd.DataFrame, output: Path) -> None:
    mean_limit = float(np.max(np.abs(surfaces["mean_signed_error_bp"])))
    q95_max = float(np.max(surfaces["q95_absolute_error_bp"]))
    pdf = PdfFigure(7.2, 5.15)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, "Implied-volatility error surfaces across 100 samples", size=10.0, align="center", bold=True)
    left, right, bottom, top, gap_x, gap_y = 64.0, 18.0, 80.0, 58.0, 20.0, 48.0
    panel_w = (pdf.width - left - right - 2.0 * gap_x) / 3.0
    panel_h = (pdf.height - bottom - top - gap_y) / 2.0
    for column, scenario in enumerate(NOISY_SCENARIOS):
        local = surfaces.loc[surfaces["scenario"] == scenario]
        maturities, moneyness, mean_matrix = _surface_matrix(local, "mean_signed_error_bp")
        _, _, q95_matrix = _surface_matrix(local, "q95_absolute_error_bp")
        x_labels = [_fmt_decimal(value) for value in moneyness]
        y_labels = [_maturity_label(value).replace(" months", "m").replace(" month", "m") for value in maturities]
        x = left + column * (panel_w + gap_x)
        _draw_heatmap(
            pdf,
            (x, bottom + panel_h + gap_y, panel_w, panel_h),
            mean_matrix,
            x_labels=x_labels,
            y_labels=y_labels,
            title=SHORT_SCENARIO_LABELS[scenario],
            gray_for_value=lambda value, limit=max(mean_limit, 1.0e-12): 0.90 - 0.78 * ((value + limit) / (2.0 * limit)),
            value_label=lambda value: f"{value:.2f}",
            ylabel="mean error" if column == 0 else "",
        )
        _draw_heatmap(
            pdf,
            (x, bottom, panel_w, panel_h),
            q95_matrix,
            x_labels=x_labels,
            y_labels=y_labels,
            title="",
            gray_for_value=lambda value, maximum=max(q95_max, 1.0e-12): 0.92 - 0.76 * (value / maximum),
            value_label=lambda value: f"{value:.1f}",
            xlabel="log-forward moneyness" if column == 1 else "",
            ylabel="95% quantile" if column == 0 else "",
        )
    _draw_colorbar(
        pdf,
        80.0,
        35.0,
        135.0,
        8.0,
        min_label=f"-{mean_limit:.2f}",
        mid_label="0",
        max_label=f"{mean_limit:.2f}",
        title="mean signed error, IV bp",
        gray_fn=lambda value: 0.90 - 0.78 * value,
    )
    _draw_colorbar(
        pdf,
        315.0,
        35.0,
        135.0,
        8.0,
        min_label="0",
        mid_label=None,
        max_label=f"{q95_max:.1f}",
        title="95% quantile of absolute error, IV bp",
        gray_fn=lambda value: 0.92 - 0.76 * value,
    )
    pdf.save(output)


def make_noise_mechanics_figure(diagnostics: pd.DataFrame, output: Path) -> None:
    pooled = diagnostics.loc[diagnostics["scenario"] == "pooled_noisy"]
    maturities, moneyness, raw_matrix = _surface_matrix(pooled, "q95_absolute_raw_error_bp")
    _, _, lower_matrix = _surface_matrix(pooled, "lower_bound_adjustment_percent")
    x_labels = [_fmt_decimal(value) for value in moneyness]
    y_labels = [_maturity_label(value) for value in maturities]
    raw_max = float(np.max(raw_matrix))
    lower_max = float(np.max(lower_matrix))

    pdf = PdfFigure(7.2, 3.35)
    pdf.text(
        pdf.width / 2.0,
        pdf.height - 18.0,
        "Noise size and price-rounding effects across 100 samples",
        size=10.0,
        align="center",
        bold=True,
    )
    left, right, bottom, top, gap = 66.0, 18.0, 54.0, 58.0, 30.0
    panel_w = (pdf.width - left - right - gap) / 2.0
    panel_h = pdf.height - bottom - top
    specifications = (
        (
            raw_matrix,
            "95% quantile of absolute raw IV error",
            lambda value: 0.94 - 0.72 * value / max(raw_max, 1.0e-12),
            lambda value: f"{value:.1f}",
        ),
        (
            lower_matrix,
            "Lower-bound adjustment share (%)",
            lambda value: 0.94 - 0.72 * value / max(lower_max, 1.0e-12),
            lambda value: f"{value:.1f}",
        ),
    )
    for index, (matrix, title, gray_fn, value_fn) in enumerate(specifications):
        x = left + index * (panel_w + gap)
        _draw_heatmap(
            pdf,
            (x, bottom, panel_w, panel_h),
            matrix,
            x_labels=x_labels,
            y_labels=y_labels if index == 0 else [""] * len(y_labels),
            title=title,
            gray_for_value=gray_fn,
            value_label=value_fn,
            xlabel="log-forward moneyness",
        )
    pdf.text(left - 5.0, bottom + panel_h + 10.0, "Maturity", size=6.5, align="right", bold=True)
    pdf.save(output)


def _box_statistics(values: Iterable[float]) -> tuple[float, float, float, float, float]:
    array = np.asarray(list(values), dtype=float)
    return tuple(float(np.quantile(array, probability)) for probability in (0.10, 0.25, 0.50, 0.75, 0.90))


def _draw_boxplot_panel(
    pdf: PdfFigure,
    rect: tuple[float, float, float, float],
    groups: list[tuple[str, np.ndarray]],
    *,
    title: str,
    ylabel: str,
    zero_line: bool,
) -> None:
    x0, y0, width, height = rect
    statistics = [_box_statistics(values) for _, values in groups]
    low = min(item[0] for item in statistics)
    high = max(item[-1] for item in statistics)
    padding = max(0.05 * (high - low), 1.0e-6)
    low -= padding
    high += padding
    if zero_line:
        low, high = min(low, 0.0), max(high, 0.0)
    _, sy = _draw_axes(
        pdf,
        rect,
        xmin=0.5,
        xmax=len(groups) + 0.5,
        ymin=low,
        ymax=high,
        title=title,
        ylabel="",
        x_ticks=[],
    )
    if ylabel:
        pdf.text(x0 - 10.0, y0 + height / 2.0, ylabel, size=7.2, align="right")
    if zero_line:
        pdf.line(x0, sy(0.0), x0 + width, sy(0.0), gray=0.55, width=0.4, dash=(2.0, 2.0))
    for index, ((label, _), (q10, q25, median, q75, q90)) in enumerate(zip(groups, statistics), start=1):
        x = x0 + (index - 0.5) / len(groups) * width
        box_width = min(30.0, 0.55 * width / len(groups))
        pdf.line(x, sy(q10), x, sy(q90), gray=0.0, width=0.7)
        pdf.line(x - box_width / 3.0, sy(q10), x + box_width / 3.0, sy(q10), gray=0.0, width=0.7)
        pdf.line(x - box_width / 3.0, sy(q90), x + box_width / 3.0, sy(q90), gray=0.0, width=0.7)
        pdf.rect(x - box_width / 2.0, sy(q25), box_width, sy(q75) - sy(q25), fill_gray=0.84, stroke_gray=0.0, line_width=0.6)
        pdf.line(x - box_width / 2.0, sy(median), x + box_width / 2.0, sy(median), gray=0.0, width=1.1)
        pdf.text(x, y0 - 14.0, label, size=6.8, align="center")


def make_state_rmse_figure(state_per_replication: pd.DataFrame, output: Path) -> None:
    groups = [
        (SHORT_SCENARIO_LABELS[scenario], state_per_replication.loc[state_per_replication["scenario"] == scenario, "rmse"].to_numpy(dtype=float) * 100.0)
        for scenario in SCENARIOS
    ]
    pdf = PdfFigure(7.2, 3.4)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, "Implied-variance RMSE across 100 samples", size=10.0, align="center", bold=True)
    _draw_boxplot_panel(pdf, (82.0, 62.0, 404.0, 135.0), groups, title="", ylabel="variance x 100", zero_line=False)
    pdf.text(pdf.width - 18.0, 16.0, "Boxes show the interquartile range; whiskers show the 10th and 90th percentiles.", size=6.8, align="right")
    pdf.save(output)


def make_risk_error_figure(risk_per_replication: pd.DataFrame, output: Path) -> None:
    pdf = PdfFigure(7.2, 3.8)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, "Average one-month risk-premium errors across 100 samples", size=10.0, align="center", bold=True)
    for index, (quantity, title, ylabel) in enumerate(
        (("erp", "Equity risk premium", "annualised percentage points"), ("vrp", "Variance risk premium", "annualised variance x 100"))
    ):
        groups = [
            (
                SHORT_SCENARIO_LABELS[scenario],
                risk_per_replication.loc[
                    (risk_per_replication["scenario"] == scenario)
                    & (risk_per_replication["quantity"] == quantity),
                    "mean_error",
                ].to_numpy(dtype=float)
                * 100.0,
            )
            for scenario in SCENARIOS
        ]
        _draw_boxplot_panel(
            pdf,
            (82.0 + index * 218.0, 68.0, 190.0, 150.0),
            groups,
            title=title,
            ylabel="percentage points" if index == 0 else "",
            zero_line=True,
        )
    pdf.text(pdf.width - 18.0, 16.0, "Boxes show the interquartile range; whiskers show the 10th and 90th percentiles.", size=6.8, align="right")
    pdf.save(output)


def make_noise_dependence_figure(dependence: pd.DataFrame, output: Path) -> None:
    pdf = PdfFigure(7.2, 3.8)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, "Observation-error dependence across 100 samples", size=10.0, align="center", bold=True)
    specifications = (
        ("mean_off_diagonal_cross_contract_correlation", "Pairwise correlation across contracts"),
        ("representative_contract_lag1_correlation", "First-order autocorrelation"),
    )
    for index, (metric, title) in enumerate(specifications):
        groups = [
            (
                SHORT_SCENARIO_LABELS[scenario],
                dependence.loc[dependence["scenario"] == scenario, metric].to_numpy(dtype=float),
            )
            for scenario in NOISY_SCENARIOS
        ]
        _draw_boxplot_panel(
            pdf,
            (58.0 + index * 245.0, 68.0, 215.0, 150.0),
            groups,
            title=title,
            ylabel="correlation" if index == 0 else "",
            zero_line=True,
        )
    pdf.text(pdf.width - 18.0, 16.0, "Temporal measure: m=0 and maturity=3 months. Whiskers show the 10th and 90th percentiles.", size=6.8, align="right")
    pdf.save(output)


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, float_format="%.12g", lineterminator="\n")


def _parameter_format(parameter: str, value: float) -> str:
    decimals = 4 if parameter == "vbar" else 3
    return f"{value:.{decimals}f}"


def write_structural_table(summary: pd.DataFrame, path: Path) -> None:
    statistic_labels = (
        ("true_value", "True value"),
        ("mean", "Mean"),
        ("median", "Median"),
        ("bias", "Bias"),
        ("empirical_sd", "Empirical SD"),
        ("rmse", "RMSE"),
        ("q10", r"10\% quantile"),
        ("q90", r"90\% quantile"),
    )
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        "Statistic & " + " & ".join(PARAMETER_LABELS[parameter] for parameter in PARAMETERS) + r" \\",
        r"\midrule",
    ]
    panel_letters = ("A", "B", "C", "D")
    for panel_letter, scenario in zip(panel_letters, SCENARIOS):
        lines.append(
            rf"\multicolumn{{7}}{{l}}{{\textit{{Panel {panel_letter}: {SCENARIO_LABELS[scenario]} ($n=100$)}}}} \\"
        )
        local = summary.loc[summary["scenario"] == scenario]
        for statistic, label in statistic_labels:
            values = []
            for parameter in PARAMETERS:
                value = float(
                    local.loc[
                        (local["statistic"] == statistic) & (local["parameter"] == parameter),
                        "value",
                    ].iloc[0]
                )
                values.append(_parameter_format(parameter, value))
            lines.append(label + " & " + " & ".join(values) + r" \\")
        if scenario != SCENARIOS[-1]:
            lines.append(r"\addlinespace")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\caption{Distribution of the first-step IS-CGMM parameter estimates across the 100 samples. The risk-neutral parameter is $\kappa_{\mathbb Q}=\kappa-\eta_V$. Bias, empirical standard deviation and RMSE are calculated separately for each design.}",
            r"\label{tab:mc_structural_parameters}",
            r"\end{table}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_noise_table(summary: pd.DataFrame, path: Path) -> None:
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\small",
        r"\setlength{\tabcolsep}{5pt}",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r"Design & Mean & MAE & RMSE & Median abs. & 95\% quantile & Maximum abs. \\",
        r"\midrule",
    ]
    for scenario in SCENARIOS:
        row = summary.loc[summary["scenario"] == scenario].iloc[0]
        values = [
            float(row[column])
            for column in (
                "mean_error_bp",
                "mae_bp",
                "rmse_bp",
                "median_absolute_error_bp",
                "q95_absolute_error_bp",
                "maximum_absolute_error_bp",
            )
        ]
        lines.append(SHORT_SCENARIO_LABELS[scenario] + " & " + " & ".join(f"{value:.2f}" for value in values) + r" \\")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\caption{Implied-volatility observation errors across the 100 samples, measured in implied-volatility basis points. Signed error is observed minus model implied volatility, and the 95\% quantile is calculated from the absolute errors. For each design, the statistics combine the errors from all 100 samples, 526 dates and 15 contracts, giving 789,000 observations. The clean benchmark is zero by construction.}",
            r"\label{tab:mc_noise_summary}",
            r"\end{table}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def _summary_mean(summary: pd.DataFrame, *, scenario: str, metric: str, quantity: str | None = None) -> float:
    mask = (summary["scenario"] == scenario) & (summary["metric"] == metric)
    if quantity is not None:
        mask &= summary["quantity"] == quantity
    return float(summary.loc[mask, "mean"].iloc[0])


def write_state_table(summary: pd.DataFrame, path: Path) -> None:
    metrics = (
        ("mean_error", "Mean error"),
        ("mae", "MAE"),
        ("rmse", "RMSE"),
        ("median_absolute_error", "Median abs."),
        ("q95_absolute_error", r"95\% quantile"),
    )
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{lrrrrr}",
        r"\toprule",
        "Design & " + " & ".join(label for _, label in metrics) + r" \\",
        r"\midrule",
    ]
    for scenario in SCENARIOS:
        values = []
        for metric, _ in metrics:
            value = _summary_mean(summary, scenario=scenario, metric=metric)
            values.append(f"{100.0 * value:.3f}")
        lines.append(SHORT_SCENARIO_LABELS[scenario] + " & " + " & ".join(values) + r" \\")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\caption{Accuracy of the final implied variance sequence. The error measures are first calculated within each sample using 526 matched dates and are then averaged across the 100 samples. The 95\% quantile is calculated from the absolute errors. Variance errors are multiplied by 100.}",
            r"\label{tab:mc_state_accuracy}",
            r"\end{table}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_risk_table(summary: pd.DataFrame, path: Path) -> None:
    metrics = (
        ("path_average_truth", "Truth"),
        ("path_average_estimate", "Estimate"),
        ("mean_error", "Mean error"),
        ("mae", "MAE"),
        ("rmse", "RMSE"),
    )
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\small",
        r"\setlength{\tabcolsep}{5pt}",
        r"\begin{tabular}{lrrrrr}",
        r"\toprule",
        "Design & " + " & ".join(label for _, label in metrics) + r" \\",
        r"\midrule",
    ]
    panels = (
        ("A", "erp", "ERP: annualised percentage points"),
        ("B", "vrp", "VRP: annualised variance multiplied by 100"),
    )
    for panel, quantity, label in panels:
        lines.append(rf"\multicolumn{{6}}{{l}}{{\textit{{Panel {panel}: {label}}}}} \\")
        for scenario in SCENARIOS:
            values = []
            for metric, _ in metrics:
                value = _summary_mean(summary, scenario=scenario, metric=metric, quantity=quantity)
                values.append(f"{100.0 * value:.4f}")
            lines.append(SHORT_SCENARIO_LABELS[scenario] + " & " + " & ".join(values) + r" \\")
        if quantity == "erp":
            lines.append(r"\addlinespace")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\caption{Accuracy of the one-month equity and variance risk premia. Each measure is first calculated from the matched estimated and true paths within a sample and is then averaged across the 100 samples. The variance risk premium is defined as physical minus risk-neutral expected average variance.}",
            r"\label{tab:mc_risk_premium_accuracy}",
            r"\end{table}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_bound_contact_table(contacts: pd.DataFrame, path: Path) -> None:
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{llrrrrrr}",
        r"\toprule",
        "Design & Bound & " + " & ".join(PARAMETER_LABELS[parameter] for parameter in PARAMETERS) + r" \\",
        r"\midrule",
    ]
    for scenario in SCENARIOS:
        local = contacts.loc[contacts["scenario"] == scenario]
        for side, column in (("Lower", "lower_contact_percent"), ("Upper", "upper_contact_percent")):
            values = [
                float(local.loc[local["parameter"] == parameter, column].iloc[0])
                for parameter in PARAMETERS
            ]
            lines.append(
                (SHORT_SCENARIO_LABELS[scenario] if side == "Lower" else "")
                + " & "
                + side
                + " & "
                + " & ".join(f"{value:.1f}" for value in values)
                + r" \\"
            )
        if scenario != SCENARIOS[-1]:
            lines.append(r"\addlinespace")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\caption{Parameter-bound contact percentages included as an additional numerical diagnostic. An estimate is counted as being close to a bound when its distance from that bound is no more than $10^{-4}$ times the width of the corresponding parameter interval.}",
            r"\label{tab:mc_bound_contacts}",
            r"\end{table}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def _parse_generation_timestamp(value: str | None) -> str:
    if value is None:
        return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise ValueError("generation timestamp must include a timezone")
    return parsed.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def generate_assets(
    *,
    source_root: Path,
    output_dir: Path,
    source_label: str,
    generation_timestamp: str | None,
) -> dict[str, Any]:
    replications, identifiers = validate_snapshot(source_root)
    output = output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)

    estimates, structural_summary, contacts = build_structural_outputs(replications)
    state_per_replication, state_summary, risk_per_replication, risk_summary = build_state_and_risk_outputs(replications)
    noise_summary, noise_levels, noise_surfaces, noise_dependence, noise_mechanics = build_noise_outputs(replications)

    csv_outputs = {
        "structural_parameter_estimates.csv": estimates,
        "structural_parameter_summary.csv": structural_summary,
        "structural_parameter_mcse.csv": structural_summary.loc[
            structural_summary["statistic"].isin(["bias", "rmse"]),
            ["scenario", "statistic", "parameter", "value", "mcse", "n"],
        ],
        "bound_contacts.csv": contacts,
        "noise_summary.csv": noise_summary,
        "noise_surface_means.csv": noise_levels,
        "noise_error_surfaces.csv": noise_surfaces,
        "noise_mechanics_diagnostics.csv": noise_mechanics,
        "noise_dependence_per_replication.csv": noise_dependence,
        "implied_state_per_replication.csv": state_per_replication,
        "implied_state_summary.csv": state_summary,
        "risk_premium_per_replication.csv": risk_per_replication,
        "risk_premium_summary.csv": risk_summary,
    }
    for filename, frame in csv_outputs.items():
        _write_csv(frame, output / filename)

    table_outputs = {
        "structural_parameter_table.tex": (write_structural_table, structural_summary),
        "bound_contacts_table.tex": (write_bound_contact_table, contacts),
        "noise_summary_table.tex": (write_noise_table, noise_summary),
        "implied_state_table.tex": (write_state_table, state_summary),
        "risk_premium_table.tex": (write_risk_table, risk_summary),
    }
    for filename, (writer, frame) in table_outputs.items():
        writer(frame, output / filename)

    figure_outputs = {
        "noise_mechanics_diagnostics.pdf": (make_noise_mechanics_figure, noise_mechanics),
        "noise_dependence.pdf": (make_noise_dependence_figure, noise_dependence),
        "implied_state_rmse_boxplot.pdf": (make_state_rmse_figure, state_per_replication),
        "risk_premium_error_boxplots.pdf": (make_risk_error_figure, risk_per_replication),
    }
    for filename, (writer, frame) in figure_outputs.items():
        writer(frame, output / filename)

    generated_files = sorted(
        [*csv_outputs, *table_outputs, *figure_outputs]
    )
    manifest = {
        **identifiers,
        "source_result_location": source_label,
        "generation_timestamp_utc": _parse_generation_timestamp(generation_timestamp),
        "included_sample_ids_by_scenario": {
            scenario: [f"{replication.sample_id:03d}" for replication in replications]
            for scenario in SCENARIOS
        },
        "included_replication_count_by_scenario": {scenario: len(replications) for scenario in SCENARIOS},
        "inclusion_rule": (
            "Include a replication when its sample metadata match the production run, the recorded path and panel files "
            "match their byte counts and SHA-256 checksums, all four scenario records and required dimensions are "
            "present, and the required parameter, criterion, true-state and implied-state fields are finite. Internal "
            "search flags are not used for inclusion."
        ),
        "state_matching_rule": "Match the 526 final implied-variance values to the 526 simulated weekly variance values by stored week index within sample.",
        "risk_premium_rule": (
            "Apply the one-month nonlinear transformations separately to every replication and date; define VRP as physical minus risk-neutral average variance."
        ),
        "noise_rule": "Signed implied-volatility error is observed_iv minus model_iv for noisy designs; the clean benchmark uses model_iv minus model_iv.",
        "noise_mechanics_rule": (
            "Within each maturity and moneyness cell, raw IV error is raw_noisy_iv minus model_iv before price rounding. "
            "The lower-bound adjustment share is the percentage of observations with cap_direction equal to lower after "
            "rounding. The plotted cells pool the three matched-marginal noisy designs."
        ),
        "representative_contract": {
            "log_forward_moneyness": REPRESENTATIVE_CONTRACT[0],
            "maturity_years": REPRESENTATIVE_CONTRACT[1],
        },
        "generated_assets": {
            filename: {"sha256": sha256_file(output / filename), "bytes": (output / filename).stat().st_size}
            for filename in generated_files
        },
    }
    manifest_path = output / "reporting_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate deterministic IS-CGMM Monte Carlo report assets from a validated snapshot.")
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--source-label", default="private production snapshot run_001-a9e17ee")
    parser.add_argument("--generation-timestamp", default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    manifest = generate_assets(
        source_root=Path(arguments.source_root),
        output_dir=Path(arguments.output_dir),
        source_label=str(arguments.source_label),
        generation_timestamp=arguments.generation_timestamp,
    )
    print(f"validated samples: {manifest['included_replication_count_by_scenario']['clean']}")
    print(f"configuration hash: {manifest['configuration_hash']}")
    print(f"generated assets: {len(manifest['generated_assets'])}")
    print(f"manifest: {Path(arguments.output_dir) / 'reporting_manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
