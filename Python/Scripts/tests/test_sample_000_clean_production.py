from __future__ import annotations

import json
import math
import os
import time
from pathlib import Path

for key in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(key, "1")

import numpy as np
import pandas as pd
import pytest

from DGPSimulation.io import load_heston_path_npz
from Scripts.experiment_config import load_experiment_config
from Scripts.sample_run import (
    PERSISTENT_FACTOR_COLUMNS,
    SCENARIO_ORDER,
    VARIANCE_LINKED_FACTOR_COLUMN,
    run_sample,
    sha256_file,
)

PRODUCTION_CONFIG = (
    Path(__file__).resolve().parents[1]
    / "configs"
    / "heston_experiment_run_002.json"
)


@pytest.mark.production_integration
@pytest.mark.slow
def test_sample_000_clean_production(tmp_path):
    config = load_experiment_config(PRODUCTION_CONFIG)
    output_root = tmp_path / "run_002"

    generation_started = time.perf_counter()
    generation_record = run_sample(
        config=config,
        sample_id=0,
        output_root=output_root,
        overwrite=True,
        generation_only=True,
        log_level="INFO",
        allow_dirty=True,
    )
    generation_elapsed = time.perf_counter() - generation_started
    print(f"production generation phase runtime: {generation_elapsed:.3f}s")

    sample_dir = output_root / "sample_000"
    path_file = sample_dir / "path.npz"
    panel_file = sample_dir / "panels.parquet"
    record_file = sample_dir / "record.json"
    log_file = sample_dir / "sample.log"
    assert {path.name for path in output_root.iterdir()} == {
        "run.json",
        "sample_000",
    }
    assert {path.name for path in sample_dir.iterdir()} == {
        "path.npz",
        "panels.parquet",
        "record.json",
        "sample.log",
    }
    for legacy in (
        "config",
        "logs",
        "paths",
        "panels_clean",
        "panels_observed",
        "noise_factors",
        "status",
        "validation",
    ):
        assert not (output_root / legacy).exists()
        assert not (sample_dir / legacy).exists()
    assert not list(output_root.rglob("*.metadata.json"))
    assert not list(output_root.rglob("*_SUCCESS"))
    assert not list(output_root.rglob("*_FAILED"))
    assert not [
        path
        for path in output_root.rglob("*")
        if path.is_file() and "factor" in path.name.lower()
    ]
    assert not [
        path
        for path in output_root.rglob("*")
        if ".tmp" in path.name
    ]

    run_payload = json.loads((output_root / "run.json").read_text())
    assert run_payload["run_id"] == "run_002"
    assert run_payload["configuration_hash"] == config.experiment_config_hash
    assert run_payload["configuration"] == config.raw_config
    assert run_payload["git_sha"]
    assert all(run_payload["environment"][name] for name in (
        "python",
        "numpy",
        "scipy",
        "pandas",
        "pyarrow",
    ))

    path, _, simulation = load_heston_path_npz(path_file)
    assert path.seed == 1234500
    assert simulation.seed == 1234500
    assert len(path.t_week) == 526
    assert len(path.logS_week) == 526
    assert len(path.V_week) == 526
    assert len(path.dlogS_week) == 525
    assert path.logS_daily is not None
    assert path.V_daily is not None
    for values in (
        path.t_week,
        path.logS_week,
        path.V_week,
        path.dlogS_week,
        path.logS_daily,
        path.V_daily,
    ):
        assert np.all(np.isfinite(values))
    assert np.all(path.V_week >= 0.0)
    assert np.all(path.V_daily >= 0.0)

    frame = pd.read_parquet(panel_file)
    assert tuple(frame["scenario"].drop_duplicates()) == SCENARIO_ORDER
    assert len(frame) == 39_450
    assert frame.groupby("scenario", sort=False).size().to_dict() == {
        scenario: 7_890 for scenario in SCENARIO_ORDER
    }
    assert frame.groupby("scenario", sort=False)["week_index"].nunique().to_dict() == {
        scenario: 526 for scenario in SCENARIO_ORDER
    }
    assert (
        frame.groupby(["scenario", "week_index"], sort=False).size() == 15
    ).all()

    scenario_rank = {name: index for index, name in enumerate(SCENARIO_ORDER)}
    order_keys = [
        (
            scenario_rank[row.scenario],
            int(row.week_index),
            float(row.maturity_years),
            float(row.log_moneyness),
        )
        for row in frame.itertuples(index=False)
    ]
    assert order_keys == sorted(order_keys)

    clean = frame[frame["scenario"] == "clean"]
    np.testing.assert_array_equal(clean["estimation_price"], clean["model_price"])
    np.testing.assert_array_equal(clean["estimation_iv"], clean["model_iv"])
    expected_seeds = {
        "path": 1_234_500,
        "low_iid": 9_900_101,
        "spatial_corr": 9_900_202,
        "persistent_factor": 9_900_303,
        "variance_linked_factor": 9_900_404,
    }
    assert generation_record["seeds"] == expected_seeds
    for scenario in SCENARIO_ORDER[1:]:
        noisy = frame[frame["scenario"] == scenario]
        np.testing.assert_array_equal(
            noisy["estimation_price"],
            noisy["observed_price"],
        )
        np.testing.assert_array_equal(
            noisy["estimation_iv"],
            noisy["observed_iv"],
        )
        assert set(noisy["noise_seed"].astype(int)) == {expected_seeds[scenario]}

    persistent = frame[frame["scenario"] == "persistent_factor"]
    assert persistent[list(PERSISTENT_FACTOR_COLUMNS)].notna().all().all()
    assert persistent[VARIANCE_LINKED_FACTOR_COLUMN].isna().all()
    assert persistent["week_index"].nunique() == 526
    assert (
        persistent.groupby("week_index")[list(PERSISTENT_FACTOR_COLUMNS)]
        .nunique()
        .to_numpy()
        == 1
    ).all()
    assert np.all(
        np.isfinite(
            persistent[list(PERSISTENT_FACTOR_COLUMNS)]
            .astype(float)
            .to_numpy()
        )
    )
    variance_linked = frame[frame["scenario"] == "variance_linked_factor"]
    assert variance_linked[list(PERSISTENT_FACTOR_COLUMNS)].notna().all().all()
    assert variance_linked[VARIANCE_LINKED_FACTOR_COLUMN].notna().all()
    assert (
        variance_linked.groupby("week_index")[
            [*PERSISTENT_FACTOR_COLUMNS, VARIANCE_LINKED_FACTOR_COLUMN]
        ]
        .nunique()
        .to_numpy()
        == 1
    ).all()
    for scenario in ("clean", "low_iid", "spatial_corr"):
        selected = frame[frame["scenario"] == scenario]
        assert selected[list(PERSISTENT_FACTOR_COLUMNS)].isna().all().all()
        assert selected[VARIANCE_LINKED_FACTOR_COLUMN].isna().all()

    assert generation_record["status"] == "generated"
    assert generation_record["current_stage"] == "validated"
    assert generation_record["completed_at_utc"] is None
    assert generation_record["validation"]["passed"] is True
    assert generation_record["artifacts"]["path"]["sha256"] == sha256_file(path_file)
    assert generation_record["artifacts"]["panels"]["sha256"] == sha256_file(panel_file)
    assert generation_record["estimation"] == {}

    path_hash = generation_record["artifacts"]["path"]["sha256"]
    panel_hash = generation_record["artifacts"]["panels"]["sha256"]
    path_mtime = path_file.stat().st_mtime_ns
    panel_mtime = panel_file.stat().st_mtime_ns

    estimation_started = time.perf_counter()
    clean_record = run_sample(
        config=config,
        sample_id=0,
        output_root=output_root,
        resume=True,
        scenario="clean",
        log_level="INFO",
        allow_dirty=True,
    )
    estimation_elapsed = time.perf_counter() - estimation_started
    clean_result = clean_record["estimation"]["clean"]
    print(f"production clean estimation runtime: {estimation_elapsed:.3f}s")
    print(f"production clean optimiser success: {clean_result['success']}")
    print(f"production clean final criterion: {clean_result['final_criterion']:.17g}")
    print(f"production clean estimated parameters: {clean_result['estimated_parameters']}")

    assert clean_record["artifacts"]["path"]["sha256"] == path_hash
    assert clean_record["artifacts"]["panels"]["sha256"] == panel_hash
    assert path_file.stat().st_mtime_ns == path_mtime
    assert panel_file.stat().st_mtime_ns == panel_mtime
    log_after_estimation = log_file.read_text()
    assert "path reused" in log_after_estimation
    assert "combined panels reused" in log_after_estimation
    assert tuple(clean_record["estimation"]) == ("clean",)
    assert clean_record["status"] == "estimating"
    assert clean_record["current_stage"] == "estimating"

    assert clean_result["candidate_starts"]
    assert clean_result["selected_start"]
    assert math.isfinite(clean_result["initial_criterion"])
    assert math.isfinite(clean_result["final_criterion"])
    assert clean_result["final_criterion"] <= clean_result["initial_criterion"]
    assert len(clean_result["free_parameters"]) == 6
    assert np.all(np.isfinite(clean_result["free_parameters"]))
    assert clean_result["function_evaluations"] >= 1
    assert clean_result["penalty_evaluations"] >= 0
    assert math.isfinite(clean_result["total_runtime_seconds"])
    assert clean_result["total_runtime_seconds"] > 0.0
    assert [item["pass_name"] for item in clean_result["powell_passes"]] == [
        "coarse",
        "refinement",
    ]

    parameters = clean_result["estimated_parameters"]
    assert all(math.isfinite(value) for value in parameters.values())
    natural = (
        parameters["eta"],
        parameters["kappa"],
        parameters["vbar"],
        parameters["sigma_v"],
        parameters["rho"],
        parameters["kappa"] - parameters["eta_v"],
    )
    for value, (lower, upper) in zip(natural, config.powell_config.natural_bounds):
        assert lower <= value <= upper

    diagnostics = clean_result["final_diagnostics"]
    assert diagnostics["n_dates"] == 526
    assert diagnostics["n_contracts"] == 7_890
    assert diagnostics["criterion_value"] == clean_result["final_criterion"]
    implied_state = diagnostics["implied_state"]
    failed = np.asarray(implied_state["failed"], dtype=bool)
    assert failed.shape == (526,)
    assert float(np.mean(~failed)) == 1.0
    assert implied_state["solver_name"] == "bounded_brent"
    assert math.isfinite(diagnostics["true_variance_rmse"])
    assert math.isfinite(diagnostics["true_variance_mae"])
    assert diagnostics["true_variance_rmse"] >= 0.0
    assert diagnostics["true_variance_mae"] >= 0.0
    assert not any(
        scenario in clean_record["estimation"]
        for scenario in SCENARIO_ORDER[1:]
    )

    evaluation_count = clean_result["function_evaluations"]
    final_criterion = clean_result["final_criterion"]
    estimated_parameters = dict(clean_result["estimated_parameters"])
    estimate_runtime = clean_result["total_runtime_seconds"]
    resumed_record = run_sample(
        config=config,
        sample_id=0,
        output_root=output_root,
        resume=True,
        scenario="clean",
        log_level="INFO",
        allow_dirty=True,
    )
    resumed_clean = resumed_record["estimation"]["clean"]
    assert path_file.stat().st_mtime_ns == path_mtime
    assert panel_file.stat().st_mtime_ns == panel_mtime
    assert resumed_clean["function_evaluations"] == evaluation_count
    assert resumed_clean["final_criterion"] == final_criterion
    assert resumed_clean["estimated_parameters"] == estimated_parameters
    assert resumed_clean["total_runtime_seconds"] == estimate_runtime
    assert "completed estimate reused" in log_file.read_text()
    assert tuple(resumed_record["estimation"]) == ("clean",)
    assert {
        path.name
        for path in output_root.iterdir()
        if path.is_dir() and path.name.startswith("sample_")
    } == {"sample_000"}
