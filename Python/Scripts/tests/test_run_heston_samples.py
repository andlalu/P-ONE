from __future__ import annotations

import json
from pathlib import Path

import pytest

from Scripts.run_heston_samples import RangeSummary, SampleOutcome, main, run_range
from Scripts.s3_persistence import S3PersistenceError

PRODUCTION_CONFIG = (
    Path(__file__).resolve().parents[1]
    / "configs"
    / "heston_experiment_run_001.json"
)


def _config(tmp_path: Path, *, n_samples: int = 8) -> Path:
    payload = json.loads(PRODUCTION_CONFIG.read_text())
    payload["run"].update(
        run_id="launcher_test",
        n_samples=n_samples,
        output_root="unused",
        workers=None,
    )
    path = tmp_path / "experiment.json"
    path.write_text(json.dumps(payload))
    return path


def _estimate(success: bool = True) -> dict:
    return {
        "success": success,
        "final_criterion": 0.25,
        "estimated_parameters": {
            "eta": 4.8,
            "kappa": 6.8,
            "vbar": 0.024,
            "sigma_v": 0.42,
            "rho": -0.52,
            "eta_v": 4.8,
            "r": 0.02,
            "q": 0.0,
        },
        "free_parameters": [4.8, 1.9, -3.7, -0.87, -0.57, 0.69],
        "function_evaluations": 12,
        "penalty_evaluations": 0,
        "powell_passes": [
            {"pass_name": "coarse"},
            {"pass_name": "refinement"},
        ],
        "final_diagnostics": {"criterion_value": 0.25},
        "total_runtime_seconds": 1.0,
    }


def _outcome(task, *, success=True):
    estimation = (
        {}
        if task.scenario is None
        else {task.scenario: _estimate(success=success)}
    )
    return SampleOutcome(
        sample_id=task.sample_id,
        record={"status": "estimating", "estimation": estimation},
        elapsed_seconds=0.1,
    )


def test_exclusive_range_and_serial_argument_routing(tmp_path, monkeypatch):
    config = _config(tmp_path)
    tasks = []

    def run(task):
        tasks.append(task)
        return _outcome(task)

    monkeypatch.setattr("Scripts.run_heston_samples._run_one_sample", run)
    summary = run_range(
        config_path=str(config),
        output_root=tmp_path / "run",
        sample_start=0,
        sample_end=1,
        sample_workers=1,
        scenario="clean",
        resume=True,
        max_dates=12,
        log_level="WARNING",
        allow_dirty=True,
    )

    assert [task.sample_id for task in tasks] == [0]
    task = tasks[0]
    assert task.resume is True
    assert task.overwrite is False
    assert task.scenario == "clean"
    assert task.max_dates == 12
    assert task.log_level == "WARNING"
    assert task.allow_dirty is True
    assert summary.estimate_counts() == (1, 0)
    assert summary.exit_code == 0


def test_range_zero_to_eight_submits_each_sample_once_with_bounded_pool(
    tmp_path,
    monkeypatch,
    capsys,
):
    config = _config(tmp_path)
    submitted = []
    pool_workers = []

    class Future:
        def __init__(self, task):
            self.task = task

        def result(self):
            return _outcome(self.task, success=self.task.sample_id != 7)

    class Executor:
        def __init__(self, max_workers):
            pool_workers.append(max_workers)

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def submit(self, function, task):
            assert function.__name__ == "_run_one_sample"
            submitted.append(task.sample_id)
            return Future(task)

    monkeypatch.setattr("Scripts.run_heston_samples.ProcessPoolExecutor", Executor)
    monkeypatch.setattr(
        "Scripts.run_heston_samples.as_completed",
        lambda futures: reversed(list(futures)),
    )
    summary = run_range(
        config_path=str(config),
        output_root=tmp_path / "run",
        sample_start=0,
        sample_end=8,
        sample_workers=3,
        scenario="clean",
        resume=True,
        allow_dirty=True,
    )

    assert pool_workers == [3]
    assert submitted == list(range(8))
    assert len(set(submitted)) == 8
    assert summary.estimate_counts() == (7, 1)
    completion_lines = [
        line for line in capsys.readouterr().out.splitlines()
        if line.startswith("sample=")
    ]
    assert completion_lines[0].startswith("sample=007")
    assert completion_lines[-1].startswith("sample=000")


def test_restore_requested_range_finishes_before_serial_workers(
    tmp_path,
    monkeypatch,
):
    config = _config(tmp_path, n_samples=6)
    events = []

    monkeypatch.setattr(
        "Scripts.run_heston_samples.restore_run_record",
        lambda **_kwargs: events.append("restore_run"),
    )
    monkeypatch.setattr(
        "Scripts.run_heston_samples.restore_sample",
        lambda *, sample_id, **_kwargs: events.append(f"restore_{sample_id}"),
    )
    monkeypatch.setattr(
        "Scripts.run_heston_samples.initialise_or_verify_run",
        lambda *_args, **_kwargs: events.append("initialise"),
    )
    monkeypatch.setattr(
        "Scripts.run_heston_samples.sync_run_record",
        lambda **_kwargs: events.append("sync_run"),
    )
    monkeypatch.setattr(
        "Scripts.run_heston_samples.sync_sample",
        lambda *, sample_id, **_kwargs: events.append(f"sync_{sample_id}"),
    )

    def run(task):
        events.append(f"worker_{task.sample_id}")
        return _outcome(task)

    monkeypatch.setattr("Scripts.run_heston_samples._run_one_sample", run)
    run_range(
        config_path=str(config),
        output_root=tmp_path / "run",
        sample_start=2,
        sample_end=4,
        sample_workers=1,
        generation_only=True,
        resume=True,
        s3_uri="s3://bucket/run",
        restore_from_s3=True,
    )

    assert events[:4] == ["restore_run", "restore_2", "restore_3", "initialise"]
    assert events.index("restore_3") < events.index("worker_2")
    assert "restore_0" not in events
    assert "restore_1" not in events
    assert "restore_4" not in events


def test_restore_requires_resume_and_s3(tmp_path):
    config = _config(tmp_path)
    with pytest.raises(ValueError, match="requires --resume"):
        run_range(
            config_path=str(config),
            output_root=tmp_path / "run_a",
            sample_start=0,
            sample_end=1,
            sample_workers=1,
            restore_from_s3=True,
            s3_uri="s3://bucket/run",
        )
    with pytest.raises(ValueError, match="requires --s3-uri"):
        run_range(
            config_path=str(config),
            output_root=tmp_path / "run_b",
            sample_start=0,
            sample_end=1,
            sample_workers=1,
            restore_from_s3=True,
            resume=True,
        )


def test_sample_persistence_failure_is_separate_and_returns_nonzero(
    tmp_path,
    monkeypatch,
):
    config = _config(tmp_path)
    root = tmp_path / "run"
    sentinel = root / "sample_000" / "record.json"
    sentinel.parent.mkdir(parents=True)
    sentinel.write_text('{"local": true}')
    before = sentinel.read_bytes()
    monkeypatch.setattr(
        "Scripts.run_heston_samples._run_one_sample",
        lambda task: _outcome(task),
    )
    monkeypatch.setattr(
        "Scripts.run_heston_samples.sync_run_record",
        lambda **_kwargs: None,
    )
    monkeypatch.setattr(
        "Scripts.run_heston_samples.sync_sample",
        lambda **_kwargs: (_ for _ in ()).throw(
            S3PersistenceError("upload failed")
        ),
    )

    summary = run_range(
        config_path=str(config),
        output_root=root,
        sample_start=0,
        sample_end=1,
        sample_workers=1,
        scenario="clean",
        resume=True,
        allow_dirty=True,
        s3_uri="s3://bucket/run",
    )

    assert summary.hard_failures == {}
    assert len(summary.persistence_failures) == 1
    assert summary.exit_code == 1
    assert sentinel.read_bytes() == before


def test_hard_computation_failure_is_not_uploaded_as_successful(
    tmp_path,
    monkeypatch,
):
    config = _config(tmp_path)
    sample_syncs = []
    monkeypatch.setattr(
        "Scripts.run_heston_samples._run_one_sample",
        lambda _task: (_ for _ in ()).throw(RuntimeError("numerical failure")),
    )
    monkeypatch.setattr(
        "Scripts.run_heston_samples.sync_run_record",
        lambda **_kwargs: None,
    )
    monkeypatch.setattr(
        "Scripts.run_heston_samples.sync_sample",
        lambda *, sample_id, **_kwargs: sample_syncs.append(sample_id),
    )

    summary = run_range(
        config_path=str(config),
        output_root=tmp_path / "run",
        sample_start=0,
        sample_end=1,
        sample_workers=1,
        scenario="clean",
        resume=True,
        allow_dirty=True,
        s3_uri="s3://bucket/run",
    )

    assert summary.exit_code == 1
    assert summary.hard_failures == {0: "RuntimeError: numerical failure"}
    assert summary.persistence_failures == []
    assert sample_syncs == []


def test_launcher_main_returns_nonzero_for_persistence_failure(monkeypatch):
    summary = RangeSummary(
        requested_sample_ids=(0,),
        requested_scenarios=("clean",),
        persistence_failures=["sample_000: upload failed"],
    )
    monkeypatch.setattr(
        "Scripts.run_heston_samples.run_range",
        lambda **_kwargs: summary,
    )

    assert main(
        [
            "--config",
            "experiment.json",
            "--output-root",
            "run",
            "--sample-start",
            "0",
            "--sample-end",
            "1",
        ]
    ) == 1
