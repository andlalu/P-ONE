from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

for key in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(key, "1")

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from Scripts.experiment_config import ExperimentConfig, load_experiment_config
from Scripts.s3_persistence import (
    S3PersistenceError,
    restore_run_record,
    restore_sample,
    sync_run_record,
    sync_sample,
    validate_s3_uri,
)
from Scripts.sample_run import (
    SCENARIO_ORDER,
    _estimate_is_complete,
    initialise_or_verify_run,
    run_sample,
)


@dataclass(frozen=True)
class SampleTask:
    config_path: str
    sample_id: int
    output_root: str
    resume: bool
    overwrite: bool
    generation_only: bool
    state_only: bool
    scenario: str | None
    max_dates: int | None
    log_level: str
    allow_dirty: bool


@dataclass(frozen=True)
class SampleOutcome:
    sample_id: int
    record: dict[str, Any]
    elapsed_seconds: float


@dataclass
class RangeSummary:
    requested_sample_ids: tuple[int, ...]
    requested_scenarios: tuple[str, ...]
    outcomes: dict[int, SampleOutcome] = field(default_factory=dict)
    hard_failures: dict[int, str] = field(default_factory=dict)
    persistence_failures: list[str] = field(default_factory=list)

    def estimate_counts(self) -> tuple[int, int]:
        succeeded = 0
        unsuccessful = 0
        for outcome in self.outcomes.values():
            estimates = outcome.record.get("estimation", {})
            for scenario in self.requested_scenarios:
                result = estimates.get(scenario)
                if _estimate_is_complete(result):
                    if result["success"]:
                        succeeded += 1
                    else:
                        unsuccessful += 1
        return succeeded, unsuccessful

    @property
    def exit_code(self) -> int:
        return int(bool(self.hard_failures or self.persistence_failures))


def _run_one_sample(task: SampleTask) -> SampleOutcome:
    started = time.perf_counter()
    config = load_experiment_config(task.config_path)
    record = run_sample(
        config=config,
        sample_id=task.sample_id,
        output_root=task.output_root,
        resume=task.resume,
        overwrite=task.overwrite,
        generation_only=task.generation_only,
        state_only=task.state_only,
        scenario=task.scenario,
        max_dates=task.max_dates,
        log_level=task.log_level,
        allow_dirty=task.allow_dirty,
    )
    return SampleOutcome(
        sample_id=task.sample_id,
        record=record,
        elapsed_seconds=time.perf_counter() - started,
    )


def _validate_options(
    *,
    config: ExperimentConfig,
    sample_start: int,
    sample_end: int,
    sample_workers: int,
    resume: bool,
    overwrite: bool,
    generation_only: bool,
    state_only: bool,
    scenario: str | None,
    max_dates: int | None,
    s3_uri: str | None,
    restore_from_s3: bool,
) -> None:
    if not 0 <= sample_start < sample_end <= config.n_samples:
        raise ValueError(
            "sample range must satisfy "
            f"0 <= sample-start < sample-end <= configured n_samples={config.n_samples}"
        )
    if sample_workers <= 0:
        raise ValueError("--sample-workers must be positive")
    if resume and overwrite:
        raise ValueError("--resume and --overwrite cannot be combined")
    if generation_only and state_only:
        raise ValueError("--generation-only and --state-only cannot be combined")
    if generation_only and scenario is not None:
        raise ValueError("--generation-only cannot be combined with --scenario")
    if max_dates is not None and max_dates <= 0:
        raise ValueError("--max-dates must be positive")
    if restore_from_s3 and not resume:
        raise ValueError("--restore-from-s3 requires --resume")
    if restore_from_s3 and overwrite:
        raise ValueError("--restore-from-s3 cannot be combined with --overwrite")
    if restore_from_s3 and s3_uri is None:
        raise ValueError("--restore-from-s3 requires --s3-uri")


def _task_for_sample(
    *,
    config_path: str,
    sample_id: int,
    output_root: Path,
    resume: bool,
    overwrite: bool,
    generation_only: bool,
    state_only: bool,
    scenario: str | None,
    max_dates: int | None,
    log_level: str,
    allow_dirty: bool,
) -> SampleTask:
    return SampleTask(
        config_path=str(Path(config_path).expanduser().resolve()),
        sample_id=sample_id,
        output_root=str(output_root),
        resume=resume,
        overwrite=overwrite,
        generation_only=generation_only,
        state_only=state_only,
        scenario=scenario,
        max_dates=max_dates,
        log_level=log_level,
        allow_dirty=allow_dirty,
    )


def _requested_scenarios(
    *, generation_only: bool, state_only: bool, scenario: str | None
) -> tuple[str, ...]:
    if generation_only or state_only:
        return ()
    return (scenario,) if scenario is not None else SCENARIO_ORDER


def _scenario_label(outcome: SampleOutcome, scenario: str | None) -> str:
    if scenario is None:
        return "estimates=all" if outcome.record.get("estimation") else "estimates=pending"
    result = outcome.record.get("estimation", {}).get(scenario)
    if not _estimate_is_complete(result):
        return f"{scenario}=pending"
    return f"{scenario}={'success' if result['success'] else 'success_false'}"


def _publish_outcome(
    *,
    outcome: SampleOutcome,
    summary: RangeSummary,
    output_root: Path,
    s3_uri: str | None,
    scenario: str | None,
) -> None:
    summary.outcomes[outcome.sample_id] = outcome
    persistence = "disabled"
    if s3_uri is not None:
        try:
            sync_sample(
                output_root=output_root,
                s3_uri=s3_uri,
                sample_id=outcome.sample_id,
            )
            persistence = "ok"
        except S3PersistenceError as exception:
            persistence = "failed"
            summary.persistence_failures.append(
                f"sample_{outcome.sample_id:03d}: {exception}"
            )
    print(
        f"sample={outcome.sample_id:03d} "
        f"status={outcome.record.get('status', 'unknown')} "
        f"{_scenario_label(outcome, scenario)} "
        f"persistence={persistence} "
        f"elapsed={outcome.elapsed_seconds:.3f}s"
    )


def _record_hard_failure(
    *, sample_id: int, exception: BaseException, summary: RangeSummary
) -> None:
    detail = f"{type(exception).__name__}: {exception}"
    summary.hard_failures[sample_id] = detail
    print(
        f"sample={sample_id:03d} status=failed "
        f"persistence=not_attempted error={detail}",
        file=sys.stderr,
    )


def run_range(
    *,
    config_path: str,
    output_root: str | Path,
    sample_start: int,
    sample_end: int,
    sample_workers: int | None = None,
    resume: bool = False,
    overwrite: bool = False,
    generation_only: bool = False,
    state_only: bool = False,
    scenario: str | None = None,
    max_dates: int | None = None,
    log_level: str = "INFO",
    allow_dirty: bool = False,
    s3_uri: str | None = None,
    restore_from_s3: bool = False,
) -> RangeSummary:
    config = load_experiment_config(config_path)
    workers = sample_workers if sample_workers is not None else (config.workers or 1)
    prefix = None if s3_uri is None else validate_s3_uri(s3_uri)
    _validate_options(
        config=config,
        sample_start=sample_start,
        sample_end=sample_end,
        sample_workers=workers,
        resume=resume,
        overwrite=overwrite,
        generation_only=generation_only,
        state_only=state_only,
        scenario=scenario,
        max_dates=max_dates,
        s3_uri=prefix,
        restore_from_s3=restore_from_s3,
    )

    root = Path(output_root).expanduser().resolve()
    sample_ids = tuple(range(sample_start, sample_end))
    summary = RangeSummary(
        requested_sample_ids=sample_ids,
        requested_scenarios=_requested_scenarios(
            generation_only=generation_only,
            state_only=state_only,
            scenario=scenario,
        ),
    )

    if restore_from_s3:
        assert prefix is not None
        restore_run_record(s3_uri=prefix, output_root=root)
        for sample_id in sample_ids:
            restore_sample(
                s3_uri=prefix,
                output_root=root,
                sample_id=sample_id,
            )

    initialise_or_verify_run(
        root,
        config,
        allow_git_mismatch=allow_dirty,
    )
    if prefix is not None:
        try:
            sync_run_record(output_root=root, s3_uri=prefix)
        except S3PersistenceError as exception:
            summary.persistence_failures.append(f"run.json before launch: {exception}")

    tasks = [
        _task_for_sample(
            config_path=config_path,
            sample_id=sample_id,
            output_root=root,
            resume=resume,
            overwrite=overwrite,
            generation_only=generation_only,
            state_only=state_only,
            scenario=scenario,
            max_dates=max_dates,
            log_level=log_level,
            allow_dirty=allow_dirty,
        )
        for sample_id in sample_ids
    ]

    if workers == 1:
        for task in tasks:
            try:
                outcome = _run_one_sample(task)
            except Exception as exception:
                _record_hard_failure(
                    sample_id=task.sample_id,
                    exception=exception,
                    summary=summary,
                )
            else:
                _publish_outcome(
                    outcome=outcome,
                    summary=summary,
                    output_root=root,
                    s3_uri=prefix,
                    scenario=scenario,
                )
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_run_one_sample, task): task.sample_id
                for task in tasks
            }
            for future in as_completed(futures):
                sample_id = futures[future]
                try:
                    outcome = future.result()
                except Exception as exception:
                    _record_hard_failure(
                        sample_id=sample_id,
                        exception=exception,
                        summary=summary,
                    )
                else:
                    _publish_outcome(
                        outcome=outcome,
                        summary=summary,
                        output_root=root,
                        s3_uri=prefix,
                        scenario=scenario,
                    )

    if prefix is not None:
        try:
            sync_run_record(output_root=root, s3_uri=prefix)
        except S3PersistenceError as exception:
            summary.persistence_failures.append(f"run.json after launch: {exception}")

    successful, unsuccessful = summary.estimate_counts()
    print(f"requested samples: {len(summary.requested_sample_ids)}")
    print(f"completed sample runs: {len(summary.outcomes)}")
    print(f"completed estimates with success=true: {successful}")
    print(f"completed estimates with success=false: {unsuccessful}")
    print(f"numerical hard failures: {len(summary.hard_failures)}")
    print(f"persistence failures: {len(summary.persistence_failures)}")
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run a range of complete Heston samples with sample-level process "
            "parallelism and optional parent-owned S3 persistence."
        )
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--sample-start", required=True, type=int)
    parser.add_argument("--sample-end", required=True, type=int)
    parser.add_argument("--sample-workers", type=int, default=None)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--resume", action="store_true")
    mode.add_argument("--overwrite", action="store_true")
    diagnostic = parser.add_mutually_exclusive_group()
    diagnostic.add_argument("--generation-only", action="store_true")
    diagnostic.add_argument("--state-only", action="store_true")
    parser.add_argument("--scenario", choices=SCENARIO_ORDER, default=None)
    parser.add_argument("--max-dates", type=int, default=None)
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    parser.add_argument("--allow-dirty", action="store_true")
    parser.add_argument("--s3-uri", default=None)
    parser.add_argument("--restore-from-s3", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        summary = run_range(
            config_path=arguments.config,
            output_root=arguments.output_root,
            sample_start=arguments.sample_start,
            sample_end=arguments.sample_end,
            sample_workers=arguments.sample_workers,
            resume=arguments.resume,
            overwrite=arguments.overwrite,
            generation_only=arguments.generation_only,
            state_only=arguments.state_only,
            scenario=arguments.scenario,
            max_dates=arguments.max_dates,
            log_level=arguments.log_level,
            allow_dirty=arguments.allow_dirty,
            s3_uri=arguments.s3_uri,
            restore_from_s3=arguments.restore_from_s3,
        )
    except ValueError as exception:
        parser.error(str(exception))
    except Exception as exception:
        print(
            f"sample range failed: {type(exception).__name__}: {exception}",
            file=sys.stderr,
        )
        return 1
    return summary.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
