from __future__ import annotations

import subprocess
from pathlib import Path
from urllib.parse import urlsplit


class S3PersistenceError(RuntimeError):
    """An AWS CLI persistence operation failed."""


def validate_s3_uri(uri: str) -> str:
    value = uri.strip().rstrip("/")
    parsed = urlsplit(value)
    if (
        parsed.scheme != "s3"
        or not parsed.netloc
        or parsed.query
        or parsed.fragment
        or parsed.username is not None
        or parsed.password is not None
        or parsed.port is not None
    ):
        raise ValueError(f"invalid S3 URI: {uri!r}")
    return value


def _missing_remote(stderr: str, stdout: str) -> bool:
    message = f"{stderr}\n{stdout}".lower()
    return any(
        marker in message
        for marker in (
            "nosuchkey",
            "not found",
            "does not exist",
            "status code: 404",
            "(404)",
        )
    )


def _run_aws(command: list[str], *, missing_ok: bool = False) -> bool:
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as exception:
        raise S3PersistenceError(
            f"could not run AWS CLI for {' '.join(command[2:4])}: {exception}"
        ) from exception
    if completed.returncode == 0:
        return True
    if missing_ok and _missing_remote(completed.stderr, completed.stdout):
        return False
    detail = completed.stderr.strip() or completed.stdout.strip() or "no error detail"
    raise S3PersistenceError(
        f"AWS CLI failed with exit code {completed.returncode}: {detail}"
    )


def restore_run_record(*, s3_uri: str, output_root: Path) -> bool:
    prefix = validate_s3_uri(s3_uri)
    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    return _run_aws(
        [
            "aws",
            "s3",
            "cp",
            f"{prefix}/run.json",
            str(root / "run.json"),
            "--only-show-errors",
        ],
        missing_ok=True,
    )


def restore_sample(*, s3_uri: str, output_root: Path, sample_id: int) -> bool:
    if sample_id < 0:
        raise ValueError("sample_id must be non-negative")
    prefix = validate_s3_uri(s3_uri)
    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    return _run_aws(
        [
            "aws",
            "s3",
            "sync",
            f"{prefix}/sample_{sample_id:03d}/",
            f"{root / f'sample_{sample_id:03d}'}/",
            "--exclude",
            "*.tmp*",
            "--only-show-errors",
        ],
        missing_ok=True,
    )


def sync_run_record(*, output_root: Path, s3_uri: str) -> None:
    prefix = validate_s3_uri(s3_uri)
    root = Path(output_root).expanduser().resolve()
    run_record = root / "run.json"
    if not run_record.is_file():
        raise S3PersistenceError(f"local run record does not exist: {run_record}")
    _run_aws(
        [
            "aws",
            "s3",
            "cp",
            str(run_record),
            f"{prefix}/run.json",
            "--only-show-errors",
        ]
    )


def sync_sample(*, output_root: Path, s3_uri: str, sample_id: int) -> None:
    if sample_id < 0:
        raise ValueError("sample_id must be non-negative")
    prefix = validate_s3_uri(s3_uri)
    root = Path(output_root).expanduser().resolve()
    sample = root / f"sample_{sample_id:03d}"
    if not sample.is_dir():
        raise S3PersistenceError(f"local sample directory does not exist: {sample}")
    _run_aws(
        [
            "aws",
            "s3",
            "sync",
            f"{sample}/",
            f"{prefix}/sample_{sample_id:03d}/",
            "--exclude",
            "*.tmp*",
            "--only-show-errors",
        ]
    )
