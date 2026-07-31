from __future__ import annotations

import subprocess

import pytest

from Scripts.s3_persistence import (
    S3PersistenceError,
    restore_run_record,
    restore_sample,
    sync_run_record,
    sync_sample,
    validate_s3_uri,
)


def _completed(command, returncode=0, stdout="", stderr=""):
    return subprocess.CompletedProcess(command, returncode, stdout, stderr)


def test_validate_s3_uri_normalises_prefix_and_rejects_non_s3_values():
    assert validate_s3_uri(" s3://bucket/p-one/run_001/ ") == (
        "s3://bucket/p-one/run_001"
    )
    for value in ("", "https://bucket/key", "s3:///key", "s3://bucket/key?x=1"):
        with pytest.raises(ValueError, match="invalid S3 URI"):
            validate_s3_uri(value)


def test_upload_commands_use_cp_and_sync_without_delete(tmp_path, monkeypatch):
    root = tmp_path / "run_001"
    sample = root / "sample_000"
    sample.mkdir(parents=True)
    (root / "run.json").write_text("{}")
    (sample / "record.json").write_text("{}")
    commands = []

    def record(command, **_kwargs):
        commands.append(command)
        return _completed(command)

    monkeypatch.setattr("Scripts.s3_persistence.subprocess.run", record)
    sync_run_record(output_root=root, s3_uri="s3://bucket/p-one/run_001")
    sync_sample(
        output_root=root,
        s3_uri="s3://bucket/p-one/run_001",
        sample_id=0,
    )

    assert commands == [
        [
            "aws",
            "s3",
            "cp",
            str(root.resolve() / "run.json"),
            "s3://bucket/p-one/run_001/run.json",
            "--only-show-errors",
        ],
        [
            "aws",
            "s3",
            "sync",
            f"{sample.resolve()}/",
            "s3://bucket/p-one/run_001/sample_000/",
            "--exclude",
            "*.tmp*",
            "--only-show-errors",
        ],
    ]
    assert all("--delete" not in command for command in commands)


def test_restore_commands_target_one_run_record_and_one_sample(tmp_path, monkeypatch):
    commands = []

    def record(command, **_kwargs):
        commands.append(command)
        return _completed(command)

    monkeypatch.setattr("Scripts.s3_persistence.subprocess.run", record)
    assert restore_run_record(
        s3_uri="s3://bucket/p-one/run_001",
        output_root=tmp_path,
    )
    assert restore_sample(
        s3_uri="s3://bucket/p-one/run_001",
        output_root=tmp_path,
        sample_id=7,
    )
    assert commands == [
        [
            "aws",
            "s3",
            "cp",
            "s3://bucket/p-one/run_001/run.json",
            str(tmp_path.resolve() / "run.json"),
            "--only-show-errors",
        ],
        [
            "aws",
            "s3",
            "sync",
            "s3://bucket/p-one/run_001/sample_007/",
            f"{tmp_path.resolve() / 'sample_007'}/",
            "--exclude",
            "*.tmp*",
            "--only-show-errors",
        ],
    ]
    assert all("--delete" not in command for command in commands)


def test_missing_remote_restore_is_absent_but_access_and_network_errors_fail(
    tmp_path,
    monkeypatch,
):
    responses = iter(
        [
            _completed([], 1, stderr="fatal error: (404) Not Found"),
            _completed([], 1, stderr="An error occurred (AccessDenied)"),
            _completed([], 1, stderr="Could not connect to the endpoint URL"),
        ]
    )
    monkeypatch.setattr(
        "Scripts.s3_persistence.subprocess.run",
        lambda *_args, **_kwargs: next(responses),
    )
    assert not restore_run_record(s3_uri="s3://bucket/run", output_root=tmp_path)
    with pytest.raises(S3PersistenceError, match="AccessDenied"):
        restore_sample(
            s3_uri="s3://bucket/run",
            output_root=tmp_path,
            sample_id=0,
        )
    with pytest.raises(S3PersistenceError, match="Could not connect"):
        restore_sample(
            s3_uri="s3://bucket/run",
            output_root=tmp_path,
            sample_id=1,
        )


def test_upload_failure_leaves_local_ebs_files_untouched(tmp_path, monkeypatch):
    sample = tmp_path / "sample_000"
    sample.mkdir()
    record = sample / "record.json"
    record.write_text('{"status": "generated"}')
    before = record.read_bytes()
    monkeypatch.setattr(
        "Scripts.s3_persistence.subprocess.run",
        lambda command, **_kwargs: _completed(
            command,
            1,
            stderr="network failure",
        ),
    )

    with pytest.raises(S3PersistenceError, match="network failure"):
        sync_sample(
            output_root=tmp_path,
            s3_uri="s3://bucket/run",
            sample_id=0,
        )
    assert record.read_bytes() == before
