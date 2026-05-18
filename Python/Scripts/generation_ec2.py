#!/usr/bin/env python3
"""
EC2 master wrapper for Chapter 3 clean generation.

This script is intentionally thin. The actual implementation should live in
Python/experiments/run_ch3_generation_ec2.py and use multiprocessing across
sample_id values.

Example:
    PYTHONPATH=Python python scripts/run_ch3_generation_ec2.py \
        --config Python/experiments/configs/ch3_clean_generation_run_001.json \
        --workers 15 \
        --output-root outputs/ch3/run_001

Optional S3 sync:
    aws s3 sync outputs/ch3/run_001 s3://<bucket>/p-one/ch3/run_001
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def _set_thread_env() -> None:
    for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(key, "1")


def _run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--output-root", default=None)
    parser.add_argument("--sample-start", type=int, default=0)
    parser.add_argument("--sample-end", type=int, default=None)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--s3-prefix", default=None)
    args = parser.parse_args()

    _set_thread_env()

    repo_root = Path.cwd()
    python_path = str(repo_root / "Python")
    os.environ["PYTHONPATH"] = python_path + os.pathsep + os.environ.get("PYTHONPATH", "")

    cmd = [
        sys.executable,
        "-m",
        "experiments.run_ch3_generation_ec2",
        "--config",
        args.config,
        "--sample-start",
        str(args.sample_start),
    ]

    if args.workers is not None:
        cmd += ["--workers", str(args.workers)]
    if args.output_root is not None:
        cmd += ["--output-root", args.output_root]
    if args.sample_end is not None:
        cmd += ["--sample-end", str(args.sample_end)]
    if args.skip_existing:
        cmd += ["--skip-existing"]

    _run(cmd)

    run_root = args.output_root or "outputs/ch3/run_001"
    _run([
        sys.executable,
        "-m",
        "experiments.validate_ch3_generation",
        "--run-root",
        run_root,
        "--expected-samples",
        "100",
    ])

    if args.s3_prefix:
        _run(["aws", "s3", "sync", run_root, args.s3_prefix])


if __name__ == "__main__":
    main()
