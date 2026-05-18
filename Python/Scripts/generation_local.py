#!/usr/bin/env python3
"""
Local smoke-run wrapper for Chapter 3 clean generation.

Intended use after Codex implements the experiment modules specified in
codex_ch3_clean_generation.md.

Example:
    PYTHONPATH=Python python scripts/run_ch3_generation_local.py \
        --config Python/experiments/configs/ch3_clean_generation_run_001.json \
        --n-samples 2 \
        --workers 2 \
        --output-root outputs/ch3/local_smoke
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
    parser.add_argument("--n-samples", type=int, default=2)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--output-root", default="outputs/ch3/local_smoke")
    args = parser.parse_args()

    _set_thread_env()

    repo_root = Path.cwd()
    python_path = str(repo_root / "Python")
    os.environ["PYTHONPATH"] = python_path + os.pathsep + os.environ.get("PYTHONPATH", "")

    _run([
        sys.executable,
        "-m",
        "experiments.run_ch3_generation_local",
        "--config",
        args.config,
        "--n-samples",
        str(args.n_samples),
        "--workers",
        str(args.workers),
        "--output-root",
        args.output_root,
    ])

    _run([
        sys.executable,
        "-m",
        "experiments.validate_ch3_generation",
        "--run-root",
        args.output_root,
        "--expected-samples",
        str(args.n_samples),
    ])


if __name__ == "__main__":
    main()
