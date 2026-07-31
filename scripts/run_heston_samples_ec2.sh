#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if [ -d ".venv" ]; then
  # shellcheck disable=SC1091
  source ".venv/bin/activate"
fi

export PYTHONPATH="$ROOT/Python${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

PYTHON_BIN="${PYTHON_BIN:-python3}"
CONFIG="${CONFIG:-Python/Scripts/configs/heston_experiment_run_001.json}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/data/p-one/outputs/run_001}"
: "${SAMPLE_START:?SAMPLE_START is required}"
: "${SAMPLE_END:?SAMPLE_END is required}"

arguments=(
  --config "$CONFIG"
  --output-root "$OUTPUT_ROOT"
  --sample-start "$SAMPLE_START"
  --sample-end "$SAMPLE_END"
  --resume
)

if [ -n "${SAMPLE_WORKERS:-}" ]; then
  arguments+=(--sample-workers "$SAMPLE_WORKERS")
fi
if [ -n "${SCENARIO:-}" ]; then
  arguments+=(--scenario "$SCENARIO")
fi
if [ -n "${S3_URI:-}" ]; then
  arguments+=(--s3-uri "$S3_URI")
fi
if [ "${RESTORE_FROM_S3:-0}" = "1" ]; then
  arguments+=(--restore-from-s3)
fi
if [ "${GENERATION_ONLY:-0}" = "1" ]; then
  arguments+=(--generation-only)
fi

PYTHON_BIN="$PYTHON_BIN" bash scripts/build_lets_be_rational.sh
"$PYTHON_BIN" Python/Scripts/run_heston_samples.py "${arguments[@]}" "$@"
