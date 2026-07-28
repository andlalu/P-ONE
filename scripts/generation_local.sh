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
OUTPUT_ROOT="${OUTPUT_ROOT:-outputs/run_001}"
SAMPLE_ID="${SAMPLE_ID:-0}"

PYTHON_BIN="$PYTHON_BIN" bash scripts/build_lets_be_rational.sh
"$PYTHON_BIN" Python/Scripts/run_heston_sample.py \
  --config "$CONFIG" \
  --sample-id "$SAMPLE_ID" \
  --output-root "$OUTPUT_ROOT" \
  "$@"
