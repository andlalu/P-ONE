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
OUTPUT_ROOT="${OUTPUT_ROOT:-outputs/generation/local_run}"
N_SAMPLES="${N_SAMPLES:-2}"
WORKERS="${WORKERS:-2}"

mkdir -p "$OUTPUT_ROOT/logs"
PYTHON_BIN="$PYTHON_BIN" bash scripts/build_lets_be_rational.sh 2>&1 | tee "$OUTPUT_ROOT/logs/build_lets_be_rational.log"
"$PYTHON_BIN" -m Scripts.generation_local \
  --config "$CONFIG" \
  --n-samples "$N_SAMPLES" \
  --workers "$WORKERS" \
  --output-root "$OUTPUT_ROOT" \
  "$@" 2>&1 | tee "$OUTPUT_ROOT/logs/generation_local.log"
