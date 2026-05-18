#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$ROOT/Python/ImpliedVolatility/jaeckel/lets_be_rational"
OUT="$ROOT/Python/ImpliedVolatility/native"
PYTHON_BIN="${PYTHON_BIN:-python3}"
CXX_BIN="${CXX:-c++}"
export PYTHONPATH="$ROOT/Python${PYTHONPATH:+:$PYTHONPATH}"

mkdir -p "$OUT"

case "$(uname -s)" in
  Darwin)
    LIB="$OUT/liblets_be_rational.dylib"
    ARCH_FLAGS=()
    PY_ARCH="$("$PYTHON_BIN" -c 'import platform; print(platform.machine())')"
    if [[ "$PY_ARCH" == "x86_64" ]]; then
      ARCH_FLAGS=(-arch x86_64)
    elif [[ "$PY_ARCH" == "arm64" ]]; then
      ARCH_FLAGS=(-arch arm64)
    fi
    ;;
  *)
    LIB="$OUT/liblets_be_rational.so"
    ARCH_FLAGS=()
    ;;
esac

"$CXX_BIN" "${ARCH_FLAGS[@]}" \
  -std=c++17 -O3 -DNDEBUG -DNO_XL_API -shared -fPIC \
  "$SRC/erf_cody.cpp" \
  "$SRC/lets_be_rational.cpp" \
  "$SRC/normaldistribution.cpp" \
  "$SRC/rationalcubic.cpp" \
  -o "$LIB"

"$PYTHON_BIN" - <<'PY'
from ImpliedVolatility._lbr_backend import black, implied_black_volatility

price = black(1.0, 1.0, 1.0, 1.0, 1.0)
vol = implied_black_volatility(price, 1.0, 1.0, 1.0, 1.0)
if abs(vol - 1.0) > 1e-14:
    raise SystemExit(f"LetsBeRational smoke test failed: {vol}")
print(f"LetsBeRational smoke test ok: price={price:.17g}, vol={vol:.17g}")
PY
