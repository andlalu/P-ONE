from __future__ import annotations

import ctypes
import platform
from functools import lru_cache
from pathlib import Path


def _library_name() -> str:
    system = platform.system()
    if system == "Darwin":
        return "liblets_be_rational.dylib"
    if system == "Windows":
        return "lets_be_rational.dll"
    return "liblets_be_rational.so"


def _library_path() -> Path:
    return Path(__file__).resolve().parent / "native" / _library_name()


@lru_cache(maxsize=1)
def load_library() -> ctypes.CDLL:
    library_path = _library_path()
    if not library_path.exists():
        raise ImportError(
            "Peter Jaeckel LetsBeRational native library is not built. "
            "Run `bash scripts/build_lets_be_rational.sh` from the repo root."
        )

    lib = ctypes.CDLL(str(library_path))
    for name in ("Black", "ImpliedBlackVolatility"):
        fn = getattr(lib, name)
        fn.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        fn.restype = ctypes.c_double
    return lib


def black(forward: float, strike: float, vol: float, tau: float, q: float) -> float:
    return float(load_library().Black(forward, strike, vol, tau, q))


def implied_black_volatility(price: float, forward: float, strike: float, tau: float, q: float) -> float:
    return float(load_library().ImpliedBlackVolatility(price, forward, strike, tau, q))
