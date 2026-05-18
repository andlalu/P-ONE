from __future__ import annotations

import numpy as np

from ImpliedVolatility._lbr_backend import implied_black_volatility


def _bounds(forward: float, strike: float, discount_factor: float, option_type: str) -> tuple[float, float]:
    if option_type == "call":
        return discount_factor * max(forward - strike, 0.0), discount_factor * forward
    if option_type == "put":
        return discount_factor * max(strike - forward, 0.0), discount_factor * strike
    raise ValueError("option_type must be 'call' or 'put'")


def _option_sign(option_type: str) -> float:
    if option_type == "call":
        return 1.0
    if option_type == "put":
        return -1.0
    raise ValueError("option_type must be 'call' or 'put'")


def implied_vol_black76(
    *,
    price: float,
    forward: float,
    strike: float,
    tau: float,
    discount_factor: float,
    option_type: str,
    on_bounds: str = "raise",
) -> float:
    """Invert Black-76 prices via Peter Jaeckel's LetsBeRational implementation."""
    kind = option_type.lower()
    if on_bounds not in {"raise", "nan", "clip"}:
        raise ValueError("on_bounds must be one of 'raise', 'nan', or 'clip'")
    if forward <= 0.0:
        raise ValueError("forward must be strictly positive")
    if strike <= 0.0:
        raise ValueError("strike must be strictly positive")
    if tau <= 0.0:
        raise ValueError("tau must be strictly positive for IV inversion")
    if discount_factor <= 0.0:
        raise ValueError("discount_factor must be strictly positive")

    lower, upper = _bounds(forward, strike, discount_factor, kind)
    eps = 1e-12 * max(1.0, upper)
    price_in = float(price)

    if price_in < lower - eps or price_in > upper + eps:
        if on_bounds == "raise":
            raise ValueError(
                f"price {price_in} outside Black-76 bounds [{lower}, {upper}] for {kind}"
            )
        if on_bounds == "nan":
            return float(np.nan)
        price_in = min(max(price_in, lower + eps), upper - eps)
    elif on_bounds == "clip":
        price_in = min(max(price_in, lower + eps), upper - eps)

    if price_in <= lower + eps:
        return 0.0

    return implied_black_volatility(
        price_in / discount_factor,
        forward,
        strike,
        tau,
        _option_sign(kind),
    )
