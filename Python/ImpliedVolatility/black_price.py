from __future__ import annotations

from ImpliedVolatility._lbr_backend import black


def _option_sign(option_type: str) -> float:
    kind = option_type.lower()
    if kind == "call":
        return 1.0
    if kind == "put":
        return -1.0
    raise ValueError("option_type must be 'call' or 'put'")


def black76_price(
    *,
    forward: float,
    strike: float,
    tau: float,
    vol: float,
    discount_factor: float,
    option_type: str,
) -> float:
    if forward <= 0.0:
        raise ValueError("forward must be strictly positive")
    if strike <= 0.0:
        raise ValueError("strike must be strictly positive")
    if tau < 0.0:
        raise ValueError("tau must be non-negative")
    if vol < 0.0:
        raise ValueError("vol must be non-negative")
    if discount_factor <= 0.0:
        raise ValueError("discount_factor must be strictly positive")

    q = _option_sign(option_type)

    if tau == 0.0 or vol == 0.0:
        intrinsic = max(q * (forward - strike), 0.0)
        return discount_factor * intrinsic

    return discount_factor * black(forward, strike, vol, tau, q)
