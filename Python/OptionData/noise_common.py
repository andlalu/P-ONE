from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_price

NOISE_SCENARIOS = (
    "low_iid",
    "spatial_corr",
    "persistent_factor",
    "variance_linked_factor",
)
NOISY_PANEL_EXTRA_COLUMNS = [
    "noise_scenario",
    "raw_noisy_iv",
    "raw_noisy_price",
    "observed_price",
    "observed_iv",
    "was_price_capped",
    "cap_direction",
    "noise_draw",
    "noise_seed",
]


@dataclass(frozen=True)
class NoiseSettings:
    """Shared controls and settings for the noise scenarios."""

    base_seed: int
    sigma_min: float
    price_epsilon: float
    scenarios: dict[str, dict[str, Any]]
    tick_size: float | None = None

    def scenario_names(self) -> tuple[str, ...]:
        return tuple(name for name in NOISE_SCENARIOS if name in self.scenarios)


def scenario_seed(base_seed: int, sample_id: int, scenario: str) -> int:
    offsets = {
        "low_iid": 101,
        "spatial_corr": 202,
        "persistent_factor": 303,
        "variance_linked_factor": 404,
    }
    if scenario not in offsets:
        raise ValueError(f"unknown noise scenario: {scenario}")
    return int(base_seed + 1000 * sample_id + offsets[scenario])


def marginal_scale(
    log_moneyness: np.ndarray,
    tau: np.ndarray,
    config: dict[str, Any],
) -> np.ndarray:
    """Set the noise scale from moneyness and maturity."""

    return float(config["alpha_0"]) * (
        1.0
        + float(config["alpha_m"]) * np.abs(log_moneyness)
        + float(config["alpha_tau"])
        / np.sqrt(np.maximum(tau, float(config["tau_min"])))
    )


def price_bounds(
    spot: float,
    strike: float,
    tau: float,
    rate: float,
    dividend_yield: float,
    option_type: str,
) -> tuple[float, float]:
    discounted_spot = spot * math.exp(-dividend_yield * tau)
    discounted_strike = strike * math.exp(-rate * tau)
    kind = option_type.lower()
    if kind == "call":
        return max(discounted_spot - discounted_strike, 0.0), discounted_spot
    if kind == "put":
        return max(discounted_strike - discounted_spot, 0.0), discounted_strike
    raise ValueError("option_type must be 'call' or 'put'")


def apply_price_projection(
    *,
    raw_price: float,
    spot: float,
    strike: float,
    tau: float,
    rate: float,
    dividend_yield: float,
    option_type: str,
    price_epsilon: float,
) -> tuple[float, bool, str]:
    """Keep a raw noisy price strictly inside its no-arbitrage bounds."""

    lower, upper = price_bounds(
        spot,
        strike,
        tau,
        rate,
        dividend_yield,
        option_type,
    )
    lower_cap = lower + price_epsilon
    upper_cap = upper - price_epsilon
    if lower_cap >= upper_cap:
        midpoint = 0.5 * (lower + upper)
        return midpoint, True, "lower" if raw_price <= midpoint else "upper"
    observed = min(upper_cap, max(lower_cap, raw_price))
    if observed <= lower_cap and raw_price < lower_cap:
        return observed, True, "lower"
    if observed >= upper_cap and raw_price > upper_cap:
        return observed, True, "upper"
    return observed, False, "none"


def apply_price_mechanics(
    *,
    raw_price: float,
    spot: float,
    strike: float,
    tau: float,
    rate: float,
    dividend_yield: float,
    option_type: str,
    tick_size: float,
    price_epsilon: float,
) -> tuple[float, float, bool, str]:
    """Run the original run_001 tick rounding followed by price projection."""

    rounded = tick_size * round(raw_price / tick_size)
    observed, was_capped, cap_direction = apply_price_projection(
        raw_price=rounded,
        spot=spot,
        strike=strike,
        tau=tau,
        rate=rate,
        dividend_yield=dividend_yield,
        option_type=option_type,
        price_epsilon=price_epsilon,
    )
    return rounded, observed, was_capped, cap_direction


def apply_noise_to_rows(
    rows: list[dict[str, Any]],
    *,
    scenario: str,
    seed: int,
    raw_noisy_iv: np.ndarray,
    noise_draw: np.ndarray,
    config: NoiseSettings,
) -> list[dict[str, Any]]:
    """Turn noisy IVs into bounded observed prices and IVs."""

    output: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        spot = float(row["S"])
        strike = float(row["strike"])
        tau = float(row["maturity_years"])
        if tau <= 0.0:
            raise ValueError("noisy panel generation requires strictly positive maturities")
        rate = float(row["r"])
        dividend_yield = float(row["q"])
        option_type = str(row["option_type"]).lower()
        forward = float(row["forward"])
        discount = math.exp(-rate * tau)
        raw_price = black76_price(
            forward=forward,
            strike=strike,
            tau=tau,
            vol=float(raw_noisy_iv[index]),
            discount_factor=discount,
            option_type=option_type,
        )
        rounded = None
        if config.tick_size is None:
            observed_price, was_capped, cap_direction = apply_price_projection(
                raw_price=raw_price,
                spot=spot,
                strike=strike,
                tau=tau,
                rate=rate,
                dividend_yield=dividend_yield,
                option_type=option_type,
                price_epsilon=config.price_epsilon,
            )
        else:
            rounded, observed_price, was_capped, cap_direction = apply_price_mechanics(
                raw_price=raw_price,
                spot=spot,
                strike=strike,
                tau=tau,
                rate=rate,
                dividend_yield=dividend_yield,
                option_type=option_type,
                tick_size=config.tick_size,
                price_epsilon=config.price_epsilon,
            )
        observed_iv = implied_vol_black76(
            price=observed_price,
            forward=forward,
            strike=strike,
            tau=tau,
            discount_factor=discount,
            option_type=option_type,
            on_bounds="clip",
        )
        item = dict(row)
        item.update(
            {
                "noise_scenario": scenario,
                "raw_noisy_iv": float(raw_noisy_iv[index]),
                "raw_noisy_price": float(raw_price),
                "raw_price_before_rounding": None if rounded is None else float(raw_price),
                "price_after_rounding": None if rounded is None else float(rounded),
                "observed_price": float(observed_price),
                "observed_iv": max(float(observed_iv), config.sigma_min),
                "was_price_capped": bool(was_capped),
                "cap_direction": cap_direction,
                "noise_draw": float(noise_draw[index]),
                "noise_seed": seed,
            }
        )
        output.append(item)
    return output
