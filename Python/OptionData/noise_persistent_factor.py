from __future__ import annotations

import math
from typing import Any, Iterable

import numpy as np

from OptionData.noise_common import NoiseSettings, marginal_scale

PERSISTENT_FACTOR_DIM = 3
PERSISTENT_FACTOR_Q_RTOL = 1e-8
PERSISTENT_FACTOR_Q_ATOL = 1e-14
PERSISTENT_FACTOR_RESIDUAL_POLICIES = (
    "match_total_marginal_scale",
    "legacy_multiplier",
)


def compute_q_diag_from_stationary_std(
    a_diag: np.ndarray | Iterable[float],
    stationary_factor_std: np.ndarray | Iterable[float],
) -> np.ndarray:
    a = np.asarray(tuple(a_diag), dtype=float)
    stationary_std = np.asarray(tuple(stationary_factor_std), dtype=float)
    return (1.0 - a * a) * stationary_std * stationary_std


def compute_stationary_std_from_q_diag(
    a_diag: np.ndarray | Iterable[float],
    q_diag: np.ndarray | Iterable[float],
) -> np.ndarray:
    a = np.asarray(tuple(a_diag), dtype=float)
    q = np.asarray(tuple(q_diag), dtype=float)
    return np.sqrt(q / (1.0 - a * a))


def validate_persistent_factor_settings(config: dict[str, Any]) -> None:
    a_diag = np.asarray(config["a_diag"], dtype=float)
    if a_diag.shape != (PERSISTENT_FACTOR_DIM,):
        raise ValueError("persistent_factor.a_diag must have length 3")
    if not np.all(np.isfinite(a_diag)):
        raise ValueError("persistent_factor.a_diag must contain finite values")
    if np.any(np.abs(a_diag) >= 1.0):
        raise ValueError(
            "persistent_factor.a_diag entries must have absolute value strictly below 1"
        )

    stationary_std = None
    if config.get("stationary_factor_std") is not None:
        stationary_std = np.asarray(config["stationary_factor_std"], dtype=float)
        if stationary_std.shape != (PERSISTENT_FACTOR_DIM,):
            raise ValueError(
                "persistent_factor.stationary_factor_std must have length 3"
            )
        if not np.all(np.isfinite(stationary_std)):
            raise ValueError(
                "persistent_factor.stationary_factor_std must contain finite values"
            )
        if np.any(stationary_std < 0.0):
            raise ValueError(
                "persistent_factor.stationary_factor_std entries must be non-negative"
            )

    q_diag = None
    if config.get("q_diag") is not None:
        q_diag = np.asarray(config["q_diag"], dtype=float)
        if q_diag.shape != (PERSISTENT_FACTOR_DIM,):
            raise ValueError("persistent_factor.q_diag must have length 3")
        if not np.all(np.isfinite(q_diag)):
            raise ValueError("persistent_factor.q_diag must contain finite values")
        if np.any(q_diag < 0.0):
            raise ValueError(
                "persistent_factor.q_diag entries must be non-negative innovation variances"
            )

    if stationary_std is None and q_diag is None:
        raise ValueError(
            "persistent_factor requires stationary_factor_std or q_diag"
        )
    if stationary_std is not None:
        implied_q_diag = compute_q_diag_from_stationary_std(a_diag, stationary_std)
        if q_diag is None:
            q_diag = implied_q_diag
        elif not np.allclose(
            q_diag,
            implied_q_diag,
            rtol=PERSISTENT_FACTOR_Q_RTOL,
            atol=PERSISTENT_FACTOR_Q_ATOL,
        ):
            raise ValueError(
                "persistent_factor.q_diag must be the innovation covariance "
                "diagonal implied by a_diag and stationary_factor_std"
            )

    residual_policy = str(
        config.get("residual_policy", "match_total_marginal_scale")
    )
    if residual_policy not in PERSISTENT_FACTOR_RESIDUAL_POLICIES:
        raise ValueError(
            "persistent_factor.residual_policy must be one of "
            f"{', '.join(PERSISTENT_FACTOR_RESIDUAL_POLICIES)}"
        )
    if residual_policy == "legacy_multiplier":
        multiplier = config.get("residual_scale_multiplier")
        if multiplier is None:
            raise ValueError(
                "persistent_factor.residual_scale_multiplier is required "
                "for legacy_multiplier policy"
            )
        if not math.isfinite(float(multiplier)) or float(multiplier) < 0.0:
            raise ValueError(
                "persistent_factor.residual_scale_multiplier must be non-negative"
            )


def persistent_factor_residual_scale(
    log_moneyness: np.ndarray,
    tau: np.ndarray,
    config: dict[str, Any],
) -> np.ndarray:
    scale = marginal_scale(log_moneyness, tau, config)
    if config["residual_policy"] == "match_total_marginal_scale":
        stationary_std = np.array(config["stationary_factor_std"], dtype=float)
        factor_variance = (
            stationary_std[0] * stationary_std[0]
            + log_moneyness
            * log_moneyness
            * stationary_std[1]
            * stationary_std[1]
            + tau * tau * stationary_std[2] * stationary_std[2]
        )
        return np.sqrt(np.maximum(scale * scale - factor_variance, 0.0))
    if config["residual_policy"] == "legacy_multiplier":
        return float(config["residual_scale_multiplier"]) * scale
    raise ValueError(
        f"unsupported persistent_factor residual_policy: {config['residual_policy']}"
    )


def persistent_factor_noise(
    rows: list[dict[str, Any]],
    rng: np.random.Generator,
    config: NoiseSettings,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Add persistent three-factor and residual IV noise."""

    clean_iv = np.array([float(row["model_iv"]) for row in rows])
    noise = np.zeros(len(rows), dtype=float)
    factors: list[dict[str, Any]] = []
    factor_config = config.scenarios["persistent_factor"]
    validate_persistent_factor_settings(factor_config)
    if factor_config["factor_initialization"] != "zero":
        raise NotImplementedError("only zero factor initialization is implemented")
    a_diag = np.array(factor_config["a_diag"], dtype=float)
    q_diag = np.array(factor_config["q_diag"], dtype=float)
    factor = np.zeros(3, dtype=float)
    for week in sorted({int(row["week_index"]) for row in rows}):
        # The three factors follow f_t = A f_(t-1) + innovation_t.
        factor = a_diag * factor + rng.normal(
            loc=0.0,
            scale=np.sqrt(q_diag),
            size=3,
        )
        factors.append(
            {
                "week_index": week,
                "factor_0": factor[0],
                "factor_1": factor[1],
                "factor_2": factor[2],
            }
        )
        indices = np.array(
            [index for index, row in enumerate(rows) if int(row["week_index"]) == week],
            dtype=int,
        )
        log_moneyness = np.array(
            [float(rows[index]["log_moneyness"]) for index in indices]
        )
        maturities = np.array(
            [float(rows[index]["maturity_years"]) for index in indices]
        )
        residual_scale = persistent_factor_residual_scale(
            log_moneyness,
            maturities,
            factor_config,
        )
        factor_component = (
            factor[0] + log_moneyness * factor[1] + maturities * factor[2]
        )
        residual = residual_scale * rng.standard_normal(len(indices))
        noise[indices] = factor_component + residual
    return np.maximum(config.sigma_min, clean_iv + noise), noise, factors
