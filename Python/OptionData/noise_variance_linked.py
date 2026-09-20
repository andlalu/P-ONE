from __future__ import annotations

import math
from typing import Any

import numpy as np

from Models.Heston.parameters import HestonPhysicalParameters
from OptionData.noise_common import NoiseSettings, marginal_scale
from OptionData.noise_persistent_factor import validate_persistent_factor_settings

RESIDUAL_VARIANCE_TOLERANCE = 1e-15


def validate_variance_linked_settings(config: dict[str, Any]) -> None:
    validate_persistent_factor_settings(config)
    share = float(config["latent_variance_share"])
    if not math.isfinite(share) or not 0.0 < share < 1.0:
        raise ValueError(
            "variance_linked_factor.latent_variance_share must be finite and in (0, 1)"
        )
    if config["residual_policy"] != "match_total_marginal_scale":
        raise ValueError(
            "variance_linked_factor requires residual_policy='match_total_marginal_scale'"
        )


def variance_linked_gamma(
    log_moneyness: np.ndarray,
    tau: np.ndarray,
    config: dict[str, Any],
) -> float:
    contracts = np.unique(
        np.column_stack(
            (
                np.asarray(log_moneyness, dtype=float),
                np.asarray(tau, dtype=float),
            )
        ),
        axis=0,
    )
    scale = marginal_scale(contracts[:, 0], contracts[:, 1], config)
    share = float(config["latent_variance_share"])
    return float(math.sqrt(share / float(np.mean(1.0 / (scale * scale)))))


def variance_linked_residual_scale(
    log_moneyness: np.ndarray,
    tau: np.ndarray,
    config: dict[str, Any],
    gamma_v: float,
) -> np.ndarray:
    scale = marginal_scale(log_moneyness, tau, config)
    stationary_std = np.asarray(config["stationary_factor_std"], dtype=float)
    factor_variance = (
        stationary_std[0] * stationary_std[0]
        + log_moneyness
        * log_moneyness
        * stationary_std[1]
        * stationary_std[1]
        + tau * tau * stationary_std[2] * stationary_std[2]
    )
    residual_variance = scale * scale - factor_variance - gamma_v * gamma_v
    if np.any(residual_variance < -RESIDUAL_VARIANCE_TOLERANCE):
        raise ValueError(
            "variance_linked_factor calibration produces negative residual variance"
        )
    return np.sqrt(np.maximum(residual_variance, 0.0))


def variance_linked_factor_noise(
    rows: list[dict[str, Any]],
    rng: np.random.Generator,
    config: NoiseSettings,
    params_p: HestonPhysicalParameters,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Add persistent factors, latent-variance linkage and residual IV noise."""

    factor_config = config.scenarios["variance_linked_factor"]
    validate_variance_linked_settings(factor_config)
    if factor_config["factor_initialization"] != "zero":
        raise NotImplementedError("only zero factor initialization is implemented")
    variance_state_std = math.sqrt(
        params_p.vbar * params_p.sigma_v * params_p.sigma_v / (2.0 * params_p.kappa)
    )
    if not math.isfinite(variance_state_std) or variance_state_std <= 0.0:
        raise ValueError("variance_linked_factor requires valid physical Heston parameters")

    clean_iv = np.array([float(row["model_iv"]) for row in rows])
    all_log_moneyness = np.array([float(row["log_moneyness"]) for row in rows])
    all_maturities = np.array([float(row["maturity_years"]) for row in rows])
    gamma_v = variance_linked_gamma(
        all_log_moneyness,
        all_maturities,
        factor_config,
    )
    noise = np.zeros(len(rows), dtype=float)
    factors: list[dict[str, Any]] = []
    a_diag = np.array(factor_config["a_diag"], dtype=float)
    q_diag = np.array(factor_config["q_diag"], dtype=float)
    factor = np.zeros(3, dtype=float)
    for week in sorted({int(row["week_index"]) for row in rows}):
        factor = a_diag * factor + rng.normal(
            loc=0.0,
            scale=np.sqrt(q_diag),
            size=3,
        )
        indices = np.array(
            [index for index, row in enumerate(rows) if int(row["week_index"]) == week],
            dtype=int,
        )
        variance = np.array([float(rows[index]["V"]) for index in indices])
        if not np.allclose(variance, variance[0], rtol=0.0, atol=1e-15):
            raise ValueError(
                "variance_linked_factor requires latent variance to be constant within date"
            )
        variance_factor = gamma_v * (variance[0] - params_p.vbar) / variance_state_std
        factors.append(
            {
                "week_index": week,
                "factor_0": factor[0],
                "factor_1": factor[1],
                "factor_2": factor[2],
                "variance_linked_factor": variance_factor,
            }
        )
        log_moneyness = all_log_moneyness[indices]
        maturities = all_maturities[indices]
        residual_scale = variance_linked_residual_scale(
            log_moneyness,
            maturities,
            factor_config,
            gamma_v,
        )
        factor_component = (
            factor[0]
            + log_moneyness * factor[1]
            + maturities * factor[2]
            + variance_factor
        )
        residual = residual_scale * rng.standard_normal(len(indices))
        noise[indices] = factor_component + residual
    return np.maximum(config.sigma_min, clean_iv + noise), noise, factors
