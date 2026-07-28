from __future__ import annotations

from typing import Any

import numpy as np

from OptionData.noise_common import NoiseSettings, marginal_scale


def low_iid_noise(
    rows: list[dict[str, Any]],
    rng: np.random.Generator,
    config: NoiseSettings,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Add independent Gaussian IV noise."""

    log_moneyness = np.array([float(row["log_moneyness"]) for row in rows])
    maturities = np.array([float(row["maturity_years"]) for row in rows])
    clean_iv = np.array([float(row["model_iv"]) for row in rows])
    scale = marginal_scale(
        log_moneyness,
        maturities,
        config.scenarios["low_iid"],
    )
    noise = scale * rng.standard_normal(len(rows))
    return np.maximum(config.sigma_min, clean_iv + noise), noise, []
