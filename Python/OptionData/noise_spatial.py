from __future__ import annotations

from typing import Any

import numpy as np

from OptionData.noise_common import NoiseSettings, marginal_scale


def spatial_noise(
    rows: list[dict[str, Any]],
    rng: np.random.Generator,
    config: NoiseSettings,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Add spatially correlated IV noise within each date."""

    raw_noisy_iv = np.array([float(row["model_iv"]) for row in rows])
    noise = np.zeros(len(rows), dtype=float)
    for week in sorted({int(row["week_index"]) for row in rows}):
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
        spatial = config.scenarios["spatial_corr"]
        scale = marginal_scale(log_moneyness, maturities, spatial)
        distance = (
            np.abs(log_moneyness[:, None] - log_moneyness[None, :])
            / float(spatial["ell_m"])
            + np.abs(maturities[:, None] - maturities[None, :])
            / float(spatial["ell_tau"])
        )
        correlation = np.exp(-distance)
        jitter = float(spatial["correlation_jitter"])
        while True:
            try:
                factor = np.linalg.cholesky(
                    correlation + jitter * np.eye(len(indices))
                )
                break
            except np.linalg.LinAlgError:
                jitter *= 10.0
                if jitter > float(spatial["max_correlation_jitter"]):
                    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
                    factor = eigenvectors @ np.diag(
                        np.sqrt(np.maximum(eigenvalues, 0.0))
                    )
                    break
        noise[indices] = scale * (
            factor @ rng.standard_normal(len(indices))
        )
    return np.maximum(config.sigma_min, raw_noisy_iv + noise), noise, []
