from __future__ import annotations

from typing import Any

import numpy as np

from OptionData.noise_common import (
    NOISE_SCENARIOS,
    NoiseSettings,
    apply_noise_to_rows,
)
from OptionData.noise_low_iid import low_iid_noise
from OptionData.noise_persistent_factor import (
    persistent_factor_noise,
    validate_persistent_factor_settings,
)
from OptionData.noise_spatial import spatial_noise
from OptionData.noise_variance_linked import (
    validate_variance_linked_settings,
    variance_linked_factor_noise,
)


def validate_noise_settings(config: NoiseSettings) -> None:
    if config.sigma_min <= 0.0 or config.price_epsilon <= 0.0:
        raise ValueError("noise sigma_min and price_epsilon must be positive")
    invalid = set(config.scenarios) - set(NOISE_SCENARIOS)
    if invalid:
        raise ValueError(f"unsupported noise scenarios: {', '.join(sorted(invalid))}")
    common_fields = ("alpha_0", "alpha_m", "alpha_tau", "tau_min")
    for name, scenario in config.scenarios.items():
        if any(float(scenario[field]) <= 0.0 for field in common_fields):
            raise ValueError(f"noise scenario {name} scale settings must be positive")
    if "spatial_corr" in config.scenarios:
        spatial = config.scenarios["spatial_corr"]
        if float(spatial["ell_m"]) <= 0.0 or float(spatial["ell_tau"]) <= 0.0:
            raise ValueError("spatial correlation lengths must be positive")
    if "persistent_factor" in config.scenarios:
        validate_persistent_factor_settings(
            config.scenarios["persistent_factor"]
        )
    if "variance_linked_factor" in config.scenarios:
        validate_variance_linked_settings(
            config.scenarios["variance_linked_factor"]
        )


def generate_noisy_panel_rows(
    clean_rows: list[dict[str, Any]],
    *,
    scenario: str,
    seed: int,
    config: NoiseSettings,
    params_p: Any | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Apply one configured noise scenario."""

    if scenario not in NOISE_SCENARIOS:
        raise ValueError(f"unsupported noise scenario: {scenario}")
    if scenario not in config.scenarios:
        raise ValueError(f"noise scenario is not configured: {scenario}")
    rows = sorted(
        clean_rows,
        key=lambda row: (
            int(row["week_index"]),
            float(row["maturity_years"]),
            float(row["log_moneyness"]),
        ),
    )
    rng = np.random.default_rng(seed)
    if scenario == "low_iid":
        raw_noisy_iv, noise, factors = low_iid_noise(rows, rng, config)
    elif scenario == "spatial_corr":
        raw_noisy_iv, noise, factors = spatial_noise(rows, rng, config)
    elif scenario == "persistent_factor":
        raw_noisy_iv, noise, factors = persistent_factor_noise(rows, rng, config)
    else:
        if params_p is None:
            raise ValueError("variance_linked_factor requires physical Heston parameters")
        raw_noisy_iv, noise, factors = variance_linked_factor_noise(
            rows,
            rng,
            config,
            params_p,
        )
    return (
        apply_noise_to_rows(
            rows,
            scenario=scenario,
            seed=seed,
            raw_noisy_iv=raw_noisy_iv,
            noise_draw=noise,
            config=config,
        ),
        factors,
    )
