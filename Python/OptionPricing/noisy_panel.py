from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_price
from OptionPricing.clean_panel import parquet_available

NOISE_SCENARIOS = ("low_iid", "spatial_corr", "persistent_factor")
PERSISTENT_FACTOR_DIM = 3
PERSISTENT_FACTOR_Q_RTOL = 1e-8
PERSISTENT_FACTOR_Q_ATOL = 1e-14
PERSISTENT_FACTOR_RESIDUAL_POLICIES = ("match_total_marginal_scale", "legacy_multiplier")
OBSERVED_PANEL_EXTRA_COLUMNS = [
    "noise_scenario",
    "raw_contaminated_iv",
    "raw_price_before_rounding",
    "price_after_rounding",
    "observed_price",
    "observed_iv",
    "was_price_capped",
    "cap_direction",
    "noise_draw",
    "noise_seed",
]


@dataclass(frozen=True)
class MarginalNoiseScaleConfig:
    alpha_0: float
    alpha_m: float
    alpha_tau: float
    tau_min: float


@dataclass(frozen=True)
class LowIidNoiseConfig(MarginalNoiseScaleConfig):
    enabled: bool = True


@dataclass(frozen=True)
class SpatialCorrNoiseConfig(MarginalNoiseScaleConfig):
    ell_m: float = 0.10
    ell_tau: float = 0.25
    correlation_jitter: float = 1e-10
    max_correlation_jitter: float = 1e-6
    enabled: bool = True


@dataclass(frozen=True)
class PersistentFactorNoiseConfig(MarginalNoiseScaleConfig):
    a_diag: tuple[float, float, float] = (0.85, 0.65, 0.65)
    stationary_factor_std: tuple[float, float, float] | None = (0.0007, 0.0025, 0.0008)
    q_diag: tuple[float, float, float] | None = (1.35975e-7, 3.609375e-6, 3.696e-7)
    residual_policy: str = "match_total_marginal_scale"
    residual_scale_multiplier: float | None = None
    factor_initialization: str = "zero"
    enabled: bool = True


@dataclass(frozen=True)
class NoiseGenerationConfig:
    base_seed: int
    sigma_min: float
    price_epsilon: float
    tick_size: float
    low_iid: LowIidNoiseConfig
    spatial_corr: SpatialCorrNoiseConfig
    persistent_factor: PersistentFactorNoiseConfig

    def enabled_scenarios(self) -> tuple[str, ...]:
        out: list[str] = []
        if self.low_iid.enabled:
            out.append("low_iid")
        if self.spatial_corr.enabled:
            out.append("spatial_corr")
        if self.persistent_factor.enabled:
            out.append("persistent_factor")
        return tuple(out)


@dataclass(frozen=True)
class NoisyPanelResult:
    sample_id: int
    noise_scenario: str
    seed: int
    input_clean_panel: str
    output_observed_panel: str
    factor_file: str
    status: str
    elapsed_seconds: float
    n_rows: int
    n_capped_lower: int
    n_capped_upper: int
    n_capped_total: int
    error_message: str = ""


def default_noise_config() -> NoiseGenerationConfig:
    return parse_noise_config(
        {
            "base_seed": 9900000,
            "sigma_min": 0.0001,
            "price_epsilon": 1e-10,
            "tick_size": 0.01,
            "scenarios": {
                "low_iid": {
                    "enabled": True,
                    "alpha_0": 0.0025,
                    "alpha_m": 1.0,
                    "alpha_tau": 0.05,
                    "tau_min": 0.01984126984126984,
                },
                "spatial_corr": {
                    "enabled": True,
                    "alpha_0": 0.0050,
                    "alpha_m": 1.0,
                    "alpha_tau": 0.05,
                    "tau_min": 0.01984126984126984,
                    "ell_m": 0.10,
                    "ell_tau": 0.25,
                    "correlation_jitter": 1e-10,
                    "max_correlation_jitter": 1e-6,
                },
                "persistent_factor": {
                    "enabled": True,
                    "alpha_0": 0.0015,
                    "alpha_m": 2.0,
                    "alpha_tau": 0.05,
                    "tau_min": 0.01984126984126984,
                    "a_diag": [0.85, 0.65, 0.65],
                    "stationary_factor_std": [0.0007, 0.0025, 0.0008],
                    "q_diag": [1.35975e-7, 3.609375e-6, 3.696e-7],
                    "residual_policy": "match_total_marginal_scale",
                    "factor_initialization": "zero",
                },
            },
        }
    )


def _as_tuple3(value: Any, name: str) -> tuple[float, float, float]:
    try:
        values = tuple(float(x) for x in value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a length-3 sequence") from exc
    if len(values) != PERSISTENT_FACTOR_DIM:
        raise ValueError(f"{name} must have length 3")
    return values  # type: ignore[return-value]


def _as_optional_tuple3(value: Any, name: str) -> tuple[float, float, float] | None:
    if value is None:
        return None
    return _as_tuple3(value, name)


def compute_q_diag_from_stationary_std(a_diag: np.ndarray | Iterable[float], stationary_factor_std: np.ndarray | Iterable[float]) -> np.ndarray:
    a = np.asarray(tuple(a_diag), dtype=float)
    stationary_std = np.asarray(tuple(stationary_factor_std), dtype=float)
    return (1.0 - a * a) * stationary_std * stationary_std


def compute_stationary_std_from_q_diag(a_diag: np.ndarray | Iterable[float], q_diag: np.ndarray | Iterable[float]) -> np.ndarray:
    a = np.asarray(tuple(a_diag), dtype=float)
    q = np.asarray(tuple(q_diag), dtype=float)
    return np.sqrt(q / (1.0 - a * a))


def validate_persistent_factor_config(config: PersistentFactorNoiseConfig) -> PersistentFactorNoiseConfig:
    a_diag = np.asarray(config.a_diag, dtype=float)
    if a_diag.shape != (PERSISTENT_FACTOR_DIM,):
        raise ValueError("persistent_factor.a_diag must have length 3")
    if not np.all(np.isfinite(a_diag)):
        raise ValueError("persistent_factor.a_diag must contain finite values")
    if np.any(np.abs(a_diag) >= 1.0):
        raise ValueError("persistent_factor.a_diag entries must have absolute value strictly below 1")

    stationary_factor_std = None
    if config.stationary_factor_std is not None:
        stationary_factor_std = np.asarray(config.stationary_factor_std, dtype=float)
        if stationary_factor_std.shape != (PERSISTENT_FACTOR_DIM,):
            raise ValueError("persistent_factor.stationary_factor_std must have length 3")
        if not np.all(np.isfinite(stationary_factor_std)):
            raise ValueError("persistent_factor.stationary_factor_std must contain finite values")
        if np.any(stationary_factor_std < 0.0):
            raise ValueError("persistent_factor.stationary_factor_std entries must be non-negative")

    q_diag = None
    if config.q_diag is not None:
        q_diag = np.asarray(config.q_diag, dtype=float)
        if q_diag.shape != (PERSISTENT_FACTOR_DIM,):
            raise ValueError("persistent_factor.q_diag must have length 3")
        if not np.all(np.isfinite(q_diag)):
            raise ValueError("persistent_factor.q_diag must contain finite values")
        if np.any(q_diag < 0.0):
            raise ValueError("persistent_factor.q_diag entries must be non-negative innovation variances")

    if stationary_factor_std is None and q_diag is None:
        raise ValueError("persistent_factor requires stationary_factor_std or q_diag")
    if stationary_factor_std is not None:
        implied_q_diag = compute_q_diag_from_stationary_std(a_diag, stationary_factor_std)
        if q_diag is None:
            q_diag = implied_q_diag
        elif not np.allclose(q_diag, implied_q_diag, rtol=PERSISTENT_FACTOR_Q_RTOL, atol=PERSISTENT_FACTOR_Q_ATOL):
            raise ValueError(
                "persistent_factor.q_diag must be the innovation covariance diagonal implied by "
                "a_diag and stationary_factor_std"
            )
    else:
        stationary_factor_std = compute_stationary_std_from_q_diag(a_diag, q_diag)

    if config.residual_policy not in PERSISTENT_FACTOR_RESIDUAL_POLICIES:
        raise ValueError(
            "persistent_factor.residual_policy must be one of "
            f"{', '.join(PERSISTENT_FACTOR_RESIDUAL_POLICIES)}"
        )
    if config.residual_policy == "legacy_multiplier":
        if config.residual_scale_multiplier is None:
            raise ValueError("persistent_factor.residual_scale_multiplier is required for legacy_multiplier policy")
        if not math.isfinite(config.residual_scale_multiplier) or config.residual_scale_multiplier < 0.0:
            raise ValueError("persistent_factor.residual_scale_multiplier must be non-negative")
    return PersistentFactorNoiseConfig(
        enabled=config.enabled,
        alpha_0=config.alpha_0,
        alpha_m=config.alpha_m,
        alpha_tau=config.alpha_tau,
        tau_min=config.tau_min,
        a_diag=tuple(float(x) for x in a_diag),  # type: ignore[arg-type]
        stationary_factor_std=tuple(float(x) for x in stationary_factor_std),  # type: ignore[arg-type]
        q_diag=tuple(float(x) for x in q_diag),  # type: ignore[arg-type]
        residual_policy=config.residual_policy,
        residual_scale_multiplier=config.residual_scale_multiplier,
        factor_initialization=config.factor_initialization,
    )


def parse_noise_config(raw_noise: dict[str, Any] | None) -> NoiseGenerationConfig:
    if raw_noise is None:
        return default_noise_config()
    scenarios = raw_noise.get("scenarios", {})
    low_raw = scenarios.get("low_iid", {})
    spatial_raw = scenarios.get("spatial_corr", {})
    factor_raw = scenarios.get("persistent_factor", {})
    has_residual_policy = "residual_policy" in factor_raw
    residual_scale_multiplier = factor_raw.get("residual_scale_multiplier")
    residual_policy = str(factor_raw.get("residual_policy", "match_total_marginal_scale" if has_residual_policy or residual_scale_multiplier is None else "legacy_multiplier"))
    persistent_factor = validate_persistent_factor_config(
        PersistentFactorNoiseConfig(
            enabled=bool(factor_raw.get("enabled", True)),
            alpha_0=float(factor_raw.get("alpha_0", 0.0015)),
            alpha_m=float(factor_raw.get("alpha_m", 2.0)),
            alpha_tau=float(factor_raw.get("alpha_tau", 0.05)),
            tau_min=float(factor_raw.get("tau_min", 0.01984126984126984)),
            a_diag=_as_tuple3(factor_raw.get("a_diag", [0.85, 0.65, 0.65]), "persistent_factor.a_diag"),
            stationary_factor_std=_as_optional_tuple3(
                factor_raw.get("stationary_factor_std", [0.0007, 0.0025, 0.0008] if "q_diag" not in factor_raw else None),
                "persistent_factor.stationary_factor_std",
            ),
            q_diag=_as_optional_tuple3(factor_raw.get("q_diag"), "persistent_factor.q_diag"),
            residual_policy=residual_policy,
            residual_scale_multiplier=float(residual_scale_multiplier) if residual_scale_multiplier is not None else None,
            factor_initialization=str(factor_raw.get("factor_initialization", "zero")),
        )
    )
    return NoiseGenerationConfig(
        base_seed=int(raw_noise.get("base_seed", 9900000)),
        sigma_min=float(raw_noise.get("sigma_min", 0.0001)),
        price_epsilon=float(raw_noise.get("price_epsilon", 1e-10)),
        tick_size=float(raw_noise.get("tick_size", 0.01)),
        low_iid=LowIidNoiseConfig(
            enabled=bool(low_raw.get("enabled", True)),
            alpha_0=float(low_raw.get("alpha_0", 0.0025)),
            alpha_m=float(low_raw.get("alpha_m", 1.0)),
            alpha_tau=float(low_raw.get("alpha_tau", 0.05)),
            tau_min=float(low_raw.get("tau_min", 0.01984126984126984)),
        ),
        spatial_corr=SpatialCorrNoiseConfig(
            enabled=bool(spatial_raw.get("enabled", True)),
            alpha_0=float(spatial_raw.get("alpha_0", 0.0050)),
            alpha_m=float(spatial_raw.get("alpha_m", 1.0)),
            alpha_tau=float(spatial_raw.get("alpha_tau", 0.05)),
            tau_min=float(spatial_raw.get("tau_min", 0.01984126984126984)),
            ell_m=float(spatial_raw.get("ell_m", 0.10)),
            ell_tau=float(spatial_raw.get("ell_tau", 0.25)),
            correlation_jitter=float(spatial_raw.get("correlation_jitter", 1e-10)),
            max_correlation_jitter=float(spatial_raw.get("max_correlation_jitter", 1e-6)),
        ),
        persistent_factor=persistent_factor,
    )


def noise_config_to_dict(config: NoiseGenerationConfig) -> dict[str, Any]:
    payload = asdict(config)
    payload["scenarios"] = {
        "low_iid": payload.pop("low_iid"),
        "spatial_corr": payload.pop("spatial_corr"),
        "persistent_factor": payload.pop("persistent_factor"),
    }
    factor_payload = payload["scenarios"]["persistent_factor"]
    if factor_payload.get("residual_scale_multiplier") is None:
        factor_payload.pop("residual_scale_multiplier")
    if factor_payload.get("stationary_factor_std") is None:
        factor_payload.pop("stationary_factor_std")
    if factor_payload.get("q_diag") is None:
        factor_payload.pop("q_diag")
    return payload


def scenario_seed(base_seed: int, sample_id: int, scenario: str) -> int:
    scenario_offsets = {"low_iid": 101, "spatial_corr": 202, "persistent_factor": 303}
    if scenario not in scenario_offsets:
        raise ValueError(f"unknown noise scenario: {scenario}")
    return int(base_seed + 1000 * sample_id + scenario_offsets[scenario])


def read_table(path: str | Path) -> list[dict[str, Any]]:
    p = Path(path)
    if p.suffix == ".parquet":
        import pandas as pd  # type: ignore[import-not-found]

        return pd.read_parquet(p).to_dict("records")
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


def write_table(rows: list[dict[str, Any]], target_without_suffix: str | Path) -> Path:
    target = Path(target_without_suffix)
    target.parent.mkdir(parents=True, exist_ok=True)
    if parquet_available():
        import pandas as pd  # type: ignore[import-not-found]

        out = target.with_suffix(".parquet")
        pd.DataFrame(rows).to_parquet(out, index=False)
        return out
    out = target.with_suffix(".csv")
    fieldnames = list(rows[0].keys()) if rows else []
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return out


def clean_panel_file(run_root: str | Path, sample_id: int) -> Path:
    root = Path(run_root)
    stem = f"sample_{sample_id:03d}"
    for suffix in (".parquet", ".csv"):
        candidate = root / "panels_clean" / f"{stem}{suffix}"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"missing clean panel for sample {sample_id:03d} under {root / 'panels_clean'}")


def observed_panel_file(run_root: str | Path, scenario: str, sample_id: int) -> Path:
    suffix = ".parquet" if parquet_available() else ".csv"
    return Path(run_root) / "panels_observed" / scenario / f"sample_{sample_id:03d}{suffix}"


def _marginal_scale(log_moneyness: np.ndarray, tau: np.ndarray, config: MarginalNoiseScaleConfig) -> np.ndarray:
    return config.alpha_0 * (
        1.0 + config.alpha_m * np.abs(log_moneyness) + config.alpha_tau / np.sqrt(np.maximum(tau, config.tau_min))
    )


def persistent_factor_residual_scale(log_moneyness: np.ndarray, tau: np.ndarray, config: PersistentFactorNoiseConfig) -> np.ndarray:
    marginal_scale = _marginal_scale(log_moneyness, tau, config)
    if config.residual_policy == "match_total_marginal_scale":
        stationary_factor_std = np.array(config.stationary_factor_std, dtype=float)
        factor_variance = (
            stationary_factor_std[0] * stationary_factor_std[0]
            + log_moneyness * log_moneyness * stationary_factor_std[1] * stationary_factor_std[1]
            + tau * tau * stationary_factor_std[2] * stationary_factor_std[2]
        )
        return np.sqrt(np.maximum(marginal_scale * marginal_scale - factor_variance, 0.0))
    if config.residual_policy == "legacy_multiplier":
        return float(config.residual_scale_multiplier) * marginal_scale
    raise ValueError(f"unsupported persistent_factor residual_policy: {config.residual_policy}")


def _price_bounds(S: float, K: float, tau: float, r: float, q: float, option_type: str) -> tuple[float, float]:
    discounted_spot = S * math.exp(-q * tau)
    discounted_strike = K * math.exp(-r * tau)
    kind = option_type.lower()
    if kind == "call":
        return max(discounted_spot - discounted_strike, 0.0), discounted_spot
    if kind == "put":
        return max(discounted_strike - discounted_spot, 0.0), discounted_strike
    raise ValueError("option_type must be 'call' or 'put'")


def _apply_price_mechanics(
    *,
    raw_price: float,
    S: float,
    K: float,
    tau: float,
    r: float,
    q: float,
    option_type: str,
    tick_size: float,
    price_epsilon: float,
) -> tuple[float, float, bool, str]:
    rounded = tick_size * round(raw_price / tick_size)
    lower, upper = _price_bounds(S, K, tau, r, q, option_type)
    lower_cap = lower + price_epsilon
    upper_cap = upper - price_epsilon
    if lower_cap >= upper_cap:
        midpoint = 0.5 * (lower + upper)
        return rounded, midpoint, True, "lower" if rounded <= midpoint else "upper"
    observed = min(upper_cap, max(lower_cap, rounded))
    if observed <= lower_cap and rounded < lower_cap:
        return rounded, observed, True, "lower"
    if observed >= upper_cap and rounded > upper_cap:
        return rounded, observed, True, "upper"
    return rounded, observed, False, "none"


def _contaminate_rows(
    rows: list[dict[str, Any]],
    *,
    scenario: str,
    seed: int,
    raw_iv: np.ndarray,
    noise_draw: np.ndarray,
    config: NoiseGenerationConfig,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for idx, row in enumerate(rows):
        S = float(row["S"])
        K = float(row["strike"])
        tau = float(row["maturity_years"])
        if tau <= 0.0:
            raise ValueError("contaminated panel generation requires strictly positive maturities")
        r = float(row["r"])
        q = float(row["q"])
        option_type = str(row["option_type"]).lower()
        forward = float(row["forward"])
        discount = math.exp(-r * tau)
        raw_price = black76_price(
            forward=forward,
            strike=K,
            tau=tau,
            vol=float(raw_iv[idx]),
            discount_factor=discount,
            option_type=option_type,
        )
        rounded, observed_price, was_capped, cap_direction = _apply_price_mechanics(
            raw_price=raw_price,
            S=S,
            K=K,
            tau=tau,
            r=r,
            q=q,
            option_type=option_type,
            tick_size=config.tick_size,
            price_epsilon=config.price_epsilon,
        )
        observed_iv = implied_vol_black76(
            price=observed_price,
            forward=forward,
            strike=K,
            tau=tau,
            discount_factor=discount,
            option_type=option_type,
            on_bounds="clip",
        )
        item = dict(row)
        item.update(
            {
                "noise_scenario": scenario,
                "raw_contaminated_iv": float(raw_iv[idx]),
                "raw_price_before_rounding": float(raw_price),
                "price_after_rounding": float(rounded),
                "observed_price": float(observed_price),
                "observed_iv": max(float(observed_iv), config.sigma_min),
                "was_price_capped": bool(was_capped),
                "cap_direction": cap_direction,
                "noise_draw": float(noise_draw[idx]),
                "noise_seed": seed,
            }
        )
        out.append(item)
    return out


def _scenario_config(config: NoiseGenerationConfig, scenario: str) -> MarginalNoiseScaleConfig:
    if scenario == "low_iid":
        return config.low_iid
    if scenario == "spatial_corr":
        return config.spatial_corr
    if scenario == "persistent_factor":
        return config.persistent_factor
    raise ValueError(f"unknown noise scenario: {scenario}")


def _low_iid_noise(rows: list[dict[str, Any]], rng: np.random.Generator, config: NoiseGenerationConfig) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    lm = np.array([float(row["log_moneyness"]) for row in rows])
    tau = np.array([float(row["maturity_years"]) for row in rows])
    clean_iv = np.array([float(row["model_iv"]) for row in rows])
    scale = _marginal_scale(lm, tau, config.low_iid)
    shocks = rng.standard_normal(len(rows))
    noise = scale * shocks
    return np.maximum(config.sigma_min, clean_iv + noise), noise, []


def _spatial_corr_noise(rows: list[dict[str, Any]], rng: np.random.Generator, config: NoiseGenerationConfig) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    raw_iv = np.array([float(row["model_iv"]) for row in rows])
    noise = np.zeros(len(rows), dtype=float)
    week_values = sorted({int(row["week_index"]) for row in rows})
    for week in week_values:
        idx = np.array([i for i, row in enumerate(rows) if int(row["week_index"]) == week], dtype=int)
        lm = np.array([float(rows[i]["log_moneyness"]) for i in idx])
        tau = np.array([float(rows[i]["maturity_years"]) for i in idx])
        scale = _marginal_scale(lm, tau, config.spatial_corr)
        distance = (
            np.abs(lm[:, None] - lm[None, :]) / config.spatial_corr.ell_m
            + np.abs(tau[:, None] - tau[None, :]) / config.spatial_corr.ell_tau
        )
        corr = np.exp(-distance)
        jitter = config.spatial_corr.correlation_jitter
        while True:
            try:
                chol = np.linalg.cholesky(corr + jitter * np.eye(len(idx)))
                break
            except np.linalg.LinAlgError:
                jitter *= 10.0
                if jitter > config.spatial_corr.max_correlation_jitter:
                    eigvals, eigvecs = np.linalg.eigh(corr)
                    chol = eigvecs @ np.diag(np.sqrt(np.maximum(eigvals, 0.0)))
                    break
        noise[idx] = scale * (chol @ rng.standard_normal(len(idx)))
    return np.maximum(config.sigma_min, raw_iv + noise), noise, []


def _persistent_factor_noise(rows: list[dict[str, Any]], rng: np.random.Generator, config: NoiseGenerationConfig) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    clean_iv = np.array([float(row["model_iv"]) for row in rows])
    noise = np.zeros(len(rows), dtype=float)
    factors: list[dict[str, Any]] = []
    factor_config = validate_persistent_factor_config(config.persistent_factor)
    if factor_config.factor_initialization != "zero":
        raise NotImplementedError("only zero factor initialization is implemented")
    a_diag = np.array(factor_config.a_diag, dtype=float)
    q_diag = np.array(factor_config.q_diag, dtype=float)
    f = np.zeros(3, dtype=float)
    week_values = sorted({int(row["week_index"]) for row in rows})
    for week in week_values:
        f = a_diag * f + rng.normal(loc=0.0, scale=np.sqrt(q_diag), size=3)
        factors.append({"week_index": week, "factor_0": f[0], "factor_1": f[1], "factor_2": f[2]})
        idx = np.array([i for i, row in enumerate(rows) if int(row["week_index"]) == week], dtype=int)
        lm = np.array([float(rows[i]["log_moneyness"]) for i in idx])
        tau = np.array([float(rows[i]["maturity_years"]) for i in idx])
        scale = persistent_factor_residual_scale(lm, tau, factor_config)
        factor_component = f[0] + lm * f[1] + tau * f[2]
        residual = scale * rng.standard_normal(len(idx))
        noise[idx] = factor_component + residual
    return np.maximum(config.sigma_min, clean_iv + noise), noise, factors


def generate_noisy_panel_rows(
    clean_rows: list[dict[str, Any]],
    *,
    scenario: str,
    seed: int,
    config: NoiseGenerationConfig,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if scenario not in NOISE_SCENARIOS:
        raise ValueError(f"unsupported noise scenario: {scenario}")
    rows = sorted(clean_rows, key=lambda row: (int(row["week_index"]), float(row["maturity_years"]), float(row["log_moneyness"])))
    rng = np.random.default_rng(seed)
    if scenario == "low_iid":
        raw_iv, noise, factors = _low_iid_noise(rows, rng, config)
    elif scenario == "spatial_corr":
        raw_iv, noise, factors = _spatial_corr_noise(rows, rng, config)
    else:
        raw_iv, noise, factors = _persistent_factor_noise(rows, rng, config)
    return _contaminate_rows(rows, scenario=scenario, seed=seed, raw_iv=raw_iv, noise_draw=noise, config=config), factors


def write_factor_file(factors: list[dict[str, Any]], run_root: str | Path, scenario: str, sample_id: int) -> Path:
    target = Path(run_root) / "noise_factors" / scenario / f"sample_{sample_id:03d}.csv"
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["week_index", "factor_0", "factor_1", "factor_2"])
        writer.writeheader()
        writer.writerows(factors)
    return target


def generate_noisy_panel_file(
    *,
    clean_panel_path: str | Path,
    run_root: str | Path,
    sample_id: int,
    scenario: str,
    config: NoiseGenerationConfig,
    skip_existing: bool = False,
) -> NoisyPanelResult:
    import time
    import traceback

    started = time.perf_counter()
    seed = scenario_seed(config.base_seed, sample_id, scenario)
    output = observed_panel_file(run_root, scenario, sample_id)
    factor_file = Path(run_root) / "noise_factors" / scenario / f"sample_{sample_id:03d}.csv"
    if skip_existing and output.exists() and (scenario != "persistent_factor" or factor_file.exists()):
        rows = read_table(output)
        lower = sum(1 for row in rows if str(row.get("cap_direction")) == "lower")
        upper = sum(1 for row in rows if str(row.get("cap_direction")) == "upper")
        return NoisyPanelResult(
            sample_id=sample_id,
            noise_scenario=scenario,
            seed=seed,
            input_clean_panel=str(clean_panel_path),
            output_observed_panel=str(output),
            factor_file=str(factor_file if scenario == "persistent_factor" else ""),
            status="skipped",
            elapsed_seconds=time.perf_counter() - started,
            n_rows=len(rows),
            n_capped_lower=lower,
            n_capped_upper=upper,
            n_capped_total=lower + upper,
        )
    try:
        clean_rows = read_table(clean_panel_path)
        rows, factors = generate_noisy_panel_rows(clean_rows, scenario=scenario, seed=seed, config=config)
        output = write_table(rows, Path(run_root) / "panels_observed" / scenario / f"sample_{sample_id:03d}")
        factor_out = ""
        if scenario == "persistent_factor":
            factor_out = str(write_factor_file(factors, run_root, scenario, sample_id))
        lower = sum(1 for row in rows if row["cap_direction"] == "lower")
        upper = sum(1 for row in rows if row["cap_direction"] == "upper")
        return NoisyPanelResult(
            sample_id=sample_id,
            noise_scenario=scenario,
            seed=seed,
            input_clean_panel=str(clean_panel_path),
            output_observed_panel=str(output),
            factor_file=factor_out,
            status="ok",
            elapsed_seconds=time.perf_counter() - started,
            n_rows=len(rows),
            n_capped_lower=lower,
            n_capped_upper=upper,
            n_capped_total=lower + upper,
        )
    except Exception as exc:
        return NoisyPanelResult(
            sample_id=sample_id,
            noise_scenario=scenario,
            seed=seed,
            input_clean_panel=str(clean_panel_path),
            output_observed_panel=str(output),
            factor_file=str(factor_file if scenario == "persistent_factor" else ""),
            status="error",
            elapsed_seconds=time.perf_counter() - started,
            n_rows=0,
            n_capped_lower=0,
            n_capped_upper=0,
            n_capped_total=0,
            error_message=f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
        )


def write_noise_config(run_root: str | Path, config: NoiseGenerationConfig) -> None:
    config_dir = Path(run_root) / "config"
    config_dir.mkdir(parents=True, exist_ok=True)
    with (config_dir / "noise_scenarios.json").open("w") as fh:
        json.dump(noise_config_to_dict(config), fh, indent=2)


def write_noisy_manifest(run_root: str | Path, results: Iterable[NoisyPanelResult]) -> None:
    config_dir = Path(run_root) / "config"
    config_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "noise_scenario",
        "seed",
        "input_clean_panel",
        "output_observed_panel",
        "factor_file",
        "status",
        "elapsed_seconds",
        "n_rows",
        "n_capped_lower",
        "n_capped_upper",
        "n_capped_total",
        "error_message",
    ]
    with (config_dir / "manifest_noisy_panels.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for result in sorted(results, key=lambda item: (item.sample_id, item.noise_scenario)):
            writer.writerow(asdict(result))
