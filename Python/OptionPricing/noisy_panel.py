from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_price
from OptionData.io import panel_metadata_path, write_panel_metadata
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
class NoiseSettings:
    """Common noise controls and the three explicit scenario mappings."""

    base_seed: int
    sigma_min: float
    price_epsilon: float
    tick_size: float
    scenarios: dict[str, dict[str, Any]]

    def scenario_names(self) -> tuple[str, ...]:
        return tuple(name for name in NOISE_SCENARIOS if name in self.scenarios)


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


def default_noise_settings() -> NoiseSettings:
    """Small deterministic fixture used by unit tests and research scripts."""

    settings = NoiseSettings(
        base_seed=9900000,
        sigma_min=0.0001,
        price_epsilon=1e-10,
        tick_size=0.01,
        scenarios={
            "low_iid": {
                "alpha_0": 0.0025,
                "alpha_m": 1.0,
                "alpha_tau": 0.05,
                "tau_min": 0.01984126984126984,
            },
            "spatial_corr": {
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
    )
    validate_noise_settings(settings)
    return settings


def compute_q_diag_from_stationary_std(a_diag: np.ndarray | Iterable[float], stationary_factor_std: np.ndarray | Iterable[float]) -> np.ndarray:
    a = np.asarray(tuple(a_diag), dtype=float)
    stationary_std = np.asarray(tuple(stationary_factor_std), dtype=float)
    return (1.0 - a * a) * stationary_std * stationary_std


def compute_stationary_std_from_q_diag(a_diag: np.ndarray | Iterable[float], q_diag: np.ndarray | Iterable[float]) -> np.ndarray:
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
        raise ValueError("persistent_factor.a_diag entries must have absolute value strictly below 1")

    stationary_factor_std = None
    if config.get("stationary_factor_std") is not None:
        stationary_factor_std = np.asarray(config["stationary_factor_std"], dtype=float)
        if stationary_factor_std.shape != (PERSISTENT_FACTOR_DIM,):
            raise ValueError("persistent_factor.stationary_factor_std must have length 3")
        if not np.all(np.isfinite(stationary_factor_std)):
            raise ValueError("persistent_factor.stationary_factor_std must contain finite values")
        if np.any(stationary_factor_std < 0.0):
            raise ValueError("persistent_factor.stationary_factor_std entries must be non-negative")

    q_diag = None
    if config.get("q_diag") is not None:
        q_diag = np.asarray(config["q_diag"], dtype=float)
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

    residual_policy = str(config.get("residual_policy", "match_total_marginal_scale"))
    if residual_policy not in PERSISTENT_FACTOR_RESIDUAL_POLICIES:
        raise ValueError(
            "persistent_factor.residual_policy must be one of "
            f"{', '.join(PERSISTENT_FACTOR_RESIDUAL_POLICIES)}"
        )
    if residual_policy == "legacy_multiplier":
        residual_scale_multiplier = config.get("residual_scale_multiplier")
        if residual_scale_multiplier is None:
            raise ValueError("persistent_factor.residual_scale_multiplier is required for legacy_multiplier policy")
        if not math.isfinite(float(residual_scale_multiplier)) or float(residual_scale_multiplier) < 0.0:
            raise ValueError("persistent_factor.residual_scale_multiplier must be non-negative")


def validate_noise_settings(config: NoiseSettings) -> None:
    if config.sigma_min <= 0.0 or config.price_epsilon <= 0.0 or config.tick_size <= 0.0:
        raise ValueError("noise sigma_min, price_epsilon and tick_size must be positive")
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
        validate_persistent_factor_settings(config.scenarios["persistent_factor"])


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


def write_table(
    rows: list[dict[str, Any]],
    target_without_suffix: str | Path,
    *,
    metadata: dict[str, Any] | None = None,
    panel_format: str,
) -> Path:
    target = Path(target_without_suffix)
    target.parent.mkdir(parents=True, exist_ok=True)
    if panel_format == "parquet":
        if not parquet_available():
            raise RuntimeError("panel_format='parquet' requires pandas and pyarrow or fastparquet")
        import pandas as pd  # type: ignore[import-not-found]

        out = target.with_suffix(".parquet")
        temporary = out.with_name(out.stem + ".tmp.parquet")
        pd.DataFrame(rows).to_parquet(temporary, index=False)
        if len(pd.read_parquet(temporary)) != len(rows):
            temporary.unlink(missing_ok=True)
            raise RuntimeError("atomic noisy Parquet validation failed before publication")
        os.replace(temporary, out)
        if metadata is not None:
            write_panel_metadata(out, metadata)
        return out
    if panel_format != "csv":
        raise ValueError("panel_format must be 'parquet' or 'csv'")
    out = target.with_suffix(".csv")
    temporary = out.with_name(out.stem + ".tmp.csv")
    fieldnames = list(rows[0].keys()) if rows else []
    with temporary.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(temporary, out)
    if metadata is not None:
        write_panel_metadata(out, metadata)
    return out


def clean_panel_file(run_root: str | Path, sample_id: int) -> Path:
    root = Path(run_root)
    stem = f"sample_{sample_id:03d}"
    for suffix in (".parquet", ".csv"):
        candidate = root / "panels_clean" / f"{stem}{suffix}"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"missing clean panel for sample {sample_id:03d} under {root / 'panels_clean'}")


def observed_panel_file(run_root: str | Path, scenario: str, sample_id: int, panel_format: str) -> Path:
    if panel_format not in {"parquet", "csv"}:
        raise ValueError("panel_format must be 'parquet' or 'csv'")
    suffix = f".{panel_format}"
    return Path(run_root) / "panels_observed" / scenario / f"sample_{sample_id:03d}{suffix}"


def _marginal_scale(log_moneyness: np.ndarray, tau: np.ndarray, config: dict[str, Any]) -> np.ndarray:
    return float(config["alpha_0"]) * (
        1.0
        + float(config["alpha_m"]) * np.abs(log_moneyness)
        + float(config["alpha_tau"]) / np.sqrt(np.maximum(tau, float(config["tau_min"])))
    )


def persistent_factor_residual_scale(
    log_moneyness: np.ndarray,
    tau: np.ndarray,
    config: dict[str, Any],
) -> np.ndarray:
    marginal_scale = _marginal_scale(log_moneyness, tau, config)
    if config["residual_policy"] == "match_total_marginal_scale":
        stationary_factor_std = np.array(config["stationary_factor_std"], dtype=float)
        factor_variance = (
            stationary_factor_std[0] * stationary_factor_std[0]
            + log_moneyness * log_moneyness * stationary_factor_std[1] * stationary_factor_std[1]
            + tau * tau * stationary_factor_std[2] * stationary_factor_std[2]
        )
        return np.sqrt(np.maximum(marginal_scale * marginal_scale - factor_variance, 0.0))
    if config["residual_policy"] == "legacy_multiplier":
        return float(config["residual_scale_multiplier"]) * marginal_scale
    raise ValueError(f"unsupported persistent_factor residual_policy: {config['residual_policy']}")


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
    config: NoiseSettings,
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


def _low_iid_noise(rows: list[dict[str, Any]], rng: np.random.Generator, config: NoiseSettings) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    lm = np.array([float(row["log_moneyness"]) for row in rows])
    tau = np.array([float(row["maturity_years"]) for row in rows])
    clean_iv = np.array([float(row["model_iv"]) for row in rows])
    scale = _marginal_scale(lm, tau, config.scenarios["low_iid"])
    shocks = rng.standard_normal(len(rows))
    noise = scale * shocks
    return np.maximum(config.sigma_min, clean_iv + noise), noise, []


def _spatial_corr_noise(rows: list[dict[str, Any]], rng: np.random.Generator, config: NoiseSettings) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    raw_iv = np.array([float(row["model_iv"]) for row in rows])
    noise = np.zeros(len(rows), dtype=float)
    week_values = sorted({int(row["week_index"]) for row in rows})
    for week in week_values:
        idx = np.array([i for i, row in enumerate(rows) if int(row["week_index"]) == week], dtype=int)
        lm = np.array([float(rows[i]["log_moneyness"]) for i in idx])
        tau = np.array([float(rows[i]["maturity_years"]) for i in idx])
        spatial = config.scenarios["spatial_corr"]
        scale = _marginal_scale(lm, tau, spatial)
        distance = (
            np.abs(lm[:, None] - lm[None, :]) / float(spatial["ell_m"])
            + np.abs(tau[:, None] - tau[None, :]) / float(spatial["ell_tau"])
        )
        corr = np.exp(-distance)
        jitter = float(spatial["correlation_jitter"])
        while True:
            try:
                chol = np.linalg.cholesky(corr + jitter * np.eye(len(idx)))
                break
            except np.linalg.LinAlgError:
                jitter *= 10.0
                if jitter > float(spatial["max_correlation_jitter"]):
                    eigvals, eigvecs = np.linalg.eigh(corr)
                    chol = eigvecs @ np.diag(np.sqrt(np.maximum(eigvals, 0.0)))
                    break
        noise[idx] = scale * (chol @ rng.standard_normal(len(idx)))
    return np.maximum(config.sigma_min, raw_iv + noise), noise, []


def _persistent_factor_noise(rows: list[dict[str, Any]], rng: np.random.Generator, config: NoiseSettings) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    clean_iv = np.array([float(row["model_iv"]) for row in rows])
    noise = np.zeros(len(rows), dtype=float)
    factors: list[dict[str, Any]] = []
    factor_config = config.scenarios["persistent_factor"]
    validate_persistent_factor_settings(factor_config)
    if factor_config["factor_initialization"] != "zero":
        raise NotImplementedError("only zero factor initialization is implemented")
    a_diag = np.array(factor_config["a_diag"], dtype=float)
    q_diag = np.array(factor_config["q_diag"], dtype=float)
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
    config: NoiseSettings,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if scenario not in NOISE_SCENARIOS:
        raise ValueError(f"unsupported noise scenario: {scenario}")
    if scenario not in config.scenarios:
        raise ValueError(f"noise scenario is not configured: {scenario}")
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
    temporary = target.with_name(target.stem + ".tmp.csv")
    with temporary.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["week_index", "factor_0", "factor_1", "factor_2"])
        writer.writeheader()
        writer.writerows(factors)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(temporary, target)
    return target


def generate_noisy_panel_file(
    *,
    clean_panel_path: str | Path,
    run_root: str | Path,
    sample_id: int,
    scenario: str,
    config: NoiseSettings,
    panel_format: str,
    skip_existing: bool = False,
) -> NoisyPanelResult:
    import time
    import traceback

    started = time.perf_counter()
    seed = scenario_seed(config.base_seed, sample_id, scenario)
    output = observed_panel_file(run_root, scenario, sample_id, panel_format)
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
        clean_metadata_path = panel_metadata_path(clean_panel_path)
        if not clean_metadata_path.exists():
            raise ValueError("clean panel is missing required COS-basis metadata sidecar")
        with clean_metadata_path.open() as fh:
            inherited_metadata = json.load(fh)
        lower = sum(1 for row in rows if row["cap_direction"] == "lower")
        upper = sum(1 for row in rows if row["cap_direction"] == "upper")
        inherited_metadata.update(
            {
                "scenario": scenario,
                "sample_id": sample_id,
                "noise_seed": seed,
                "n_capped_lower": lower,
                "n_capped_upper": upper,
                "n_capped_total": lower + upper,
            }
        )
        output = write_table(
            rows,
            Path(run_root) / "panels_observed" / scenario / f"sample_{sample_id:03d}",
            metadata=inherited_metadata,
            panel_format=panel_format,
        )
        factor_out = ""
        if scenario == "persistent_factor":
            factor_out = str(write_factor_file(factors, run_root, scenario, sample_id))
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
    target = config_dir / "manifest_noisy_panels.csv"
    temporary = target.with_name(target.stem + ".tmp.csv")
    with temporary.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for result in sorted(results, key=lambda item: (item.sample_id, item.noise_scenario)):
            writer.writerow(asdict(result))
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(temporary, target)
