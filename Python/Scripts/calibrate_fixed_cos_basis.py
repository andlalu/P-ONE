from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_vega
from Models.Heston.parameters import HestonRiskNeutralParameters
from OptionPricing.cos_basis import FixedCosBasisConfig
from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.heston_fourier_reference import heston_option_price_fourier_reference

LOGGER = logging.getLogger(__name__)
MATURITIES = (1.0 / 12.0, 0.25, 0.5)
LOG_MONEYNESS = (-0.15, -0.075, 0.0, 0.075, 0.15)
VARIANCE_STATES = (1e-6, 1e-4, 0.0025, 0.01, 0.0225, 0.04, 0.08, 0.15)
CANDIDATE_WIDTHS = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
CANDIDATE_TERMS = (64, 128, 256, 512, 576, 640, 704, 768, 1024, 1536, 2048, 3072, 4096)
SPOT = 100.0
RATE = 0.02
DIVIDEND_YIELD = 0.0
DESIGN_SEED = 21072026
REFERENCE_WIDTH = 5.0
REFERENCE_TERMS = 8192
CALIBRATION_MIN_VEGA = 5e-6


@dataclass(frozen=True)
class CalibrationParameterPoint:
    point_id: str
    kappa_q: float
    vbar_q: float
    sigma_v: float
    rho: float
    source: str

    def risk_neutral_parameters(self) -> HestonRiskNeutralParameters:
        return HestonRiskNeutralParameters(
            kappa=self.kappa_q,
            vbar=self.vbar_q,
            sigma_v=self.sigma_v,
            rho=self.rho,
            r=RATE,
            q=DIVIDEND_YIELD,
        )


def calibration_parameter_design() -> list[CalibrationParameterPoint]:
    """Return a deterministic LHS over the production-induced Q-domain plus anchors."""

    from scipy.stats import qmc  # type: ignore[import-not-found]

    lower = np.array([1.5, 0.012, 0.20, -0.80])
    upper = np.array([4.0, 0.24, 0.65, -0.20])
    sample = qmc.scale(qmc.LatinHypercube(4, seed=DESIGN_SEED).random(20), lower, upper)
    anchors = [
        ("true_dgp", 2.0, 0.07875, 0.40, -0.50, "true DGP"),
        ("poc_start", 2.0, 0.0816, 0.42, -0.52, "POC base start"),
        ("moderate_low", 1.75, 0.035, 0.30, -0.65, "moderate perturbation"),
        ("moderate_high", 3.25, 0.14, 0.55, -0.35, "moderate perturbation"),
    ]
    points = [
        CalibrationParameterPoint(
            point_id=f"lhs_{index:02d}",
            kappa_q=float(row[0]),
            vbar_q=float(row[1]),
            sigma_v=float(row[2]),
            rho=float(row[3]),
            source="seeded Latin hypercube",
        )
        for index, row in enumerate(sample)
    ]
    points.extend(CalibrationParameterPoint(*anchor) for anchor in anchors)
    for point in points:
        point.risk_neutral_parameters().validate()
    return points


def _contract(maturity: float, log_moneyness: float) -> tuple[float, float, str, float]:
    forward = SPOT * math.exp((RATE - DIVIDEND_YIELD) * maturity)
    strike = forward * math.exp(log_moneyness)
    option_type = "put" if log_moneyness < 0.0 else "call"
    return forward, strike, option_type, math.exp(-RATE * maturity)


def _implied_volatility(price: float, maturity: float, log_moneyness: float) -> float:
    forward, strike, option_type, discount = _contract(maturity, log_moneyness)
    return implied_vol_black76(
        price=price,
        forward=forward,
        strike=strike,
        tau=maturity,
        discount_factor=discount,
        option_type=option_type,
        on_bounds="raise",
    )


def _cos_prices(
    pricer: CosOptionPricer,
    parameters: HestonRiskNeutralParameters,
    maturity: float,
    variance_states: np.ndarray,
    effective_width: float,
    n_cos: int,
) -> np.ndarray:
    strikes = np.array([_contract(maturity, value)[1] for value in LOG_MONEYNESS])
    option_types = np.array(
        [[[_contract(maturity, value)[2]] for value in LOG_MONEYNESS] for _ in variance_states]
    ).reshape(variance_states.size, len(LOG_MONEYNESS), 1)
    basis = pricer.prepare_fixed_basis(
        maturity=maturity,
        effective_width=effective_width,
        n_cos=n_cos,
        model_params=parameters,
    )
    return pricer.price_matrix_fixed_basis(
        log_s=np.full(variance_states.size, math.log(SPOT)),
        variance=variance_states,
        strike_grid=strikes,
        rate=RATE,
        dividend_yield=DIVIDEND_YIELD,
        basis=basis,
        option_type=option_types,
    )[:, :, 0]


def _reference_records(points: list[CalibrationParameterPoint]) -> list[dict[str, object]]:
    pricer = CosOptionPricer()
    states = np.asarray(VARIANCE_STATES)
    records: list[dict[str, object]] = []
    for point in points:
        parameters = point.risk_neutral_parameters()
        for maturity in MATURITIES:
            prices = _cos_prices(pricer, parameters, maturity, states, REFERENCE_WIDTH, REFERENCE_TERMS)
            for variance_index, variance in enumerate(states):
                for money_index, log_moneyness in enumerate(LOG_MONEYNESS):
                    price = float(prices[variance_index, money_index])
                    implied_volatility = _implied_volatility(price, maturity, log_moneyness)
                    _, strike, option_type, discount = _contract(maturity, log_moneyness)
                    vega = black76_vega(
                        forward=SPOT * math.exp(RATE * maturity),
                        strike=strike,
                        tau=maturity,
                        vol=implied_volatility,
                        discount_factor=discount,
                    )
                    records.append(
                        {
                            "point_id": point.point_id,
                            "maturity": maturity,
                            "variance": float(variance),
                            "log_moneyness": log_moneyness,
                            "option_type": option_type,
                            "reference_price": price,
                            "reference_iv": implied_volatility,
                            "reference_vega": vega,
                        }
                    )
    return records


def _reference_crosscheck(
    points: list[CalibrationParameterPoint], reference_records: list[dict[str, object]]
) -> list[dict[str, object]]:
    record_lookup = {
        (record["point_id"], record["maturity"], record["variance"], record["log_moneyness"]): record
        for record in reference_records
    }
    selected_points = [points[0], points[7], points[-4], points[-1]]
    selected_states = (1e-6, 0.01, 0.15)
    selected_money = (-0.15, 0.0, 0.15)
    rows: list[dict[str, object]] = []
    for point in selected_points:
        parameters = point.risk_neutral_parameters()
        for maturity in MATURITIES:
            for variance in selected_states:
                for log_moneyness in selected_money:
                    _, strike, option_type, _ = _contract(maturity, log_moneyness)
                    started = time.perf_counter()
                    fourier_price = heston_option_price_fourier_reference(
                        spot=SPOT,
                        initial_variance=variance,
                        maturity=maturity,
                        strike=strike,
                        option_type=option_type,
                        parameters=parameters,
                        integration_limit=1500.0,
                        absolute_tolerance=2e-10,
                        relative_tolerance=2e-10,
                        integration_subinterval_limit=600,
                    )
                    reference = record_lookup[(point.point_id, maturity, variance, log_moneyness)]
                    rows.append(
                        {
                            "point_id": point.point_id,
                            "maturity": maturity,
                            "variance": variance,
                            "log_moneyness": log_moneyness,
                            "fourier_price": fourier_price,
                            "high_resolution_cos_price": reference["reference_price"],
                            "absolute_price_difference": abs(fourier_price - float(reference["reference_price"])),
                            "fourier_runtime_seconds": time.perf_counter() - started,
                            "integration_limit": 1500.0,
                            "absolute_tolerance": 2e-10,
                            "relative_tolerance": 2e-10,
                            "high_resolution_width": REFERENCE_WIDTH,
                            "high_resolution_n_cos": REFERENCE_TERMS,
                        }
                    )
    return rows


def _candidate_summary(
    points: list[CalibrationParameterPoint], reference_records: list[dict[str, object]]
) -> list[dict[str, object]]:
    pricer = CosOptionPricer()
    states = np.asarray(VARIANCE_STATES)
    summaries: list[dict[str, object]] = []
    for maturity in MATURITIES:
        maturity_records = [record for record in reference_records if record["maturity"] == maturity]
        for n_cos in CANDIDATE_TERMS:
            for effective_width in CANDIDATE_WIDTHS:
                started = time.perf_counter()
                price_errors: list[float] = []
                relative_errors: list[float] = []
                iv_errors: list[float] = []
                regular_vega_iv_errors: list[float] = []
                vega_errors: list[float] = []
                invalid_prices = 0
                invalid_ivs = 0
                record_offset = 0
                for point in points:
                    candidate_prices = _cos_prices(
                        pricer,
                        point.risk_neutral_parameters(),
                        maturity,
                        states,
                        effective_width,
                        n_cos,
                    )
                    for variance_index, _ in enumerate(states):
                        for money_index, log_moneyness in enumerate(LOG_MONEYNESS):
                            reference = maturity_records[record_offset]
                            record_offset += 1
                            candidate_price = float(candidate_prices[variance_index, money_index])
                            reference_price = float(reference["reference_price"])
                            price_error = abs(candidate_price - reference_price)
                            price_errors.append(price_error)
                            if reference_price > 1e-8:
                                relative_errors.append(price_error / reference_price)
                            reference_vega = float(reference["reference_vega"])
                            vega_errors.append(price_error / max(reference_vega, 1e-8))
                            forward, strike, option_type, discount = _contract(maturity, log_moneyness)
                            lower = discount * max(
                                forward - strike if option_type == "call" else strike - forward, 0.0
                            )
                            upper = discount * (forward if option_type == "call" else strike)
                            if not math.isfinite(candidate_price) or not lower - 1e-10 <= candidate_price <= upper + 1e-10:
                                invalid_prices += 1
                                invalid_ivs += 1
                                continue
                            try:
                                candidate_iv = _implied_volatility(candidate_price, maturity, log_moneyness)
                            except (ValueError, RuntimeError, OverflowError):
                                invalid_ivs += 1
                                continue
                            if not math.isfinite(candidate_iv):
                                invalid_ivs += 1
                            else:
                                iv_error = abs(candidate_iv - float(reference["reference_iv"]))
                                iv_errors.append(iv_error)
                                if reference_vega >= CALIBRATION_MIN_VEGA:
                                    regular_vega_iv_errors.append(iv_error)
                iv_array = np.asarray(iv_errors or [math.inf])
                regular_vega_iv_array = np.asarray(regular_vega_iv_errors or [math.inf])
                price_array = np.asarray(price_errors or [math.inf])
                relative_array = np.asarray(relative_errors or [math.inf])
                vega_array = np.asarray(vega_errors or [math.inf])
                summaries.append(
                    {
                        "maturity": maturity,
                        "effective_width": effective_width,
                        "n_cos": n_cos,
                        "observation_count": len(maturity_records),
                        "invalid_price_count": invalid_prices,
                        "invalid_iv_count": invalid_ivs,
                        "max_absolute_price_error": float(np.max(price_array)),
                        "p99_absolute_price_error": float(np.quantile(price_array, 0.99)),
                        "median_absolute_price_error": float(np.median(price_array)),
                        "max_relative_price_error": float(np.max(relative_array)),
                        "p99_relative_price_error": float(np.quantile(relative_array, 0.99)),
                        "median_relative_price_error": float(np.median(relative_array)),
                        "max_absolute_iv_error": float(np.max(iv_array)),
                        "p99_absolute_iv_error": float(np.quantile(iv_array, 0.99)),
                        "median_absolute_iv_error": float(np.median(iv_array)),
                        "low_vega_count": len(iv_errors) - len(regular_vega_iv_errors),
                        "max_absolute_iv_error_above_vega_threshold": float(np.max(regular_vega_iv_array)),
                        "p99_absolute_iv_error_above_vega_threshold": float(
                            np.quantile(regular_vega_iv_array, 0.99)
                        ),
                        "max_vega_scaled_price_error": float(np.max(vega_array)),
                        "p99_vega_scaled_price_error": float(np.quantile(vega_array, 0.99)),
                        "median_vega_scaled_price_error": float(np.median(vega_array)),
                        "runtime_seconds": time.perf_counter() - started,
                    }
                )
                LOGGER.info(
                    "calibrated maturity %.8f, width %.2f, n_cos %d: max IV error %.3g",
                    maturity,
                    effective_width,
                    n_cos,
                    summaries[-1]["max_absolute_iv_error"],
                )
    return summaries


def _select_configuration(summaries: list[dict[str, object]]) -> dict[str, object]:
    def select_for_tolerance(tolerance: float) -> tuple[int, list[dict[str, object]], float]:
        alternatives: list[tuple[float, int, list[dict[str, object]]]] = []
        for n_cos in CANDIDATE_TERMS:
            selected_rows: list[dict[str, object]] = []
            for maturity in MATURITIES:
                eligible = [
                    row
                    for row in summaries
                    if row["maturity"] == maturity
                    and row["n_cos"] == n_cos
                    and row["invalid_price_count"] == 0
                    and row["invalid_iv_count"] == 0
                    and float(row["max_absolute_iv_error_above_vega_threshold"]) <= tolerance
                    and float(row["max_absolute_iv_error"]) <= 2.5e-4
                ]
                if not eligible:
                    break
                selected_rows.append(min(eligible, key=lambda row: float(row["effective_width"])))
            if len(selected_rows) == len(MATURITIES):
                runtime = sum(float(row["runtime_seconds"]) for row in selected_rows)
                alternatives.append((runtime, n_cos, selected_rows))
        if not alternatives:
            raise RuntimeError(f"no candidate configuration meets IV tolerance {tolerance}")
        runtime, n_cos, rows = min(alternatives, key=lambda item: (item[1], item[0]))
        return n_cos, rows, runtime

    generation_n_cos, generation_rows, generation_runtime = select_for_tolerance(2e-6)
    estimation_n_cos, estimation_rows, estimation_runtime = select_for_tolerance(1e-5)
    widths = [float(row["effective_width"]) for row in generation_rows]
    for maturity, width in zip(MATURITIES, widths):
        estimation_row = next(row for row in estimation_rows if row["maturity"] == maturity)
        compatible = next(
            (
                row
                for row in summaries
                if row["maturity"] == maturity
                and row["n_cos"] == estimation_n_cos
                and row["effective_width"] == width
            ),
            None,
        )
        if (
            compatible is None
            or float(compatible["max_absolute_iv_error_above_vega_threshold"]) > 1e-5
            or float(compatible["max_absolute_iv_error"]) > 2.5e-4
        ):
            raise RuntimeError("lowest-cost estimation choice is incompatible with generation widths")
    compatible_estimation_rows = [
        next(
            row
            for row in summaries
            if row["maturity"] == maturity
            and row["n_cos"] == estimation_n_cos
            and row["effective_width"] == width
        )
        for maturity, width in zip(MATURITIES, widths)
    ]
    return {
        "selection_rule": "lowest term count, then narrowest interval, meeting all zero-failure IV-error constraints",
        "maturities": list(MATURITIES),
        "effective_widths": widths,
        "generation_n_cos": generation_n_cos,
        "estimation_n_cos": estimation_n_cos,
        "generation_max_absolute_iv_error": max(float(row["max_absolute_iv_error"]) for row in generation_rows),
        "generation_p99_absolute_iv_error": max(float(row["p99_absolute_iv_error"]) for row in generation_rows),
        "estimation_max_absolute_iv_error": max(
            float(row["max_absolute_iv_error"]) for row in compatible_estimation_rows
        ),
        "estimation_p99_absolute_iv_error": max(
            float(row["p99_absolute_iv_error"]) for row in compatible_estimation_rows
        ),
        "generation_max_absolute_iv_error_above_vega_threshold": max(
            float(row["max_absolute_iv_error_above_vega_threshold"]) for row in generation_rows
        ),
        "estimation_max_absolute_iv_error_above_vega_threshold": max(
            float(row["max_absolute_iv_error_above_vega_threshold"])
            for row in compatible_estimation_rows
        ),
        "low_vega_exception": {
            "black_vega_threshold": CALIBRATION_MIN_VEGA,
            "maximum_absolute_iv_error_allowed": 2.5e-4,
            "justification": "At option prices near 1e-10, price differences near floating-point and quadrature accuracy map to unstable IV differences. The full-domain cap remains six times smaller than the 0.0015 observation-error scale.",
        },
        "generation_candidate_runtime_seconds": generation_runtime,
        "estimation_candidate_runtime_seconds": sum(
            float(row["runtime_seconds"]) for row in compatible_estimation_rows
        ),
        "accepted_variance_state_range": [min(VARIANCE_STATES), max(VARIANCE_STATES)],
        "accepted_risk_neutral_parameter_domain": {
            "kappa_q": [1.5, 4.0],
            "vbar_q": [0.012, 0.24],
            "sigma_v": [0.20, 0.65],
            "rho": [-0.80, -0.20],
        },
        "production_optimizer_natural_bounds": {
            "eta": [2.0, 8.0],
            "kappa": [4.0, 9.0],
            "vbar": [0.012, 0.04],
            "sigma_v": [0.20, 0.65],
            "rho": [-0.80, -0.20],
            "kappa_q": [1.5, 4.0],
        },
        "excluded_parameter_combinations": [
            "Cartesian corners were excluded because simultaneous independent extremes need not be induced by one physical parameter vector.",
            "Points violating kappa_q > 0, vbar_q > 0, sigma_v > 0, or rho in (-1,1) are outside the accepted domain.",
        ],
        "design_seed": DESIGN_SEED,
        "latin_hypercube_points": 20,
        "anchor_points": 4,
        "variance_states": list(VARIANCE_STATES),
        "reference": {
            "method": "Carr-Madan damped Fourier inversion validated high-resolution COS",
            "fourier_integration_limit": 1500.0,
            "fourier_absolute_tolerance": 2e-10,
            "fourier_relative_tolerance": 2e-10,
            "high_resolution_cos_width": REFERENCE_WIDTH,
            "high_resolution_cos_n_cos": REFERENCE_TERMS,
        },
    }


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as file_handle:
        writer = csv.DictWriter(file_handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_calibration(output_directory: Path, basis_path: Path) -> dict[str, object]:
    LOGGER.info("building deterministic calibration design")
    points = calibration_parameter_design()
    output_directory.mkdir(parents=True, exist_ok=True)
    with (output_directory / "calibration_domain.json").open("w") as file_handle:
        json.dump([asdict(point) for point in points], file_handle, indent=2)
    LOGGER.info("computing high-resolution COS reference over the full domain")
    reference_records = _reference_records(points)
    LOGGER.info("cross-checking representative reference prices with independent Fourier inversion")
    crosscheck = _reference_crosscheck(points, reference_records)
    _write_csv(output_directory / "reference_crosscheck.csv", crosscheck)
    LOGGER.info("searching fixed COS candidate configurations")
    summaries = _candidate_summary(points, reference_records)
    _write_csv(output_directory / "calibration_summary.csv", summaries)
    selected = _select_configuration(summaries)
    selected["reference_crosscheck_max_absolute_price_difference"] = max(
        float(row["absolute_price_difference"]) for row in crosscheck
    )
    with (output_directory / "calibration_selected.json").open("w") as file_handle:
        json.dump(selected, file_handle, indent=2)
    production_basis = FixedCosBasisConfig(
        maturities=tuple(selected["maturities"]),
        effective_widths=tuple(selected["effective_widths"]),
        generation_n_cos=int(selected["generation_n_cos"]),
        estimation_n_cos=int(selected["estimation_n_cos"]),
    )
    basis_path.parent.mkdir(parents=True, exist_ok=True)
    with basis_path.open("w") as file_handle:
        json.dump(production_basis.to_dict(), file_handle, indent=2)
        file_handle.write("\n")
    LOGGER.info("selected production basis %s (hash %s)", production_basis.to_dict(), production_basis.stable_hash())
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description="Calibrate the fixed Heston COS production basis.")
    parser.add_argument("--output", type=Path, default=Path("outputs/calibration/fixed_cos"))
    parser.add_argument(
        "--basis-output",
        type=Path,
        default=Path("Python/Scripts/configs/heston_cos_basis_production.json"),
    )
    parser.add_argument("--log-level", default="INFO")
    arguments = parser.parse_args()
    logging.basicConfig(level=getattr(logging, arguments.log_level.upper()), format="%(levelname)s %(message)s")
    run_calibration(arguments.output, arguments.basis_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
