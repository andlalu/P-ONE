from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import time
from pathlib import Path

import numpy as np

from Estimation.ISCGMM.implied_state_jacobian import (
    finite_difference_iv_variance_jacobian,
    price_iv_and_initial_variance_jacobian_fixed_basis,
)
from OptionPricing.cos_basis import FixedCosBasisConfig
from OptionPricing.cos_pricer import CosOptionPricer
from Scripts.calibrate_fixed_cos_basis import (
    LOG_MONEYNESS,
    MATURITIES,
    SPOT,
    VARIANCE_STATES,
    _contract,
    calibration_parameter_design,
)

LOGGER = logging.getLogger(__name__)


def validate_jacobians(basis_path: Path, output_directory: Path) -> dict[str, object]:
    basis_config = FixedCosBasisConfig.load(basis_path)
    pricer = CosOptionPricer()
    rows: list[dict[str, object]] = []
    for point in calibration_parameter_design():
        parameters = point.risk_neutral_parameters()
        for maturity in MATURITIES:
            strikes = np.array([_contract(maturity, value)[1] for value in LOG_MONEYNESS])
            option_types = np.array([_contract(maturity, value)[2] for value in LOG_MONEYNESS])
            fixed_basis = pricer.prepare_fixed_basis(
                maturity=maturity,
                effective_width=basis_config.width_for_maturity(maturity),
                n_cos=basis_config.estimation_n_cos,
                model_params=parameters,
            )
            for variance in VARIANCE_STATES:
                analytical_started = time.perf_counter()
                analytical = price_iv_and_initial_variance_jacobian_fixed_basis(
                    log_spot=math.log(SPOT),
                    implied_variance=variance,
                    strikes=strikes,
                    option_types=option_types,
                    rate=parameters.r,
                    dividend_yield=parameters.q,
                    fixed_cos_basis=fixed_basis,
                    minimum_black_vega=5e-6,
                )
                analytical_runtime = time.perf_counter() - analytical_started
                for relative_step in (1e-3, 1e-4, 1e-5):
                    finite_difference_started = time.perf_counter()
                    try:
                        numerical = finite_difference_iv_variance_jacobian(
                            log_spot=math.log(SPOT),
                            implied_variance=variance,
                            strikes=strikes,
                            option_types=option_types,
                            rate=parameters.r,
                            dividend_yield=parameters.q,
                            fixed_cos_basis=fixed_basis,
                            variance_bounds=(min(VARIANCE_STATES), max(VARIANCE_STATES)),
                            relative_step=relative_step,
                            absolute_minimum_step=1e-7,
                            minimum_black_vega=5e-6,
                        )
                        numerical_runtime = time.perf_counter() - finite_difference_started
                        numerical_jacobian = numerical.initial_variance_jacobian
                        stencil = numerical.diagnostics.stencil_type
                        evaluation_count = numerical.diagnostics.function_evaluation_count
                        numerical_failed = numerical.diagnostics.failed_contracts
                        lower_point = numerical.diagnostics.lower_evaluation_point
                        upper_point = numerical.diagnostics.upper_evaluation_point
                        option_tensor = option_types.reshape(1, -1, 1)

                        def prices_at(candidate_variance: float) -> np.ndarray:
                            return pricer.price_matrix_fixed_basis(
                                log_s=[math.log(SPOT)],
                                variance=[candidate_variance],
                                strike_grid=strikes,
                                rate=parameters.r,
                                dividend_yield=parameters.q,
                                basis=fixed_basis,
                                option_type=option_tensor,
                            )[0, :, 0]

                        if stencil == "central":
                            numerical_price_jacobian = (
                                prices_at(upper_point) - prices_at(lower_point)
                            ) / (upper_point - lower_point)
                        elif stencil == "forward":
                            numerical_price_jacobian = (
                                prices_at(upper_point) - analytical.prices
                            ) / (upper_point - variance)
                        else:
                            numerical_price_jacobian = (
                                analytical.prices - prices_at(lower_point)
                            ) / (variance - lower_point)
                    except (ValueError, FloatingPointError, OverflowError):
                        numerical_runtime = time.perf_counter() - finite_difference_started
                        numerical_jacobian = np.full(strikes.size, np.nan)
                        stencil = "failed"
                        evaluation_count = 0
                        numerical_failed = np.ones(strikes.size, dtype=bool)
                        numerical_price_jacobian = np.full(strikes.size, np.nan)
                    for contract_index, log_moneyness in enumerate(LOG_MONEYNESS):
                        analytical_value = analytical.initial_variance_jacobian[contract_index]
                        numerical_value = numerical_jacobian[contract_index]
                        absolute_error = abs(analytical_value - numerical_value)
                        denominator = max(abs(analytical_value), abs(numerical_value), 1e-10)
                        price_absolute_error = abs(
                            analytical.price_jacobian[contract_index]
                            - numerical_price_jacobian[contract_index]
                        )
                        price_relative_error = price_absolute_error / max(
                            abs(analytical.price_jacobian[contract_index]),
                            abs(numerical_price_jacobian[contract_index]),
                            1e-10,
                        )
                        rows.append(
                            {
                                "point_id": point.point_id,
                                "maturity": maturity,
                                "variance": variance,
                                "log_moneyness": log_moneyness,
                                "relative_step": relative_step,
                                "stencil_type": stencil,
                                "analytical_iv_jacobian": analytical_value,
                                "finite_difference_iv_jacobian": numerical_value,
                                "absolute_derivative_error": absolute_error,
                                "relative_derivative_error": absolute_error / denominator,
                                "analytical_price_jacobian": analytical.price_jacobian[contract_index],
                                "finite_difference_price_jacobian": numerical_price_jacobian[contract_index],
                                "absolute_price_derivative_error": price_absolute_error,
                                "relative_price_derivative_error": price_relative_error,
                                "black_vega": analytical.vegas[contract_index],
                                "low_vega": bool(analytical.low_vega_contracts[contract_index]),
                                "boundary_case": bool(
                                    variance in {min(VARIANCE_STATES), max(VARIANCE_STATES)}
                                ),
                                "analytical_failed": bool(analytical.failed_contracts[contract_index]),
                                "finite_difference_failed": bool(numerical_failed[contract_index]),
                                "analytical_runtime_seconds": analytical_runtime,
                                "finite_difference_runtime_seconds": numerical_runtime,
                                "finite_difference_evaluation_count": evaluation_count,
                            }
                        )
    output_directory.mkdir(parents=True, exist_ok=True)
    csv_path = output_directory / "jacobian_validation.csv"
    with csv_path.open("w", newline="") as file_handle:
        writer = csv.DictWriter(file_handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    comparable = [
        row
        for row in rows
        if not row["analytical_failed"]
        and not row["finite_difference_failed"]
        and row["relative_step"] == 1e-4
    ]
    absolute_errors = np.array([float(row["absolute_derivative_error"]) for row in comparable])
    relative_errors = np.array([float(row["relative_derivative_error"]) for row in comparable])
    price_derivative_errors = np.array(
        [float(row["absolute_price_derivative_error"]) for row in comparable]
    )
    price_derivative_relative_errors = np.array(
        [float(row["relative_price_derivative_error"]) for row in comparable]
    )
    interior = [row for row in comparable if not row["boundary_case"]]
    boundary = [row for row in comparable if row["boundary_case"]]
    summary: dict[str, object] = {
        "basis_path": str(basis_path.resolve()),
        "basis_hash": basis_config.stable_hash(),
        "record_count": len(rows),
        "comparison_count_at_selected_step": len(comparable),
        "maximum_absolute_derivative_error": float(np.max(absolute_errors)),
        "p99_absolute_derivative_error": float(np.quantile(absolute_errors, 0.99)),
        "maximum_relative_derivative_error": float(np.max(relative_errors)),
        "p99_relative_derivative_error": float(np.quantile(relative_errors, 0.99)),
        "maximum_absolute_price_derivative_error": float(np.max(price_derivative_errors)),
        "p99_absolute_price_derivative_error": float(np.quantile(price_derivative_errors, 0.99)),
        "maximum_relative_price_derivative_error": float(np.max(price_derivative_relative_errors)),
        "p99_relative_price_derivative_error": float(
            np.quantile(price_derivative_relative_errors, 0.99)
        ),
        "interior_maximum_absolute_iv_derivative_error": max(
            float(row["absolute_derivative_error"]) for row in interior
        ),
        "interior_p99_absolute_iv_derivative_error": float(
            np.quantile([float(row["absolute_derivative_error"]) for row in interior], 0.99)
        ),
        "boundary_maximum_absolute_iv_derivative_error": max(
            float(row["absolute_derivative_error"]) for row in boundary
        ),
        "boundary_p99_absolute_iv_derivative_error": float(
            np.quantile([float(row["absolute_derivative_error"]) for row in boundary], 0.99)
        ),
        "low_vega_record_count": sum(bool(row["low_vega"]) for row in rows),
        "boundary_record_count": sum(bool(row["boundary_case"]) for row in rows),
        "analytical_failure_count": sum(bool(row["analytical_failed"]) for row in rows),
        "finite_difference_failure_count": sum(bool(row["finite_difference_failed"]) for row in rows),
        "selected_relative_step": 1e-4,
        "absolute_minimum_step": 1e-7,
        "minimum_black_vega": 5e-6,
        "analytical_runtime_seconds": sum(float(row["analytical_runtime_seconds"]) for row in rows)
        / (len(LOG_MONEYNESS) * 3),
        "finite_difference_runtime_seconds": sum(
            float(row["finite_difference_runtime_seconds"]) for row in rows
        )
        / len(LOG_MONEYNESS),
    }
    summary["acceptance_tolerances"] = {
        "maximum_relative_iv_derivative_error": 0.003,
        "p99_absolute_iv_derivative_error": 0.0002,
        "maximum_relative_price_derivative_error": 0.003,
        "p99_relative_price_derivative_error": 5e-5,
    }
    summary["accepted"] = bool(
        summary["maximum_relative_derivative_error"] <= 0.003
        and summary["p99_absolute_derivative_error"] <= 0.0002
        and summary["maximum_relative_price_derivative_error"] <= 0.003
        and summary["p99_relative_price_derivative_error"] <= 5e-5
    )
    with (output_directory / "jacobian_validation_summary.json").open("w") as file_handle:
        json.dump(summary, file_handle, indent=2)
        file_handle.write("\n")
    LOGGER.info("Jacobian validation summary: %s", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate semi-analytical implied-state Jacobians.")
    parser.add_argument(
        "--basis",
        type=Path,
        default=Path("Python/Scripts/configs/heston_cos_basis_production.json"),
    )
    parser.add_argument("--output", type=Path, default=Path("outputs/calibration/fixed_cos"))
    parser.add_argument("--log-level", default="INFO")
    arguments = parser.parse_args()
    logging.basicConfig(level=getattr(logging, arguments.log_level.upper()), format="%(levelname)s %(message)s")
    validate_jacobians(arguments.basis, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
