from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
import time
from pathlib import Path

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from Estimation.ISCGMM.config import ImpliedStateConfig
from Estimation.ISCGMM.implied_state import imply_heston_variance_path
from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_price
from Models.Heston.parameters import HestonParameters
from OptionData.panel import OptionPanel, OptionPanelDate
from OptionPricing.cos_basis import FixedCosBasisConfig, cos_specification_metadata
from OptionPricing.cos_pricer import CosOptionPricer
from Scripts.experiment_config import load_experiment_config

LOGGER = logging.getLogger(__name__)
SOLVERS = (
    "golden_section",
    "bounded_brent",
    "least_squares_finite_difference",
    "least_squares_analytic",
)
SCENARIOS = ("clean", "low_iid", "spatial_corr", "persistent_factor")
VARIANCE_PATH = np.array([1e-6, 0.0025, 0.0225, 0.08, 0.149, 0.04])


def _theta() -> HestonParameters:
    return HestonParameters(5.0, 7.0, 0.0225, 0.4, -0.5, 5.0, 0.02, 0.0)


def _clean_panel(basis_config: FixedCosBasisConfig) -> OptionPanel:
    theta = _theta()
    parameters = theta.to_risk_neutral()
    pricer = CosOptionPricer()
    bases = {
        maturity: pricer.prepare_fixed_basis(
            maturity=maturity,
            effective_width=basis_config.width_for_maturity(maturity),
            n_cos=basis_config.estimation_n_cos,
            model_params=parameters,
        )
        for maturity in basis_config.maturities
    }
    dates: list[OptionPanelDate] = []
    for date_index, variance in enumerate(VARIANCE_PATH):
        spot = 100.0 * math.exp(0.002 * date_index)
        strikes: list[float] = []
        maturities: list[float] = []
        option_types: list[str] = []
        implied_volatilities: list[float] = []
        prices: list[float] = []
        for maturity in basis_config.maturities:
            forward = spot * math.exp(parameters.r * maturity)
            local_strikes = forward * np.exp(np.array([-0.15, -0.075, 0.0, 0.075, 0.15]))
            local_types = np.array(["put", "put", "call", "call", "call"])
            local_prices = pricer.price_matrix_fixed_basis(
                log_s=[math.log(spot)],
                variance=[variance],
                strike_grid=local_strikes,
                rate=parameters.r,
                dividend_yield=parameters.q,
                basis=bases[maturity],
                option_type=local_types.reshape(1, -1, 1),
            )[0, :, 0]
            discount = math.exp(-parameters.r * maturity)
            for strike, option_type, price in zip(local_strikes, local_types, local_prices):
                strikes.append(float(strike))
                maturities.append(float(maturity))
                option_types.append(str(option_type))
                prices.append(float(price))
                implied_volatilities.append(
                    implied_vol_black76(
                        price=float(price),
                        forward=forward,
                        strike=float(strike),
                        tau=maturity,
                        discount_factor=discount,
                        option_type=str(option_type),
                        on_bounds="clip",
                    )
                )
        contract_count = len(strikes)
        dates.append(
            OptionPanelDate(
                date_index=date_index,
                time=date_index * 5.0 / 252.0,
                spot=spot,
                log_spot=math.log(spot),
                strikes=np.asarray(strikes),
                maturities=np.asarray(maturities),
                option_types=np.asarray(option_types),
                observed_iv=np.asarray(implied_volatilities),
                rates=np.full(contract_count, parameters.r),
                dividend_yields=np.full(contract_count, parameters.q),
                observed_price=np.asarray(prices),
                clean_iv=np.asarray(implied_volatilities),
                true_variance=float(variance),
            )
        )
    return OptionPanel(tuple(dates), metadata={"cos_basis": cos_specification_metadata(basis_config)})


def _contaminated_panel(clean_panel: OptionPanel, scenario: str) -> OptionPanel:
    if scenario == "clean":
        return clean_panel
    rng = np.random.default_rng({"low_iid": 301, "spatial_corr": 302, "persistent_factor": 303}[scenario])
    factor = 0.0
    dates: list[OptionPanelDate] = []
    for panel_date in clean_panel.dates:
        clean_iv = np.asarray(panel_date.clean_iv)
        if scenario == "low_iid":
            noise = 0.0015 * rng.standard_normal(panel_date.n_contracts)
        elif scenario == "spatial_corr":
            common = rng.normal()
            noise = 0.0015 * (0.7 * common + math.sqrt(1.0 - 0.7**2) * rng.standard_normal(panel_date.n_contracts))
        else:
            factor = 0.85 * factor + rng.normal(scale=0.0007)
            noise = factor + 0.0013 * rng.standard_normal(panel_date.n_contracts)
        contaminated_iv = np.maximum(0.0001, clean_iv + noise)
        observed_price = np.empty(panel_date.n_contracts)
        observed_iv = np.empty(panel_date.n_contracts)
        for index, (maturity, strike, option_type, volatility, rate) in enumerate(
            zip(
                panel_date.maturities,
                panel_date.strikes,
                panel_date.option_types,
                contaminated_iv,
                panel_date.rates,
            )
        ):
            forward = panel_date.spot * math.exp(rate * maturity)
            discount = math.exp(-rate * maturity)
            raw_price = black76_price(
                forward=forward,
                strike=float(strike),
                tau=float(maturity),
                vol=float(volatility),
                discount_factor=discount,
                option_type=str(option_type),
            )
            rounded = 0.01 * round(raw_price / 0.01)
            intrinsic = discount * max(
                forward - strike if option_type == "call" else strike - forward, 0.0
            )
            upper = discount * (forward if option_type == "call" else strike)
            observed_price[index] = min(upper - 1e-10, max(intrinsic + 1e-10, rounded))
            observed_iv[index] = max(
                0.0001,
                implied_vol_black76(
                    price=float(observed_price[index]),
                    forward=forward,
                    strike=float(strike),
                    tau=float(maturity),
                    discount_factor=discount,
                    option_type=str(option_type),
                    on_bounds="clip",
                ),
            )
        dates.append(
            OptionPanelDate(
                date_index=panel_date.date_index,
                time=panel_date.time,
                spot=panel_date.spot,
                log_spot=panel_date.log_spot,
                strikes=panel_date.strikes,
                maturities=panel_date.maturities,
                option_types=panel_date.option_types,
                observed_iv=observed_iv,
                rates=panel_date.rates,
                dividend_yields=panel_date.dividend_yields,
                observed_price=observed_price,
                clean_iv=panel_date.clean_iv,
                true_variance=panel_date.true_variance,
            )
        )
    return OptionPanel(tuple(dates), metadata={**clean_panel.metadata, "scenario": scenario})


def run_benchmark(config_path: Path, output_directory: Path) -> dict[str, object]:
    basis_config = load_experiment_config(config_path).cos_basis
    clean = _clean_panel(basis_config)
    panels = {scenario: _contaminated_panel(clean, scenario) for scenario in SCENARIOS}
    references: dict[str, np.ndarray] = {}
    for scenario, panel in panels.items():
        references[scenario] = imply_heston_variance_path(
            _theta(),
            panel,
            ImpliedStateConfig(
                basis_config,
                v_min=1e-6,
                v_max=0.15,
                tol=1e-10,
                max_iter=160,
                boundary_tol=1e-6,
                warm_start_window=None,
                state_solver="golden_section",
                fallback_solver=None,
                minimum_black_vega=5e-6,
            ),
        ).variance
    rows: list[dict[str, object]] = []
    for scenario, panel in panels.items():
        for solver in SOLVERS:
            for warm_start_window in (None, 0.04):
                config = ImpliedStateConfig(
                    basis_config,
                    v_min=1e-6,
                    v_max=0.15,
                    tol=1e-7,
                    max_iter=100,
                    boundary_tol=1e-5,
                    warm_start_window=warm_start_window,
                    state_solver=solver,
                    fallback_solver="golden_section" if solver.startswith("least_squares") else None,
                    finite_difference_relative_step=1e-4,
                    finite_difference_absolute_step=1e-7,
                    minimum_black_vega=5e-6,
                )
                started = time.perf_counter()
                result = imply_heston_variance_path(_theta(), panel, config)
                runtime = time.perf_counter() - started
                deviation = np.abs(result.variance - references[scenario])
                true_deviation = np.abs(result.variance - VARIANCE_PATH)
                rows.append(
                    {
                        "scenario": scenario,
                        "solver": solver,
                        "warm_start_window": "none" if warm_start_window is None else warm_start_window,
                        "maximum_reference_deviation": float(np.max(deviation)),
                        "rmse_from_reference": float(np.sqrt(np.mean(deviation**2))),
                        "rmse_from_true_variance": float(np.sqrt(np.mean(true_deviation**2))),
                        "mean_final_loss": float(np.mean(result.objective)),
                        "objective_evaluations": int(np.sum(result.nfev)),
                        "jacobian_evaluations": result.jacobian_evaluation_count,
                        "total_pricing_calls": int(np.sum(result.nfev)) + result.jacobian_evaluation_count,
                        "runtime_seconds": runtime,
                        "failure_count": int(np.sum(result.failed)),
                        "boundary_hit_count": int(np.sum(result.boundary_hit)),
                        "fallback_count": result.fallback_count,
                    }
                )
                LOGGER.info("benchmarked %s / %s / warm=%s", scenario, solver, warm_start_window)
    output_directory.mkdir(parents=True, exist_ok=True)
    with (output_directory / "state_solver_benchmark.csv").open("w", newline="") as file_handle:
        writer = csv.DictWriter(file_handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    aggregate: dict[str, dict[str, float]] = {}
    for solver in SOLVERS:
        solver_rows = [row for row in rows if row["solver"] == solver]
        aggregate[solver] = {
            "maximum_reference_deviation": max(float(row["maximum_reference_deviation"]) for row in solver_rows),
            "failure_count": sum(int(row["failure_count"]) for row in solver_rows),
            "fallback_count": sum(int(row["fallback_count"]) for row in solver_rows),
            "median_runtime_seconds": float(np.median([float(row["runtime_seconds"]) for row in solver_rows])),
            "median_pricing_calls": float(np.median([int(row["total_pricing_calls"]) for row in solver_rows])),
        }
    eligible_scalar = [solver for solver in ("bounded_brent", "golden_section") if aggregate[solver]["failure_count"] == 0]
    selected_solver = min(
        eligible_scalar,
        key=lambda solver: (
            aggregate[solver]["maximum_reference_deviation"] > 2e-6,
            aggregate[solver]["median_runtime_seconds"],
        ),
    )
    summary: dict[str, object] = {
        "selected_production_solver": selected_solver,
        "fallback_solver": "golden_section" if selected_solver != "golden_section" else None,
        "selection_reason": "The scalar solver had zero failures on clean and three rounded/capped contaminations and avoided low-vega Jacobian fallbacks; the eligible scalar method with lower median runtime was selected.",
        "deterministic_seeds": {"low_iid": 301, "spatial_corr": 302, "persistent_factor": 303},
        "variance_dates": VARIANCE_PATH.tolist(),
        "aggregate": aggregate,
    }
    with (output_directory / "state_solver_benchmark_summary.json").open("w") as file_handle:
        json.dump(summary, file_handle, indent=2)
        file_handle.write("\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark implied-state solvers on clean and noisy panels.")
    parser.add_argument(
        "--config", type=Path, default=Path("Python/Scripts/configs/heston_experiment_run_001.json")
    )
    parser.add_argument("--output", type=Path, default=Path("outputs/calibration/fixed_cos"))
    parser.add_argument("--log-level", default="INFO")
    arguments = parser.parse_args()
    logging.basicConfig(level=getattr(logging, arguments.log_level.upper()), format="%(levelname)s %(message)s")
    run_benchmark(arguments.config, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
