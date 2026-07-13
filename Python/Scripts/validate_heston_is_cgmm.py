from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from Estimation.ISCGMM.config import (
    CcfQuadratureConfig,
    CgmmConfig,
    ImpliedStateConfig,
    LoggingConfig,
    OptimizerConfig,
    PowellStageConfig,
)
from Estimation.ISCGMM.estimate import estimate_first_step
from Models.Heston.parameters import HestonParameters
from OptionData.io import load_option_panel
from OptionPricing.cos_basis import FixedCosBasisConfig

LOGGER = logging.getLogger(__name__)


def _load_json(path: str | Path) -> dict:
    with Path(path).open() as fh:
        value = json.load(fh)
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _load_basis(path: str | Path) -> FixedCosBasisConfig:
    raw = _load_json(path)
    config = FixedCosBasisConfig(
        maturities=tuple(float(value) for value in raw["maturities"]),
        effective_widths=tuple(float(value) for value in raw["effective_widths"]),
        n_cos=int(raw["n_cos"]),
        maturity_tolerance=float(raw.get("maturity_tolerance", 1e-10)),
    )
    config.validate()
    return config


def _load_profile(path: str | Path, basis: FixedCosBasisConfig) -> tuple[CgmmConfig, OptimizerConfig, LoggingConfig]:
    raw = _load_json(path)
    state = raw.get("implied_state", {})
    quadrature = raw.get("quadrature", {})
    cgmm = raw.get("cgmm", {})
    optimizer = raw["optimizer"]
    stage1 = optimizer["stage1"]
    stage2 = optimizer["stage2"]
    criterion_config = CgmmConfig(
        implied_state=ImpliedStateConfig(
            cos_basis=basis,
            v_min=float(state.get("v_min", 1e-8)),
            v_max=float(state.get("v_max", 1.0)),
            tol=float(state.get("tol", 2e-4)),
            max_iter=int(state.get("max_iter", 40)),
            boundary_tol=float(state.get("boundary_tol", 1e-5)),
            warm_start_window=state.get("warm_start_window", 0.25),
        ),
        quadrature=CcfQuadratureConfig(
            dimension=int(quadrature.get("dimension", 2)),
            order=int(quadrature.get("order", 3)),
            scale=float(quadrature.get("scale", 1.0)),
        ),
        instrument_precision=tuple(float(value) for value in cgmm.get("instrument_precision", [10.0, 50.0])),
        transition_rk_steps=int(cgmm.get("transition_rk_steps", 32)),
        transition_cf_method=str(cgmm.get("transition_cf_method", "analytic")),
        dt=None if cgmm.get("dt") is None else float(cgmm["dt"]),
        spacing_tolerance=float(cgmm.get("spacing_tolerance", 1e-10)),
    )
    optimizer_config = OptimizerConfig(
        base_start=HestonParameters(**optimizer["base_start"]),
        natural_bounds=tuple(tuple(float(value) for value in pair) for pair in optimizer["natural_bounds"]),
        candidate_relative_perturbations=tuple(
            tuple(float(value) for value in item) for item in optimizer.get("candidate_relative_perturbations", [])
        ),
        stage1=PowellStageConfig(
            max_evaluations=int(stage1["max_evaluations"]),
            xtol=float(stage1["xtol"]),
            ftol=float(stage1["ftol"]),
        ),
        stage2=PowellStageConfig(
            max_evaluations=int(stage2["max_evaluations"]),
            xtol=float(stage2["xtol"]),
            ftol=float(stage2["ftol"]),
        ),
        penalty_value=float(optimizer.get("penalty_value", 1e12)),
    )
    log_config = LoggingConfig(progress_every=int(raw.get("logging", {}).get("progress_every", 10)))
    return criterion_config, optimizer_config, log_config


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one local fixed-basis Heston first-step IS-CGMM job.")
    parser.add_argument("--panel", required=True, help="Input clean or contaminated CSV/Parquet panel.")
    parser.add_argument("--iv-column", choices=("model_iv", "clean_iv", "observed_iv"), default=None)
    parser.add_argument("--max-dates", type=int, default=None)
    parser.add_argument("--cos-basis-config", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    basis = _load_basis(args.cos_basis_config)
    criterion_config, optimizer_config, log_config = _load_profile(args.profile, basis)
    panel = load_option_panel(args.panel, iv_column=args.iv_column, max_dates=args.max_dates)
    estimate = estimate_first_step(
        panel,
        criterion_config=criterion_config,
        optimizer_config=optimizer_config,
        logging_config=log_config,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as fh:
        json.dump(estimate.to_dict(), fh, indent=2)
    LOGGER.info("result written path=%s", output)
    return 0 if estimate.success else 2


if __name__ == "__main__":
    raise SystemExit(main())
