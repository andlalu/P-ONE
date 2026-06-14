from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.types import CosPricingConfig

from Estimation.ISCGMM.cgmm_criterion import CgmmFirstStepCriterion, criterion_diagnostics_to_dict
from Estimation.ISCGMM.panel import load_option_panel_data
from Estimation.ISCGMM.types import CgmmConfig, HestonEstimationParams, ImpliedStateConfig, QuadratureConfig
from Scripts.generation import generate_one_sample, load_config, with_overrides


def _sample_panel_path(run_root: str | Path, sample_id: int, scenario: str) -> Path:
    root = Path(run_root)
    stem = f"sample_{sample_id:03d}"
    if scenario == "clean":
        base = root / "panels_clean"
    else:
        base = root / "panels_observed" / scenario
    for suffix in (".parquet", ".csv"):
        candidate = base / f"{stem}{suffix}"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"missing {scenario} panel {stem} under {base}")


def _true_theta(config) -> HestonEstimationParams:
    return HestonEstimationParams(
        eta=config.dgp.eta,
        kappa=config.dgp.kappa,
        vbar=config.dgp.vbar,
        sigma_v=config.dgp.sigma_v,
        rho=config.dgp.rho,
        eta_v=config.eta_v,
        r=config.dgp.r,
        q=config.dgp.q,
    )


def _perturb(theta: HestonEstimationParams, *, label: str) -> HestonEstimationParams:
    if label == "mild":
        kappa = theta.kappa * 0.80
        kappa_q = theta.kappa_q * 1.20
        return HestonEstimationParams(
            eta=theta.eta * 1.20,
            kappa=kappa,
            vbar=theta.vbar * 1.20,
            sigma_v=theta.sigma_v * 1.20,
            rho=max(-0.95, min(0.95, theta.rho - 0.20)),
            eta_v=kappa - kappa_q,
            r=theta.r,
            q=theta.q,
        )
    if label == "strong":
        kappa = theta.kappa * 0.70
        kappa_q = theta.kappa_q * 1.35
        return HestonEstimationParams(
            eta=theta.eta * 1.25,
            kappa=kappa,
            vbar=theta.vbar * 1.30,
            sigma_v=theta.sigma_v * 0.75,
            rho=max(-0.95, min(0.95, theta.rho - 0.20)),
            eta_v=kappa - kappa_q,
            r=theta.r,
            q=theta.q,
        )
    raise ValueError(f"unknown perturbation label: {label}")


def _panel_maturities(panel) -> np.ndarray:
    return np.unique(np.concatenate([date.maturities for date in panel.dates]))


def _resolve_effective_width(args: argparse.Namespace, panel, pricing_config: CosPricingConfig) -> tuple[float | None, dict]:
    mode = args.effective_width_mode
    if mode == "dynamic":
        return None, {"mode": mode, "value": None}
    if mode == "fixed":
        if args.effective_truncation_width is None:
            raise ValueError("--effective-truncation-width is required for --effective-width-mode fixed")
        return float(args.effective_truncation_width), {"mode": mode, "value": float(args.effective_truncation_width)}

    maturities = _panel_maturities(panel)
    if mode == "clean-panel-diagnostic":
        true_variance = panel.true_variance
        if true_variance is None:
            raise ValueError("clean-panel-diagnostic width mode requires true variance in the panel")
        width = CosOptionPricer.effective_truncation_width(true_variance, maturities, pricing_config)
        return width, {"mode": mode, "value": width, "source": "panel_true_variance"}
    if mode == "conservative":
        width = CosOptionPricer.effective_truncation_width(np.array([args.v_max], dtype=float), maturities, pricing_config)
        return width, {"mode": mode, "value": width, "source": "v_max"}
    raise ValueError(f"unsupported effective width mode: {mode}")


def run_validation(args: argparse.Namespace) -> dict:
    config = load_config(args.config)
    config = with_overrides(config, workers=1)
    try:
        panel_path = _sample_panel_path(config.output_root, args.sample_id, args.scenario)
    except FileNotFoundError:
        if not args.generate_if_missing:
            raise
        config = with_overrides(config, n_samples=max(config.n_samples, args.sample_id + 1), workers=1)
        generation_config = replace(config, noise=None) if args.scenario == "clean" else config
        result = generate_one_sample(args.sample_id, generation_config)
        if result.status != "ok" and result.status != "skipped":
            raise RuntimeError(result.error)
        panel_path = _sample_panel_path(config.output_root, args.sample_id, args.scenario)

    max_dates = None if args.full_panel else args.max_dates
    panel = load_option_panel_data(panel_path, max_dates=max_dates)
    theta_true = _true_theta(config)
    pricing_config = CosPricingConfig(n_cos=args.cos_terms, truncation_width=args.cos_truncation)
    effective_width, effective_width_info = _resolve_effective_width(args, panel, pricing_config)
    criterion_config = CgmmConfig(
        implied_state=ImpliedStateConfig(
            v_min=args.v_min,
            v_max=args.v_max,
            tol=args.state_tol,
            max_iter=args.state_max_iter,
            effective_truncation_width=effective_width,
            pricing_config=pricing_config,
        ),
        quadrature=QuadratureConfig(order=args.quad_order, scale=args.quad_scale),
        instrument_precision=(args.instrument_return_precision, args.instrument_variance_precision),
        transition_rk_steps=args.transition_rk_steps,
        transition_cf_method=args.transition_cf_method,
    )
    criterion = CgmmFirstStepCriterion(panel, criterion_config)

    evaluations = []
    for label, theta in (
        ("true", theta_true),
        ("mild_perturbation", _perturb(theta_true, label="mild")),
        ("strong_perturbation", _perturb(theta_true, label="strong")),
    ):
        diagnostics = criterion.evaluate(theta, return_diagnostics=True)
        evaluations.append(
            {
                "label": label,
                **criterion_diagnostics_to_dict(theta, diagnostics),
            }
        )

    summary = {
        "config": str(args.config),
        "panel_path": str(panel_path),
        "scenario": args.scenario,
        "sample_id": args.sample_id,
        "panel": {
            "n_dates": panel.n_dates,
            "n_contracts": panel.n_contracts,
            "full_panel": bool(args.full_panel),
            "max_dates": None if args.full_panel else args.max_dates,
        },
        "settings": {
            "cos_terms": args.cos_terms,
            "cos_truncation": args.cos_truncation,
            "effective_width": effective_width_info,
            "transition_cf_method": args.transition_cf_method,
            "transition_rk_steps": args.transition_rk_steps,
            "state_tol": args.state_tol,
            "state_max_iter": args.state_max_iter,
            "quad_order": args.quad_order,
            "quad_scale": args.quad_scale,
            "instrument_precision": [args.instrument_return_precision, args.instrument_variance_precision],
        },
        "evaluations": evaluations,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a local Heston implied-state C-GMM validation.")
    parser.add_argument("--config", default="Python/Scripts/configs/clean_generation_run_001.json")
    parser.add_argument("--sample-id", type=int, default=0)
    parser.add_argument("--scenario", default="clean")
    parser.add_argument("--generate-if-missing", action="store_true")
    parser.add_argument("--max-dates", type=int, default=30)
    parser.add_argument("--full-panel", action="store_true", help="Use all weekly dates in the selected generated panel.")
    parser.add_argument("--output", default="outputs/estimation/is_cgmm_validation_sample_000.json")
    parser.add_argument("--cos-terms", type=int, default=256)
    parser.add_argument("--cos-truncation", type=float, default=10.0)
    parser.add_argument(
        "--effective-width-mode",
        choices=("dynamic", "fixed", "clean-panel-diagnostic", "conservative"),
        default="dynamic",
    )
    parser.add_argument("--effective-truncation-width", type=float, default=None)
    parser.add_argument("--quad-order", type=int, default=3)
    parser.add_argument("--quad-scale", type=float, default=1.0)
    parser.add_argument("--transition-rk-steps", type=int, default=32)
    parser.add_argument("--transition-cf-method", choices=("analytic", "rk4"), default="analytic")
    parser.add_argument("--state-tol", type=float, default=2e-4)
    parser.add_argument("--state-max-iter", type=int, default=40)
    parser.add_argument("--v-min", type=float, default=1e-8)
    parser.add_argument("--v-max", type=float, default=1.0)
    parser.add_argument("--instrument-return-precision", type=float, default=10.0)
    parser.add_argument("--instrument-variance-precision", type=float, default=50.0)
    args = parser.parse_args()

    summary = run_validation(args)
    print(json.dumps({item["label"]: item["criterion_value"] for item in summary["evaluations"]}, indent=2))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
