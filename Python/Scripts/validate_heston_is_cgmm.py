from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from Estimation.ISCGMM.estimate import estimate_first_step
from Estimation.ISCGMM.implied_state import imply_heston_variance_path
from OptionData.io import load_option_panel
from Scripts.experiment_config import load_experiment_config

LOGGER = logging.getLogger(__name__)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one local fixed-basis Heston first-step IS-CGMM job.")
    parser.add_argument("--config", required=True, help="Unified experiment configuration.")
    parser.add_argument("--panel", required=True, help="Input clean or contaminated CSV/Parquet panel.")
    parser.add_argument("--iv-column", choices=("model_iv", "clean_iv", "observed_iv"), default=None)
    parser.add_argument("--max-dates", type=int, default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--state-only",
        action="store_true",
        help="Run only the clean implied-state inversion at the configured base start.",
    )
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    experiment = load_experiment_config(args.config)
    panel = load_option_panel(args.panel, iv_column=args.iv_column, max_dates=args.max_dates)
    if args.state_only:
        implied_state = imply_heston_variance_path(
            experiment.optimizer_config.base_start,
            panel,
            experiment.criterion_config.implied_state,
        )
        true_variance = panel.true_variance
        state_rmse = None
        if true_variance is not None:
            state_rmse = float(np.sqrt(np.mean((implied_state.variance - true_variance) ** 2)))
        payload = {
            "mode": "state_only",
            "base_start": asdict(experiment.optimizer_config.base_start),
            "state_rmse": state_rmse,
            "implied_state": implied_state.to_dict(),
        }
        success = implied_state.success_rate == 1.0
    else:
        estimate = estimate_first_step(
            panel,
            criterion_config=experiment.criterion_config,
            optimizer_config=experiment.optimizer_config,
        )
        payload = estimate.to_dict()
        success = estimate.success
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as file_handle:
        json.dump(payload, file_handle, indent=2)
        file_handle.write("\n")
    LOGGER.info("result written path=%s", output)
    return 0 if success else 2


if __name__ == "__main__":
    raise SystemExit(main())
