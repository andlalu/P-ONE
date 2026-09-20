from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

for key in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(key, "1")

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from Scripts.experiment_config import load_experiment_config
from Scripts.sample_run import SCENARIO_ORDER, run_sample


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run one complete Heston Monte Carlo sample: generation, validation "
            "and five first-step IS-CGMM estimates."
        )
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-id", required=True, type=int)
    parser.add_argument("--output-root", required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--resume", action="store_true")
    mode.add_argument("--overwrite", action="store_true")
    diagnostic = parser.add_mutually_exclusive_group()
    diagnostic.add_argument("--generation-only", action="store_true")
    diagnostic.add_argument("--state-only", action="store_true")
    parser.add_argument("--scenario", choices=SCENARIO_ORDER, default=None)
    parser.add_argument("--max-dates", type=int, default=None)
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="Allow an intentional Git-SHA mismatch with an existing run.json.",
    )
    arguments = parser.parse_args()

    try:
        config = load_experiment_config(arguments.config)
        run_sample(
            config=config,
            sample_id=arguments.sample_id,
            output_root=arguments.output_root,
            resume=arguments.resume,
            overwrite=arguments.overwrite,
            generation_only=arguments.generation_only,
            state_only=arguments.state_only,
            scenario=arguments.scenario,
            max_dates=arguments.max_dates,
            log_level=arguments.log_level,
            allow_dirty=arguments.allow_dirty,
        )
    except Exception as exception:
        print(f"sample run failed: {type(exception).__name__}: {exception}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
