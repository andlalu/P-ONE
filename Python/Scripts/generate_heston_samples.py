from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from Scripts.generation import load_config, run_samples, with_overrides, write_run_metadata


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate Heston state paths and option panels. --sample-end is exclusive."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-start", type=int, default=0)
    parser.add_argument("--sample-end", type=int, default=None, help="Exclusive sample index.")
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--paths-only", action="store_true")
    parser.add_argument("--panels-only", action="store_true")
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
    arguments = parser.parse_args()
    if arguments.paths_only and arguments.panels_only:
        parser.error("--paths-only and --panels-only cannot be combined")
    logging.basicConfig(
        level=getattr(logging, arguments.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    try:
        config = with_overrides(
            load_config(arguments.config),
            workers=arguments.workers,
            skip_existing=arguments.skip_existing,
            paths_only=arguments.paths_only,
            panels_only=arguments.panels_only,
        )
        results = run_samples(
            config=config,
            sample_start=arguments.sample_start,
            sample_end=arguments.sample_end,
        )
        write_run_metadata(config, results)
    except Exception:
        logging.exception("generation setup failed")
        return 1
    for result in sorted(results, key=lambda item: item.sample_id):
        logging.info(
            "sample=%03d status=%s elapsed=%.3fs",
            result.sample_id,
            result.status,
            result.elapsed_seconds,
        )
        if result.status == "error":
            logging.error("sample=%03d failure=%s", result.sample_id, result.error)
    return 1 if any(result.status == "error" for result in results) else 0


if __name__ == "__main__":
    sys.exit(main())
