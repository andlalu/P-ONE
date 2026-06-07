from __future__ import annotations

import argparse
import sys

from Scripts.generation import load_config, run_samples, with_overrides, write_noisy_metadata, write_run_metadata
from Scripts.validate_generation import validate_run


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--n-samples", type=int, default=2)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--output-root", default="outputs/generation/local_run")
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    config = with_overrides(
        load_config(args.config),
        n_samples=args.n_samples,
        workers=args.workers,
        output_root=args.output_root,
        skip_existing=args.skip_existing,
    )
    results = run_samples(config=config)
    write_run_metadata(config, results)
    write_noisy_metadata(config, results)
    errors = [result for result in results if result.status == "error"]
    if errors:
        for result in errors:
            print(f"sample {result.sample_id} failed: {result.error}", file=sys.stderr)
        return 1
    validate_run(run_root=args.output_root, expected_samples=args.n_samples)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
