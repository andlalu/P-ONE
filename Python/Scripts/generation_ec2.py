from __future__ import annotations

import argparse
import sys

from Scripts.generation import load_config, run_samples, with_overrides, write_noisy_metadata, write_run_metadata


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--output-root", default=None)
    parser.add_argument("--sample-start", type=int, default=0)
    parser.add_argument("--sample-end", type=int, default=None)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--paths-only", action="store_true")
    parser.add_argument("--panels-only", action="store_true")
    args = parser.parse_args()

    if args.paths_only and args.panels_only:
        parser.error("--paths-only and --panels-only cannot both be set")

    config = with_overrides(
        load_config(args.config),
        workers=args.workers,
        output_root=args.output_root,
        skip_existing=args.skip_existing,
        paths_only=args.paths_only,
        panels_only=args.panels_only,
    )
    results = run_samples(config=config, sample_start=args.sample_start, sample_end=args.sample_end)
    write_run_metadata(config, results)
    write_noisy_metadata(config, results)

    errors = [result for result in results if result.status == "error"]
    for result in sorted(results, key=lambda item: item.sample_id):
        print(
            f"sample={result.sample_id} status={result.status} "
            f"elapsed={result.elapsed_seconds:.3f}s path={result.path_file} panel={result.panel_file}"
        )
    if errors:
        for result in errors:
            print(f"sample {result.sample_id} failed: {result.error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
