from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
from pathlib import Path

from OptionPricing.noisy_panel import (
    NOISE_SCENARIOS,
    NoiseGenerationConfig,
    clean_panel_file,
    generate_noisy_panel_file,
    parse_noise_config,
    write_noise_config,
    write_noisy_manifest,
)
from Scripts.generation import default_workers, set_thread_env


def _load_noise_config(config_path: str | Path) -> NoiseGenerationConfig:
    with Path(config_path).open() as fh:
        raw = json.load(fh)
    return parse_noise_config(raw.get("noise"))


def _task(args: tuple[str, int, str, NoiseGenerationConfig, str, bool]):
    run_root, sample_id, scenario, config, panel_format, skip_existing = args
    clean_path = clean_panel_file(run_root, sample_id)
    return generate_noisy_panel_file(
        clean_panel_path=clean_path,
        run_root=run_root,
        sample_id=sample_id,
        scenario=scenario,
        config=config,
        panel_format=panel_format,
        skip_existing=skip_existing,
    )


def _clean_sample_ids(run_root: str | Path) -> list[int]:
    panel_dir = Path(run_root) / "panels_clean"
    ids = set()
    for path in sorted(panel_dir.glob("sample_*.parquet")) + sorted(panel_dir.glob("sample_*.csv")):
        ids.add(int(path.stem.split("_")[1]))
    return sorted(ids)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--sample-start", type=int, default=0)
    parser.add_argument("--sample-end", type=int, default=None)
    parser.add_argument("--scenarios", default=",".join(NOISE_SCENARIOS))
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    set_thread_env()
    config = _load_noise_config(args.config)
    with Path(args.config).open() as file_handle:
        panel_format = str(json.load(file_handle)["panel_format"])
    requested = tuple(item.strip() for item in args.scenarios.split(",") if item.strip())
    invalid = sorted(set(requested) - set(NOISE_SCENARIOS))
    if invalid:
        parser.error(f"unsupported scenarios: {', '.join(invalid)}")
    enabled = set(config.enabled_scenarios())
    scenarios = tuple(scenario for scenario in requested if scenario in enabled)
    sample_ids = [sample_id for sample_id in _clean_sample_ids(args.run_root) if sample_id >= args.sample_start]
    if args.sample_end is not None:
        sample_ids = [sample_id for sample_id in sample_ids if sample_id < args.sample_end]

    tasks = [
        (args.run_root, sample_id, scenario, config, panel_format, args.skip_existing)
        for sample_id in sample_ids
        for scenario in scenarios
    ]
    if not tasks:
        write_noise_config(args.run_root, config)
        write_noisy_manifest(args.run_root, [])
        return 0

    workers = default_workers(len(tasks), args.workers)
    if workers == 1:
        results = [_task(task) for task in tasks]
    else:
        start_method = "spawn" if "spawn" in mp.get_all_start_methods() else mp.get_all_start_methods()[0]
        ctx = mp.get_context(start_method)
        with ctx.Pool(processes=workers) as pool:
            results = list(pool.imap_unordered(_task, tasks))

    write_noise_config(args.run_root, config)
    write_noisy_manifest(args.run_root, results)
    errors = [result for result in results if result.status == "error"]
    for result in sorted(results, key=lambda item: (item.sample_id, item.noise_scenario)):
        print(
            f"sample={result.sample_id} scenario={result.noise_scenario} status={result.status} "
            f"rows={result.n_rows} capped={result.n_capped_total} output={result.output_observed_panel}"
        )
    if errors:
        for result in errors:
            print(f"sample {result.sample_id} {result.noise_scenario} failed: {result.error_message}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
