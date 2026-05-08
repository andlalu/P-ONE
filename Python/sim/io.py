from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Tuple

import numpy as np

from sim.types import HestonParamsP, HestonPath, HestonSimConfig

_FORMAT_VERSION = 1


def save_heston_path_npz(
    file_path: str | Path,
    path: HestonPath,
    params: HestonParamsP,
    config: HestonSimConfig,
    compressed: bool = True,
) -> None:
    """Persist a Heston path and its simulation metadata to .npz."""
    destination = Path(file_path)
    destination.parent.mkdir(parents=True, exist_ok=True)

    metadata = {
        "format_version": _FORMAT_VERSION,
        "seed": path.seed,
        "params": asdict(params),
        "config": asdict(config),
    }

    payload = {
        "metadata_json": np.array(json.dumps(metadata)),
        "t_week": path.t_week,
        "logS_week": path.logS_week,
        "V_week": path.V_week,
        "dlogS_week": path.dlogS_week,
    }

    if path.logS_daily is not None:
        payload["logS_daily"] = path.logS_daily
    if path.V_daily is not None:
        payload["V_daily"] = path.V_daily

    if compressed:
        np.savez_compressed(destination, **payload)
    else:
        np.savez(destination, **payload)


def load_heston_path_npz(file_path: str | Path) -> Tuple[HestonPath, HestonParamsP, HestonSimConfig]:
    """Load a Heston path and metadata previously saved by save_heston_path_npz."""
    source = Path(file_path)

    with np.load(source, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
        if metadata.get("format_version") != _FORMAT_VERSION:
            raise ValueError(
                f"Unsupported heston path format version: {metadata.get('format_version')}"
            )

        params = HestonParamsP(**metadata["params"])
        config = HestonSimConfig(**metadata["config"])

        logS_daily = data["logS_daily"] if "logS_daily" in data.files else None
        V_daily = data["V_daily"] if "V_daily" in data.files else None

        path = HestonPath(
            seed=metadata.get("seed"),
            t_week=data["t_week"],
            logS_week=data["logS_week"],
            V_week=data["V_week"],
            dlogS_week=data["dlogS_week"],
            logS_daily=logS_daily,
            V_daily=V_daily,
        )

    return path, params, config
