from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from OptionPricing.base import OptionPriceCubeStore
from OptionPricing.types import OptionPriceCube


class NpzOptionPriceCubeStore(OptionPriceCubeStore):
    def save(self, panel: OptionPriceCube, file_path: str) -> None:
        p = Path(file_path)
        np.savez_compressed(
            p,
            prices=panel.prices,
            observation_index=panel.observation_index,
            strikes=panel.strikes,
            maturities=panel.maturities,
            metadata_json=json.dumps(panel.metadata),
        )

    def load(self, file_path: str):
        data = np.load(file_path, allow_pickle=False)
        return OptionPriceCube(
            prices=data["prices"],
            observation_index=data["observation_index"],
            strikes=data["strikes"],
            maturities=data["maturities"],
            metadata=json.loads(str(data["metadata_json"])),
        )
