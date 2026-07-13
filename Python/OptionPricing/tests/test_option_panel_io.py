import numpy as np

from OptionPricing.io import NpzOptionPriceCubeStore
from OptionPricing.types import OptionPriceCube


def test_option_panel_npz_roundtrip(tmp_path):
    panel = OptionPriceCube(
        prices=np.ones((2, 3, 4)),
        observation_index=np.array([0, 1]),
        strikes=np.array([80.0, 100.0, 120.0]),
        maturities=np.array([0.25, 0.5, 1.0, 2.0]),
        metadata={"seed": 7},
    )
    store = NpzOptionPriceCubeStore()
    target = tmp_path / "panel.npz"
    store.save(panel, str(target))
    loaded = store.load(str(target))
    np.testing.assert_allclose(loaded.prices, panel.prices)
    assert loaded.metadata["seed"] == 7
