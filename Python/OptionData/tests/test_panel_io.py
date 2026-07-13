import csv

import pytest

from OptionData.io import load_option_panel, write_panel_metadata
from OptionPricing.cos_basis import FixedCosBasisConfig


def _rows():
    rows = []
    for week in range(3):
        rows.append(
            {
                "week_index": week,
                "t": week / 52.0,
                "S": 100.0 + week,
                "logS": 4.605170186 + week * 0.01,
                "V": 0.04,
                "maturity_years": 0.25,
                "strike": 100.0,
                "option_type": "call",
                "r": 0.02,
                "q": 0.0,
                "model_iv": 0.2,
                "model_price": 4.0,
            }
        )
    return rows


def _metadata():
    return {"sample_id": 0, "scenario": "clean", "cos_basis": FixedCosBasisConfig((0.25,), (1.5,), 32).to_metadata()}


def test_canonical_panel_loading_from_csv_and_metadata_roundtrip(tmp_path):
    target = tmp_path / "panel.csv"
    with target.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(_rows()[0]))
        writer.writeheader()
        writer.writerows(_rows())
    write_panel_metadata(target, _metadata())
    panel = load_option_panel(target)
    assert panel.n_dates == 3
    assert panel.n_contracts == 3
    assert panel.metadata["cos_basis"] == _metadata()["cos_basis"]
    assert panel.dates[0].clean_iv[0] == pytest.approx(0.2)


def test_canonical_panel_loading_from_parquet(tmp_path):
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pyarrow")
    target = tmp_path / "panel.parquet"
    pd.DataFrame(_rows()).to_parquet(target, index=False)
    write_panel_metadata(target, _metadata())
    panel = load_option_panel(target, max_dates=2)
    assert panel.n_dates == 2

