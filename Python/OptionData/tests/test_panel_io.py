import csv
import math

import pytest

from OptionData.io import load_option_panel, write_panel_metadata
from OptionData.noise_common import NOISE_SCENARIOS
from OptionPricing.config import FixedCosBasisConfig
from OptionPricing.cos_basis import cos_specification_metadata


def _rows():
    rows = []
    for week in range(3):
        rows.append(
            {
                "week_index": week,
                "t": week / 52.0,
                "S": 100.0 + week,
                "logS": math.log(100.0 + week),
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
    return {
        "sample_id": 0,
        "scenario": "clean",
        "cos_basis": cos_specification_metadata(FixedCosBasisConfig((0.25,), (1.5,), 32, 16)),
    }


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


def test_combined_panel_loading_supports_variance_linked_factor(tmp_path):
    target = tmp_path / "combined.csv"
    rows = []
    for scenario in ("clean",) + NOISE_SCENARIOS:
        row = dict(_rows()[0])
        row.update(
            sample_id=0,
            scenario=scenario,
            log_moneyness=0.0,
            estimation_iv=0.21,
            estimation_price=4.1,
        )
        rows.append(row)
    with target.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    panel = load_option_panel(target, scenario="variance_linked_factor")

    assert panel.n_dates == 1
    assert panel.n_contracts == 1
    assert panel.metadata["scenario"] == "variance_linked_factor"
    assert panel.dates[0].observed_iv[0] == pytest.approx(0.21)
