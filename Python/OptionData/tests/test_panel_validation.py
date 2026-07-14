import math
from dataclasses import replace

import numpy as np
import pytest

from OptionData.panel import OptionPanel, OptionPanelDate


def _date(index=0, time=0.0):
    return OptionPanelDate(
        date_index=index,
        time=time,
        spot=100.0,
        log_spot=math.log(100.0),
        strikes=np.array([95.0, 105.0]),
        maturities=np.array([0.25, 0.25]),
        option_types=np.array(["put", "call"]),
        observed_iv=np.array([0.2, 0.21]),
        rates=np.array([0.02, 0.02]),
        dividend_yields=np.array([0.0, 0.0]),
        true_variance=0.04,
    )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("observed_iv", np.array([0.2, np.nan]), "observed IVs"),
        ("option_types", np.array(["put", "straddle"]), "option types"),
        ("rates", np.array([0.02, np.inf]), "rates"),
        ("dividend_yields", np.array([0.0, np.nan]), "dividend yields"),
        ("true_variance", -0.1, "true_variance"),
    ],
)
def test_option_panel_date_rejects_invalid_contract_data(field, value, message):
    with pytest.raises(ValueError, match=message):
        replace(_date(), **{field: value})


def test_spot_and_log_spot_must_be_consistent():
    with pytest.raises(ValueError, match="inconsistent"):
        replace(_date(), log_spot=math.log(101.0))


def test_panel_times_and_indices_must_be_strictly_increasing():
    with pytest.raises(ValueError, match="times"):
        OptionPanel((_date(0, 0.0), _date(1, 0.0)))
    with pytest.raises(ValueError, match="indices"):
        OptionPanel((_date(1, 0.0), _date(0, 0.1)))
