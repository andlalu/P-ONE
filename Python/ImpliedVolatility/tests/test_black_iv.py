import math

import numpy as np
import pytest

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_price


def test_implied_vol_recovers_clean_black_price():
    price = black76_price(
        forward=100.0,
        strike=103.0,
        tau=0.25,
        vol=0.22,
        discount_factor=math.exp(-0.02 * 0.25),
        option_type="call",
    )
    iv = implied_vol_black76(
        price=price,
        forward=100.0,
        strike=103.0,
        tau=0.25,
        discount_factor=math.exp(-0.02 * 0.25),
        option_type="call",
    )
    assert iv == pytest.approx(0.22, abs=1e-8)


def test_implied_vol_bounds_modes():
    kwargs = {
        "price": 101.0,
        "forward": 100.0,
        "strike": 100.0,
        "tau": 0.25,
        "discount_factor": 1.0,
        "option_type": "call",
    }
    with pytest.raises(ValueError):
        implied_vol_black76(**kwargs, on_bounds="raise")
    assert np.isnan(implied_vol_black76(**kwargs, on_bounds="nan"))
    assert math.isfinite(implied_vol_black76(**kwargs, on_bounds="clip"))
