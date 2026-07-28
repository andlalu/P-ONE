import numpy as np

from Models.Heston.parameters import HestonRiskNeutralParameters
from OptionPricing.cos_pricer import CosOptionPricer


def _basis(n_cos: int = 64):
    parameters = HestonRiskNeutralParameters(
        kappa=3.0,
        vbar=0.04,
        sigma_v=0.4,
        rho=-0.7,
        r=0.02,
        q=0.0,
    )
    return CosOptionPricer().prepare_fixed_basis(
        maturity=0.25,
        effective_width=1.25,
        n_cos=n_cos,
        model_params=parameters,
    )


def test_fixed_basis_price_matrix_shape_and_nonnegative():
    prices = CosOptionPricer().price_matrix_fixed_basis(
        log_s=np.log(np.array([90.0, 100.0])),
        variance=np.array([0.04, 0.04]),
        strike_grid=np.array([80.0, 100.0, 120.0]),
        rate=0.02,
        dividend_yield=0.0,
        basis=_basis(),
    )
    assert prices.shape == (2, 3, 1)
    assert np.all(prices >= 0.0)


def test_heston_prices_change_with_variance():
    prices = CosOptionPricer().price_matrix_fixed_basis(
        log_s=np.log(np.array([100.0, 100.0])),
        variance=np.array([0.01, 0.09]),
        strike_grid=np.array([100.0]),
        rate=0.02,
        dividend_yield=0.0,
        basis=_basis(128),
    )
    assert prices[1, 0, 0] > prices[0, 0, 0]
