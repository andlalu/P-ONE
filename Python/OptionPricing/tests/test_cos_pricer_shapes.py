import numpy as np

from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.types import CoefficientTensor, CosPricingConfig


def test_price_matrix_shape_and_nonnegative():
    pricer = CosOptionPricer()
    cfg = CosPricingConfig(n_cos=32)
    state = np.log(np.array([90.0, 100.0]))
    strikes = np.array([80.0, 100.0, 120.0])
    maturities = np.array([0.25, 1.0])
    rates = np.array([0.01, 0.02])
    coef = CoefficientTensor(u_grid=np.array([0.0]), maturities=maturities, cf_values=np.ones((1, 2), dtype=complex))
    out = pricer.price_matrix(state, strikes, maturities, rates, coef, cfg)
    assert out.shape == (2, 3, 2)
    assert np.all(out >= 0.0)
