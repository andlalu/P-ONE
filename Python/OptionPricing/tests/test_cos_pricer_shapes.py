import numpy as np

from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.types import CoefficientTensor, CosPricingConfig


def test_price_matrix_shape_and_nonnegative():
    pricer = CosOptionPricer()
    cfg = CosPricingConfig(n_cos=32, truncation_width=10.0)
    state = np.log(np.array([90.0, 100.0]))
    strikes = np.array([80.0, 100.0, 120.0])
    maturities = np.array([0.0, 1.0])
    rates = np.array([0.01, 0.02])
    coef = CoefficientTensor(
        u_grid=np.linspace(0.0, 5.0, 32),
        maturities=maturities,
        cf_values=np.exp(-0.1 * np.outer(np.linspace(0.0, 5.0, 32) ** 2, maturities)).astype(complex),
    )
    out = pricer.price_matrix(state, strikes, maturities, rates, coef, cfg)
    assert out.shape == (2, 3, 2)
    assert np.all(out >= 0.0)
    np.testing.assert_allclose(out[:, :, 0], np.maximum(np.exp(state)[:, None] - strikes[None, :], 0.0))
