import numpy as np

from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver
from OptionPricing.types import CoefficientTensor, CosPricingConfig
from OptionPricing.types import HestonPricingParamsQ


def test_price_matrix_shape_and_nonnegative():
    pricer = CosOptionPricer()
    cfg = CosPricingConfig(n_cos=32, truncation_width=10.0)
    state = np.log(np.array([90.0, 100.0]))
    variance = np.array([0.04, 0.04])
    strikes = np.array([80.0, 100.0, 120.0])
    maturities = np.array([0.0, 1.0])
    rates = np.array([0.01, 0.02])
    width = pricer.effective_truncation_width(variance, maturities, cfg)
    u, _ = pricer._get_static_terms(cfg.n_cos, width)
    coef = CoefficientTensor(
        u_grid=u,
        maturities=maturities,
        cf_a=np.zeros((32, 2), dtype=complex),
        cf_b=np.zeros((32, 2), dtype=complex),
    )
    out = pricer.price_matrix(
        log_s=state,
        variance=variance,
        strike_grid=strikes,
        maturity_grid=maturities,
        rate_grid=rates,
        coefficients=coef,
        pricing_config=cfg,
    )
    assert out.shape == (2, 3, 2)
    assert np.all(out >= 0.0)
    np.testing.assert_allclose(out[:, :, 0], np.maximum(np.exp(state)[:, None] - strikes[None, :], 0.0))


def test_heston_prices_change_with_variance():
    pricer = CosOptionPricer()
    cfg = CosPricingConfig(n_cos=128, truncation_width=10.0)
    maturities = np.array([0.25])
    width = pricer.effective_truncation_width(np.array([0.01, 0.09]), maturities, cfg)
    u, _ = pricer._get_static_terms(cfg.n_cos, width)
    params = HestonPricingParamsQ(kappa=3.0, vbar=0.04, sigma_v=0.4, rho=-0.7, r=0.02, q=0.0)
    coef = HestonAnalyticCcfSolver().solve_coefficients(u, maturities, params)

    out = pricer.price_matrix(
        log_s=np.log(np.array([100.0, 100.0])),
        variance=np.array([0.01, 0.09]),
        strike_grid=np.array([100.0]),
        maturity_grid=maturities,
        rate_grid=np.array([0.02]),
        dividend_yield_grid=np.array([0.0]),
        coefficients=coef,
        pricing_config=cfg,
    )
    assert out[1, 0, 0] > out[0, 0, 0]


def test_price_one_reuses_matrix_logic():
    pricer = CosOptionPricer()
    cfg = CosPricingConfig(n_cos=128, truncation_width=10.0)
    params = HestonPricingParamsQ(kappa=3.0, vbar=0.04, sigma_v=0.4, rho=-0.7, r=0.02, q=0.0)
    width = pricer.effective_truncation_width(np.array([0.04]), np.array([0.25]), cfg)
    u, _ = pricer._get_static_terms(cfg.n_cos, width)
    coef = HestonAnalyticCcfSolver().solve_coefficients(u, np.array([0.25]), params)

    one = pricer.price_one(
        S=100.0,
        V=0.04,
        tau=0.25,
        K=101.0,
        option_type="call",
        model_params=params,
        coefficients=coef,
        pricing_config=cfg,
    )
    matrix = pricer.price_matrix(
        log_s=np.array([np.log(100.0)]),
        variance=np.array([0.04]),
        strike_grid=np.array([101.0]),
        maturity_grid=np.array([0.25]),
        rate_grid=np.array([0.02]),
        dividend_yield_grid=np.array([0.0]),
        coefficients=coef,
        pricing_config=cfg,
    )[0, 0, 0]
    np.testing.assert_allclose(one, matrix)
