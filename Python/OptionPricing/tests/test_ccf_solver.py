import numpy as np

from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver
from OptionPricing.types import HestonPricingParamsQ


def test_ccf_shapes_and_finiteness():
    solver = HestonAnalyticCcfSolver()
    params = HestonPricingParamsQ(kappa=2.0, vbar=0.04, sigma_v=0.5, rho=-0.7)
    u = np.linspace(0.0, 50.0, 64)
    t = np.array([0.1, 0.5, 1.0])
    coef = solver.solve_coefficients(u, t, params)
    assert coef.cf_values.shape == (64, 3)
    assert np.all(np.isfinite(coef.cf_values.real))
    assert np.all(np.isfinite(coef.cf_values.imag))
