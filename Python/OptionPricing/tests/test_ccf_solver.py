import numpy as np

from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver
from Models.Heston.parameters import HestonRiskNeutralParameters


def test_ccf_shapes_and_finiteness():
    solver = HestonAnalyticCcfSolver()
    params = HestonRiskNeutralParameters(kappa=2.0, vbar=0.04, sigma_v=0.5, rho=-0.7)
    u = np.linspace(0.0, 50.0, 64)
    t = np.array([0.1, 0.5, 1.0])
    coef = solver.solve_coefficients(u, t, params)
    assert coef.cf_a.shape == (64, 3)
    assert coef.cf_b.shape == (64, 3)
    assert np.all(np.isfinite(coef.cf_a.real))
    assert np.all(np.isfinite(coef.cf_a.imag))
    assert np.all(np.isfinite(coef.cf_b.real))
    assert np.all(np.isfinite(coef.cf_b.imag))
