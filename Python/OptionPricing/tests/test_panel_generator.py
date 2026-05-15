import numpy as np

from OptionPricing.cos_pricer import CosOptionPricer
from OptionPricing.heston_ccf_solver import HestonAnalyticCcfSolver
from OptionPricing.panel_generator import HestonOptionPanelGenerator
from OptionPricing.types import CosPricingConfig, HestonPricingParamsQ, OptionPanelConfig, PricingStack


def test_panel_generator_outputs_expected_shape():
    stack = PricingStack(ccf_solver=HestonAnalyticCcfSolver(), option_pricer=CosOptionPricer())
    panel_cfg = OptionPanelConfig(
        strikes=np.array([90.0, 100.0]),
        maturities=np.array([0.25, 1.0]),
        rates=np.array([0.01, 0.015]),
    )
    gen = HestonOptionPanelGenerator(
        pricing_stack=stack,
        panel_config=panel_cfg,
        solver_params=HestonPricingParamsQ(kappa=2.0, vbar=0.04, sigma_v=0.5, rho=-0.5),
        pricing_config=CosPricingConfig(n_cos=32),
    )
    panel = gen.generate_panel(np.log(np.array([95.0, 100.0, 105.0])), np.array([0.03, 0.04, 0.05]))
    assert panel.prices.shape == (3, 2, 2)
