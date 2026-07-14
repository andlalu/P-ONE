import math

import numpy as np
import pytest

from Estimation.ISCGMM.implied_state_jacobian import (
    finite_difference_iv_variance_jacobian,
    implied_variance_loss_and_jacobian,
    price_iv_and_initial_variance_jacobian_fixed_basis,
)
from Models.Heston.parameters import HestonRiskNeutralParameters
from OptionPricing.cos_pricer import CosOptionPricer


def _fixture():
    parameters = HestonRiskNeutralParameters(2.0, 0.07875, 0.4, -0.5, 0.02, 0.0)
    maturity = 0.25
    forward = 100.0 * math.exp(parameters.r * maturity)
    strikes = forward * np.exp(np.array([-0.075, 0.0, 0.075]))
    option_types = np.array(["put", "call", "call"])
    basis = CosOptionPricer().prepare_fixed_basis(
        maturity=maturity,
        effective_width=1.25,
        n_cos=512,
        model_params=parameters,
    )
    return parameters, strikes, option_types, basis


def test_numerical_jacobian_uses_central_and_boundary_stencils():
    parameters, strikes, option_types, basis = _fixture()
    common = dict(
        log_spot=math.log(100.0),
        strikes=strikes,
        option_types=option_types,
        rate=parameters.r,
        dividend_yield=parameters.q,
        fixed_cos_basis=basis,
        variance_bounds=(1e-6, 0.15),
        relative_step=1e-4,
        absolute_minimum_step=1e-7,
    )
    interior = finite_difference_iv_variance_jacobian(implied_variance=0.04, **common)
    boundary = finite_difference_iv_variance_jacobian(implied_variance=1e-6, **common)
    assert interior.diagnostics.stencil_type == "central"
    assert interior.diagnostics.function_evaluation_count == 3
    assert boundary.diagnostics.stencil_type == "forward"
    assert boundary.diagnostics.function_evaluation_count == 2


def test_semi_analytical_price_and_iv_derivatives_match_finite_differences():
    parameters, strikes, option_types, basis = _fixture()
    variance = 0.04
    analytical = price_iv_and_initial_variance_jacobian_fixed_basis(
        log_spot=math.log(100.0),
        implied_variance=variance,
        strikes=strikes,
        option_types=option_types,
        rate=parameters.r,
        dividend_yield=parameters.q,
        fixed_cos_basis=basis,
        minimum_black_vega=5e-6,
    )
    numerical = finite_difference_iv_variance_jacobian(
        log_spot=math.log(100.0),
        implied_variance=variance,
        strikes=strikes,
        option_types=option_types,
        rate=parameters.r,
        dividend_yield=parameters.q,
        fixed_cos_basis=basis,
        variance_bounds=(1e-6, 0.15),
        relative_step=1e-5,
        absolute_minimum_step=1e-8,
        minimum_black_vega=5e-6,
    )
    assert not np.any(analytical.failed_contracts)
    np.testing.assert_allclose(
        analytical.initial_variance_jacobian,
        numerical.initial_variance_jacobian,
        rtol=2e-5,
        atol=2e-6,
    )

    step = 1e-6
    pricer = CosOptionPricer()
    option_tensor = option_types.reshape(1, -1, 1)
    upper = pricer.price_matrix_fixed_basis(
        log_s=[math.log(100.0)], variance=[variance + step], strike_grid=strikes,
        rate=parameters.r, dividend_yield=parameters.q, basis=basis, option_type=option_tensor,
    )[0, :, 0]
    lower = pricer.price_matrix_fixed_basis(
        log_s=[math.log(100.0)], variance=[variance - step], strike_grid=strikes,
        rate=parameters.r, dividend_yield=parameters.q, basis=basis, option_type=option_tensor,
    )[0, :, 0]
    np.testing.assert_allclose(analytical.price_jacobian, (upper - lower) / (2.0 * step), rtol=1e-7, atol=1e-7)


def test_low_vega_is_reported_and_loss_jacobian_has_explicit_shape():
    parameters, _, _, basis = _fixture()
    forward = 100.0 * math.exp(parameters.r * basis.maturity)
    analytical = price_iv_and_initial_variance_jacobian_fixed_basis(
        log_spot=math.log(100.0),
        implied_variance=1e-6,
        strikes=np.array([forward * math.exp(0.5)]),
        option_types=np.array(["call"]),
        rate=parameters.r,
        dividend_yield=parameters.q,
        fixed_cos_basis=basis,
        minimum_black_vega=5e-6,
    )
    assert analytical.low_vega_contracts[0] or analytical.price_or_iv_boundary_contracts[0]
    assert analytical.failed_contracts[0]

    loss, derivative = implied_variance_loss_and_jacobian(
        np.array([0.21, 0.19]), np.array([0.20, 0.20]), np.array([2.0, 3.0])
    )
    assert loss == pytest.approx(0.0002)
    assert derivative == pytest.approx(0.02)
    with pytest.raises(ValueError, match="aligned 1D"):
        implied_variance_loss_and_jacobian(np.array([0.2]), np.array([0.2, 0.3]), np.array([1.0]))
