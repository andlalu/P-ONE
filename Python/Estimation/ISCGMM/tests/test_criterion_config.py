import numpy as np
import pytest

from Estimation.ISCGMM.cgmm_criterion import _infer_constant_transition_interval, make_cgmm_quadrature
from Estimation.ISCGMM.config import CcfQuadratureConfig, CgmmConfig, ImpliedStateConfig
from OptionData.panel import OptionPanel, OptionPanelDate
from OptionPricing.cos_basis import FixedCosBasisConfig


def _panel(times):
    basis = FixedCosBasisConfig((0.25,), (1.5,), 16)
    dates = tuple(
        OptionPanelDate(
            date_index=index,
            time=time,
            spot=100.0,
            log_spot=np.log(100.0),
            strikes=np.array([100.0]),
            maturities=np.array([0.25]),
            option_types=np.array(["call"]),
            observed_iv=np.array([0.2]),
            rates=np.array([0.0]),
            dividend_yields=np.array([0.0]),
        )
        for index, time in enumerate(times)
    )
    return OptionPanel(dates, metadata={"cos_basis": basis.to_metadata()}), CgmmConfig(ImpliedStateConfig(basis))


def test_constant_spacing_acceptance_and_irregular_rejection():
    regular, config = _panel((0.0, 0.1, 0.2, 0.3))
    assert _infer_constant_transition_interval(regular, config) == pytest.approx(0.1)
    irregular, config = _panel((0.0, 0.1, 0.21, 0.31))
    with pytest.raises(ValueError, match="equally spaced"):
        _infer_constant_transition_interval(irregular, config)


def test_gauss_hermite_quadrature_shape_determinism_and_normalisation():
    config = CcfQuadratureConfig(dimension=2, order=3, scale=1.25)
    nodes_a, weights_a = make_cgmm_quadrature(config)
    nodes_b, weights_b = make_cgmm_quadrature(config)
    assert nodes_a.shape == (9, 2)
    assert weights_a.shape == (9,)
    np.testing.assert_array_equal(nodes_a, nodes_b)
    np.testing.assert_array_equal(weights_a, weights_b)
    assert np.sum(weights_a) == pytest.approx(1.0)

