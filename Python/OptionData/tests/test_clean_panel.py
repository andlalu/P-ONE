import numpy as np

from DGPSimulation.path import HestonPath
from Models.Heston.parameters import HestonParameters, HestonPhysicalParameters
from OptionData.clean_panel import generate_clean_option_panel_rows
from OptionPricing.config import FixedCosBasisConfig
from OptionPricing.cos_pricer import CosOptionPricer


def test_clean_panel_uses_canonical_physical_parameter_construction(monkeypatch):
    physical = HestonPhysicalParameters(
        eta=1.5,
        kappa=3.0,
        vbar=0.04,
        sigma_v=0.4,
        rho=-0.7,
        r=0.02,
        q=0.0,
    )
    path = HestonPath(
        t_week=np.array([0.0]),
        logS_week=np.log(np.array([100.0])),
        V_week=np.array([0.04]),
        dlogS_week=np.array([]),
    )
    calls = []
    original = HestonParameters.from_physical

    def counting_constructor(cls, params, *, eta_v):
        calls.append((params, eta_v))
        return original(params, eta_v=eta_v)

    monkeypatch.setattr(
        HestonParameters,
        "from_physical",
        classmethod(counting_constructor),
    )
    rows = generate_clean_option_panel_rows(
        run_id="test",
        sample_id=0,
        path=path,
        params_p=physical,
        eta_v=0.0,
        maturities_years=(0.25,),
        log_moneyness=(-0.05, 0.0, 0.05),
        atm_option_type="call",
        pricing_method="COS",
        iv_method="lets_be_rational",
        pricer=CosOptionPricer(),
        cos_basis=FixedCosBasisConfig((0.25,), (1.25,), 64, 32),
    )
    assert calls == [(physical, 0.0)]
    assert len(rows) == 3
