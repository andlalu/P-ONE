import json
import logging

import numpy as np

import Estimation.ISCGMM.estimate as estimate_module
from Estimation.ISCGMM.config import (
    CgmmConfig,
    ImpliedStateConfig,
    LoggingConfig,
    OptimizerConfig,
    PowellStageConfig,
)
from Estimation.ISCGMM.results import CriterionDiagnostics, ImpliedStateResult
from Models.Heston.parameters import HestonParameters
from OptionData.panel import OptionPanel, OptionPanelDate
from OptionPricing.cos_basis import FixedCosBasisConfig


def _basis():
    return FixedCosBasisConfig((0.25,), (1.5,), 16, 16)


def _panel():
    dates = tuple(
        OptionPanelDate(
            date_index=index,
            time=index * 0.1,
            spot=100.0 * np.exp(index * 0.01),
            log_spot=np.log(100.0) + index * 0.01,
            strikes=np.array([100.0]),
            maturities=np.array([0.25]),
            option_types=np.array(["call"]),
            observed_iv=np.array([0.2]),
            rates=np.array([0.02]),
            dividend_yields=np.array([0.0]),
        )
        for index in range(3)
    )
    return OptionPanel(
        dates,
        metadata={"sample_id": 0, "scenario": "clean", "cos_basis": _basis().generation_metadata()},
    )


def _optimizer():
    return OptimizerConfig(
        base_start=HestonParameters(5.0, 7.0, 0.0225, 0.4, -0.5, 5.0, 0.02, 0.0),
        natural_bounds=((0.0, 10.0), (2.0, 12.0), (0.005, 0.08), (0.15, 0.8), (-0.9, -0.1), (0.5, 6.0)),
        candidate_relative_perturbations=((0.1, -0.1, 0.1, 0.1, 0.1, 0.1), (-0.1, 0.1, -0.1, -0.1, -0.1, -0.1)),
        stage1=PowellStageConfig(30, 1e-2, 1e-3),
        stage2=PowellStageConfig(40, 1e-4, 1e-5),
    )


class _QuadraticCriterion:
    def __init__(self, panel, config):
        self.panel = panel
        self.config = config

    def evaluate_free(self, free, *, r, q):
        target = np.array([4.5, np.log(6.5), np.log(0.025), np.log(0.38), np.arctanh(-0.55), np.log(2.2)])
        return float(np.sum((np.asarray(free) - target) ** 2))

    def evaluate(self, theta, *, return_diagnostics=False):
        value = self.evaluate_free(estimate_module.to_free(theta), r=theta.r, q=theta.q)
        implied = ImpliedStateResult(
            variance=np.full(3, theta.vbar),
            objective=np.zeros(3),
            boundary_hit=np.zeros(3, dtype=bool),
            failed=np.zeros(3, dtype=bool),
            nfev=np.ones(3, dtype=int),
            start_values=np.full(3, theta.vbar),
        )
        return CriterionDiagnostics(value, 3, 1, 3, implied, {"min": theta.vbar, "max": theta.vbar, "mean": theta.vbar, "std": 0.0}, None, None, {"total": 0.0}, {"cos_basis": self.config.implied_state.cos_basis.to_dict()})


class _PenaltyCriterion(_QuadraticCriterion):
    def evaluate_free(self, free, *, r, q):
        if float(np.asarray(free)[0]) > 5.2:
            raise ValueError("expected transformed-parameter failure")
        return super().evaluate_free(free, r=r, q=q)


def test_powell_reproducibility_bounds_start_limit_logging_and_serialisation(monkeypatch, caplog):
    monkeypatch.setattr(estimate_module, "CgmmFirstStepCriterion", _QuadraticCriterion)
    config = CgmmConfig(ImpliedStateConfig(_basis()))
    caplog.set_level(logging.INFO)
    first = estimate_module.estimate_first_step(
        _panel(), criterion_config=config, optimizer_config=_optimizer(), logging_config=LoggingConfig(progress_every=5)
    )
    second = estimate_module.estimate_first_step(
        _panel(), criterion_config=config, optimizer_config=_optimizer(), logging_config=LoggingConfig(progress_every=5)
    )
    np.testing.assert_allclose(first.free_parameters, second.free_parameters, rtol=0.0, atol=0.0)
    assert first.final_criterion == second.final_criterion
    assert len(first.candidate_starts) == 3
    natural = np.array([
        first.estimated_parameters.eta,
        first.estimated_parameters.kappa,
        first.estimated_parameters.vbar,
        first.estimated_parameters.sigma_v,
        first.estimated_parameters.rho,
        first.estimated_parameters.kappa_q,
    ])
    bounds = np.array(_optimizer().natural_bounds)
    assert np.all(natural >= bounds[:, 0]) and np.all(natural <= bounds[:, 1])
    assert any("Powell progress" in record.message for record in caplog.records)
    assert not any("implied variance" in record.message.lower() for record in caplog.records)
    json.dumps(first.to_dict())


def test_powell_expected_numerical_failures_use_finite_penalty(monkeypatch):
    monkeypatch.setattr(estimate_module, "CgmmFirstStepCriterion", _PenaltyCriterion)
    result = estimate_module.estimate_first_step(
        _panel(),
        criterion_config=CgmmConfig(ImpliedStateConfig(_basis())),
        optimizer_config=_optimizer(),
        logging_config=LoggingConfig(progress_every=100),
    )
    assert result.penalty_evaluations > 0
    assert np.isfinite(result.final_criterion)
