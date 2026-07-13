from __future__ import annotations

import math
from dataclasses import dataclass, replace


def _validate_common(*, kappa: float, vbar: float, sigma_v: float, rho: float, r: float, q: float) -> None:
    if not math.isfinite(kappa) or kappa <= 0.0:
        raise ValueError("kappa must be finite and strictly positive")
    if not math.isfinite(vbar) or vbar <= 0.0:
        raise ValueError("vbar must be finite and strictly positive")
    if not math.isfinite(sigma_v) or sigma_v <= 0.0:
        raise ValueError("sigma_v must be finite and strictly positive")
    if not math.isfinite(rho) or not -1.0 < rho < 1.0:
        raise ValueError("rho must be finite and in (-1, 1)")
    if not math.isfinite(r) or not math.isfinite(q):
        raise ValueError("r and q must be finite")


@dataclass(frozen=True)
class HestonPhysicalParameters:
    """Heston parameters under the physical measure."""

    eta: float
    kappa: float
    vbar: float
    sigma_v: float
    rho: float
    r: float = 0.0
    q: float = 0.0

    def validate(self) -> None:
        if not math.isfinite(self.eta):
            raise ValueError("eta must be finite")
        _validate_common(
            kappa=self.kappa,
            vbar=self.vbar,
            sigma_v=self.sigma_v,
            rho=self.rho,
            r=self.r,
            q=self.q,
        )


@dataclass(frozen=True)
class HestonRiskNeutralParameters:
    """Heston parameters under the risk-neutral pricing measure."""

    kappa: float
    vbar: float
    sigma_v: float
    rho: float
    r: float = 0.0
    q: float = 0.0

    def validate(self) -> None:
        _validate_common(
            kappa=self.kappa,
            vbar=self.vbar,
            sigma_v=self.sigma_v,
            rho=self.rho,
            r=self.r,
            q=self.q,
        )


@dataclass(frozen=True)
class HestonParameters:
    """Canonical joint Heston parameter vector used by estimation."""

    eta: float
    kappa: float
    vbar: float
    sigma_v: float
    rho: float
    eta_v: float
    r: float = 0.0
    q: float = 0.0

    @property
    def kappa_q(self) -> float:
        return self.kappa - self.eta_v

    @property
    def vbar_q(self) -> float:
        return self.kappa * self.vbar / self.kappa_q

    def validate(self) -> None:
        if not math.isfinite(self.eta) or not math.isfinite(self.eta_v):
            raise ValueError("eta and eta_v must be finite")
        _validate_common(
            kappa=self.kappa,
            vbar=self.vbar,
            sigma_v=self.sigma_v,
            rho=self.rho,
            r=self.r,
            q=self.q,
        )
        if not math.isfinite(self.kappa_q) or self.kappa_q <= 0.0:
            raise ValueError("kappa - eta_v must be finite and strictly positive")

    def to_physical(self) -> HestonPhysicalParameters:
        self.validate()
        return HestonPhysicalParameters(
            eta=self.eta,
            kappa=self.kappa,
            vbar=self.vbar,
            sigma_v=self.sigma_v,
            rho=self.rho,
            r=self.r,
            q=self.q,
        )

    def to_risk_neutral(self) -> HestonRiskNeutralParameters:
        self.validate()
        return HestonRiskNeutralParameters(
            kappa=self.kappa_q,
            vbar=self.vbar_q,
            sigma_v=self.sigma_v,
            rho=self.rho,
            r=self.r,
            q=self.q,
        )

    def with_rates(self, *, r: float, q: float) -> HestonParameters:
        return replace(self, r=float(r), q=float(q))

