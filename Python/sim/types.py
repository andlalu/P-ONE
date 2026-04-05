from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from sim.base import ModelParams, SimulationConfig, SimulationPath


@dataclass(frozen=True)
class HestonParamsP(ModelParams):
    eta: float
    kappa: float
    vbar: float
    sigma_v: float
    rho: float
    r: float = 0.0
    q: float = 0.0

    def validate(self) -> None:
        if self.kappa <= 0.0:
            raise ValueError("kappa must be strictly positive")
        if self.vbar < 0.0:
            raise ValueError("vbar must be non-negative")
        if self.sigma_v <= 0.0:
            raise ValueError("sigma_v must be strictly positive")
        if not -1.0 <= self.rho <= 1.0:
            raise ValueError("rho must be in [-1, 1]")


@dataclass(frozen=True)
class HestonSimConfig(SimulationConfig):
    delta: float = 1.0 / 252.0
    m_week: int = 5
    t_week: int = 525
    burnin_days: int = 756
    s0: float = 100.0
    v0: Optional[float] = None
    seed: int = 12345
    return_daily: bool = False

    def validate(self) -> None:
        if self.delta <= 0.0:
            raise ValueError("delta must be strictly positive")
        if self.m_week <= 0:
            raise ValueError("m_week must be a positive integer")
        if self.t_week <= 0:
            raise ValueError("t_week must be a positive integer")
        if self.burnin_days < 0:
            raise ValueError("burnin_days must be non-negative")
        if self.s0 <= 0.0:
            raise ValueError("s0 must be strictly positive")
        if self.v0 is not None and self.v0 < 0.0:
            raise ValueError("v0 must be non-negative when provided")


@dataclass(frozen=True)
class HestonPath(SimulationPath):
    seed: Optional[int]
    t_week: np.ndarray
    logS_week: np.ndarray
    V_week: np.ndarray
    dlogS_week: np.ndarray
    logS_daily: Optional[np.ndarray] = None
    V_daily: Optional[np.ndarray] = None
