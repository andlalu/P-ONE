from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from DGPSimulation.base import SimulationConfig, SimulationPath


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
    t_week: np.ndarray
    logS_week: np.ndarray
    V_week: np.ndarray
    dlogS_week: np.ndarray
    logS_daily: Optional[np.ndarray] = None
    V_daily: Optional[np.ndarray] = None
    seed: Optional[int] = None

