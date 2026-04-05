from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class ModelParams(ABC):
    """Base class for model parameter containers."""


class SimulationConfig(ABC):
    """Base class for simulation configuration containers."""


class SimulationPath(ABC):
    """Base class for path output containers."""


class VarianceDrawer(ABC):
    """Interface for one-step variance propagation."""

    @abstractmethod
    def draw_next_variance(self, v_n: float, delta: float, rng: Any, params: Any) -> float:
        """Draw or compute the next variance value V_{n+1} from V_n."""


class PathSimulator(ABC):
    """Interface for path simulators."""

    @abstractmethod
    def simulate(self):
        """Run one simulation and return a SimulationPath."""
