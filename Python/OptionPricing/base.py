from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class CcfSolver(ABC):
    @abstractmethod
    def solve_coefficients(self, u_grid: Any, maturity_grid: Any, model_params: Any, measure_params: Any | None = None):
        """Return coefficient container over transform and maturity grids."""


class OptionPricer(ABC):
    @abstractmethod
    def price_matrix_fixed_basis(
        self,
        *,
        log_s: Any,
        variance: Any,
        strike_grid: Any,
        rate: float,
        dividend_yield: float,
        basis: Any,
        option_type: Any = "call",
    ):
        """Return prices using a prepared variance-independent COS basis."""


class OptionPriceCubeGenerator(ABC):
    @abstractmethod
    def generate_panel(self, log_s_week: Any, v_week: Any):
        """Generate an option panel from state paths."""


class OptionPriceCubeStore(ABC):
    @abstractmethod
    def save(self, panel: Any, file_path: str) -> None:
        """Persist panel."""

    @abstractmethod
    def load(self, file_path: str):
        """Load persisted panel."""
