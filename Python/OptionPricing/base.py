from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class CcfSolver(ABC):
    @abstractmethod
    def solve_coefficients(self, u_grid: Any, maturity_grid: Any, model_params: Any, measure_params: Any | None = None):
        """Return coefficient container over transform and maturity grids."""


class OptionPricer(ABC):
    @abstractmethod
    def price_matrix(
        self,
        *,
        log_s: Any,
        variance: Any,
        strike_grid: Any,
        maturity_grid: Any,
        rate_grid: Any,
        coefficients: Any,
        pricing_config: Any,
        dividend_yield_grid: Any | None = None,
        option_type: Any = "call",
    ):
        """Return price matrix over strike and maturity grids."""


class OptionPanelGenerator(ABC):
    @abstractmethod
    def generate_panel(self, log_s_week: Any, v_week: Any):
        """Generate an option panel from state paths."""


class OptionPanelStore(ABC):
    @abstractmethod
    def save(self, panel: Any, file_path: str) -> None:
        """Persist panel."""

    @abstractmethod
    def load(self, file_path: str):
        """Load persisted panel."""
