from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(frozen=True)
class OptionPanelDate:
    date_index: int
    time: float
    spot: float
    log_spot: float
    strikes: np.ndarray
    maturities: np.ndarray
    option_types: np.ndarray
    observed_iv: np.ndarray
    rates: np.ndarray
    dividend_yields: np.ndarray
    observed_price: np.ndarray | None = None
    clean_iv: np.ndarray | None = None
    true_variance: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        arrays = (
            self.strikes,
            self.maturities,
            self.option_types,
            self.observed_iv,
            self.rates,
            self.dividend_yields,
        )
        lengths = {np.asarray(value).reshape(-1).size for value in arrays}
        if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
            raise ValueError("all required option-date arrays must be non-empty and aligned")
        n_contracts = next(iter(lengths))
        for name, optional in (("observed_price", self.observed_price), ("clean_iv", self.clean_iv)):
            if optional is not None and np.asarray(optional).reshape(-1).size != n_contracts:
                raise ValueError(f"{name} must align with the required option-date arrays")
        if not np.isfinite(self.time) or not np.isfinite(self.spot) or self.spot <= 0.0:
            raise ValueError("time must be finite and spot must be finite and strictly positive")
        if not np.isfinite(self.log_spot):
            raise ValueError("log_spot must be finite")
        if np.any(np.asarray(self.strikes, dtype=float) <= 0.0):
            raise ValueError("strikes must be strictly positive")
        if np.any(np.asarray(self.maturities, dtype=float) <= 0.0):
            raise ValueError("maturities must be strictly positive")

    @property
    def n_contracts(self) -> int:
        return int(np.asarray(self.observed_iv).shape[0])


@dataclass(frozen=True)
class OptionPanel:
    dates: tuple[OptionPanelDate, ...]
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.dates:
            raise ValueError("option panel must contain at least one date")

    @property
    def log_spot(self) -> np.ndarray:
        return np.array([item.log_spot for item in self.dates], dtype=float)

    @property
    def spot(self) -> np.ndarray:
        return np.array([item.spot for item in self.dates], dtype=float)

    @property
    def times(self) -> np.ndarray:
        return np.array([item.time for item in self.dates], dtype=float)

    @property
    def true_variance(self) -> np.ndarray | None:
        values = [item.true_variance for item in self.dates]
        if any(value is None for value in values):
            return None
        return np.array(values, dtype=float)

    @property
    def n_contracts(self) -> int:
        return int(sum(item.n_contracts for item in self.dates))

    @property
    def n_dates(self) -> int:
        return len(self.dates)

    def first_rate_pair(self) -> tuple[float, float]:
        first = self.dates[0]
        return float(np.mean(first.rates)), float(np.mean(first.dividend_yields))

    def truncate_dates(self, max_dates: int | None) -> OptionPanel:
        if max_dates is None or max_dates >= len(self.dates):
            return self
        if max_dates <= 0:
            raise ValueError("max_dates must be positive")
        return OptionPanel(dates=self.dates[:max_dates], metadata=dict(self.metadata))

