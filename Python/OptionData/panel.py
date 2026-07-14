from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

SPOT_LOG_RELATIVE_TOLERANCE = 1e-10
SPOT_LOG_ABSOLUTE_TOLERANCE = 1e-10


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
        required_arrays = {
            "strikes": np.asarray(self.strikes),
            "maturities": np.asarray(self.maturities),
            "option_types": np.asarray(self.option_types),
            "observed_iv": np.asarray(self.observed_iv),
            "rates": np.asarray(self.rates),
            "dividend_yields": np.asarray(self.dividend_yields),
        }
        if any(contract_array.ndim != 1 for contract_array in required_arrays.values()):
            raise ValueError("all required option-date arrays must be one dimensional")
        lengths = {contract_array.size for contract_array in required_arrays.values()}
        if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
            raise ValueError("all required option-date arrays must be non-empty and aligned")
        n_contracts = next(iter(lengths))
        for name, optional in (("observed_price", self.observed_price), ("clean_iv", self.clean_iv)):
            if optional is not None:
                optional_array = np.asarray(optional, dtype=float)
                if optional_array.ndim != 1 or optional_array.size != n_contracts:
                    raise ValueError(f"{name} must be a 1D array aligned with the contract dimension")
                if not np.all(np.isfinite(optional_array)):
                    raise ValueError(f"{name} must contain only finite entries")
        if not np.isfinite(self.time) or not np.isfinite(self.spot) or self.spot <= 0.0:
            raise ValueError("time must be finite and spot must be finite and strictly positive")
        if not np.isfinite(self.log_spot):
            raise ValueError("log_spot must be finite")
        if not np.isclose(
            self.spot,
            np.exp(self.log_spot),
            rtol=SPOT_LOG_RELATIVE_TOLERANCE,
            atol=SPOT_LOG_ABSOLUTE_TOLERANCE,
        ):
            raise ValueError("spot and exp(log_spot) are inconsistent")
        strikes = np.asarray(self.strikes, dtype=float)
        maturities = np.asarray(self.maturities, dtype=float)
        observed_iv = np.asarray(self.observed_iv, dtype=float)
        rates = np.asarray(self.rates, dtype=float)
        dividend_yields = np.asarray(self.dividend_yields, dtype=float)
        if not np.all(np.isfinite(strikes)) or np.any(strikes <= 0.0):
            raise ValueError("strikes must be finite and strictly positive")
        if not np.all(np.isfinite(maturities)) or np.any(maturities <= 0.0):
            raise ValueError("maturities must be finite and strictly positive")
        if not np.all(np.isfinite(observed_iv)) or np.any(observed_iv < 0.0):
            raise ValueError("observed IVs must be finite and non-negative")
        if not np.all(np.isfinite(rates)):
            raise ValueError("rates must contain only finite entries")
        if not np.all(np.isfinite(dividend_yields)):
            raise ValueError("dividend yields must contain only finite entries")
        option_types = np.char.lower(np.asarray(self.option_types).astype(str))
        if not np.all(np.isin(option_types, ("call", "put"))):
            raise ValueError("option types must contain only 'call' or 'put'")
        if self.true_variance is not None and (
            not np.isfinite(self.true_variance) or self.true_variance < 0.0
        ):
            raise ValueError("true_variance must be finite and non-negative when present")

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
        times = np.array([panel_date.time for panel_date in self.dates], dtype=float)
        date_indices = np.array([panel_date.date_index for panel_date in self.dates], dtype=int)
        if np.any(np.diff(times) <= 0.0):
            raise ValueError("option panel times must be strictly increasing")
        if np.unique(date_indices).size != date_indices.size:
            raise ValueError("option panel date indices must be unique")
        if np.any(np.diff(date_indices) <= 0):
            raise ValueError("option panel date indices must be consistently increasing")

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
