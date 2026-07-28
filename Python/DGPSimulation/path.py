from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class HestonPath:
    """Weekly Heston states, with daily states when requested."""

    t_week: np.ndarray
    logS_week: np.ndarray
    V_week: np.ndarray
    dlogS_week: np.ndarray
    logS_daily: Optional[np.ndarray] = None
    V_daily: Optional[np.ndarray] = None
    seed: Optional[int] = None
