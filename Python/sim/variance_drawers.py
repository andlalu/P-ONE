from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from sim.base import VarianceDrawer
from sim.types import HestonParamsP


@dataclass(frozen=True)
class AndersenQeVarianceDrawer(VarianceDrawer):
    psi_c: float = 1.5
    eps: float = 1e-16

    def draw_next_variance(self, v_n: float, delta: float, rng: Any, params: HestonParamsP) -> float:
        e = math.exp(-params.kappa * delta)
        m = params.vbar + (v_n - params.vbar) * e
        s2 = (
            v_n * (params.sigma_v * params.sigma_v) * e * (1.0 - e) / params.kappa
            + params.vbar * (params.sigma_v * params.sigma_v) * ((1.0 - e) ** 2) / (2.0 * params.kappa)
        )

        m = max(m, self.eps)
        psi = s2 / (m * m)

        if psi <= self.psi_c:
            psi_eff = max(psi, self.eps)
            two_over_psi = 2.0 / psi_eff
            b2 = two_over_psi - 1.0 + math.sqrt(two_over_psi) * math.sqrt(max(two_over_psi - 1.0, 0.0))
            a = m / (1.0 + b2)
            z = rng.normal()
            v_next = a * ((math.sqrt(b2) + z) ** 2)
            return max(v_next, 0.0)

        p = (psi - 1.0) / (psi + 1.0)
        p = min(max(p, 0.0), 1.0 - self.eps)
        beta = (1.0 - p) / m
        u = rng.random()
        if u <= p:
            return 0.0

        u2 = (u - p) / (1.0 - p)
        u2 = max(u2, self.eps)
        return max(-(1.0 / beta) * math.log(u2), 0.0)


@dataclass(frozen=True)
class EulerVarianceDrawer(VarianceDrawer):
    """Euler-Maruyama with full truncation to keep variance non-negative."""

    def draw_next_variance(self, v_n: float, delta: float, rng: Any, params: HestonParamsP) -> float:
        v_pos = max(v_n, 0.0)
        z = rng.normal()
        dv = params.kappa * (params.vbar - v_pos) * delta + params.sigma_v * math.sqrt(v_pos * delta) * z
        return max(v_n + dv, 0.0)
