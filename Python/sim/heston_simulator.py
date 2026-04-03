from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from sim.base import PathSimulator, VarianceDrawer
from sim.types import HestonParamsP, HestonPath, HestonSimConfig
from sim.variance_drawers import AndersenQeVarianceDrawer


@dataclass
class HestonPathSimulator(PathSimulator):
    params: HestonParamsP
    config: HestonSimConfig
    variance_drawer: VarianceDrawer = AndersenQeVarianceDrawer()

    def simulate(self) -> HestonPath:
        self.params.validate()
        self.config.validate()

        rng = np.random.default_rng(self.config.seed)

        total_days = self.config.burnin_days + self.config.m_week * self.config.t_week

        v0 = self.params.vbar if self.config.v0 is None else self.config.v0
        log_s0 = math.log(self.config.s0)

        logS_all = np.empty(total_days + 1, dtype=float)
        V_all = np.empty(total_days + 1, dtype=float)
        logS_all[0] = log_s0
        V_all[0] = v0

        sqrt_one_minus_rho_sq = math.sqrt(max(1.0 - self.params.rho * self.params.rho, 0.0))

        for n in range(total_days):
            v_n = V_all[n]
            v_np1 = self.variance_drawer.draw_next_variance(v_n, self.config.delta, rng, self.params)
            v_np1 = max(v_np1, 0.0)

            i_hat = 0.5 * self.config.delta * (v_n + v_np1)
            i_hat = max(i_hat, 0.0)

            z_perp = rng.normal()
            dy = (self.params.r - self.params.q) * self.config.delta
            dy += (self.params.eta - 0.5) * i_hat
            dy += (self.params.rho / self.params.sigma_v) * (
                (v_np1 - v_n) - self.params.kappa * (self.params.vbar * self.config.delta - i_hat)
            )
            dy += sqrt_one_minus_rho_sq * math.sqrt(i_hat) * z_perp

            logS_all[n + 1] = logS_all[n] + dy
            V_all[n + 1] = v_np1

        week_indices = self.config.burnin_days + np.arange(
            0, self.config.m_week * self.config.t_week + 1, self.config.m_week
        )

        logS_week = logS_all[week_indices]
        V_week = V_all[week_indices]
        t_week = np.arange(self.config.t_week + 1, dtype=float) * self.config.m_week * self.config.delta
        dlogS_week = np.diff(logS_week)

        if self.config.return_daily:
            return HestonPath(
                t_week=t_week,
                logS_week=logS_week,
                V_week=V_week,
                dlogS_week=dlogS_week,
                logS_daily=logS_all,
                V_daily=V_all,
            )

        return HestonPath(t_week=t_week, logS_week=logS_week, V_week=V_week, dlogS_week=dlogS_week)
