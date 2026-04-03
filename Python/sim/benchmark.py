from __future__ import annotations

import time

import numpy as np

from sim.heston_simulator import HestonPathSimulator
from sim.types import HestonParamsP, HestonSimConfig
from sim.variance_drawers import AndersenQeVarianceDrawer


def run_benchmark(repetitions: int = 100_000) -> None:
    params = HestonParamsP(eta=1.5, kappa=3.0, vbar=0.04, sigma_v=0.4, rho=-0.7, r=0.02, q=0.0)
    config = HestonSimConfig(seed=7)

    rng = np.random.default_rng(config.seed)
    drawer = AndersenQeVarianceDrawer()

    v = params.vbar
    t0 = time.perf_counter_ns()
    for _ in range(repetitions):
        v = drawer.draw_next_variance(v, config.delta, rng, params)
    t1 = time.perf_counter_ns()
    ns_step = (t1 - t0) / repetitions

    sim = HestonPathSimulator(params=params, config=config, variance_drawer=drawer)
    t2 = time.perf_counter_ns()
    sim.simulate()
    t3 = time.perf_counter_ns()
    total_days = config.burnin_days + config.m_week * config.t_week
    ns_full_day = (t3 - t2) / total_days

    print(f"QE stepper: {ns_step:,.1f} ns/step")
    print(f"Full daily update: {ns_full_day:,.1f} ns/step")
    print(f"One replication total: {(t3 - t2) / 1e6:,.2f} ms")


if __name__ == "__main__":
    run_benchmark()
