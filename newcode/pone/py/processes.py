import numpy as np
from parameters import ParameterSet
from simulation import Simulator, SimulationType

class StochProch:
    def __init__(self, dimension: int,   time: float, state: np.array, parameters: ParameterSet):
        self.dimension = dimension
        self.time = time
        if len(state) != dimension:
            raise ValueError("Incorrect state size")
        self.state = state
        self.parameters = parameters
        self.sim_states = None

    # simulates
    def simulate(self,
                 simulator: Simulator,
                 n_timepoints: int,
                 delta_t: float,
                 n_sims_per_timepoint: int,
                 n_versions: int):
        return NotImplementedError("Not implemented")

    # loads pre-simulated states
    def load_simulated_states(self):
        return NotImplementedError("Not implemented")


class Heston(StochProch):
    def __init__(self,
                 time: float,
                 state: np.array,
                 kappa: float,
                 v_bar: float,
                 sigma_v: float,
                 rho: float,
                 eta: float,
                 eta_v: float):
        super().__init__(2,
                         time,
                         state,
                         ParameterSet(
                             4,
                             np.array([kappa, v_bar, sigma_v, rho]),
                             1,
                             np.array([eta]),
                             1,
                             np.array([eta_v])))
    def kappa(self):
        return self.parameters.param()[0]

    def v_bar(self):
        return self.parameters.param()[1]

    def sigma_v(self):
        return self.parameters.param()[2]

    def rho(self):
        return self.parameters.param()[3]

    def eta(self):
        return self.parameters.p_param()[0]

    def eta_v(self):
        return self.parameters.q_param()[0]

    def simulate(self,
                 simulator: Simulator,
                 n_timepoints: int,
                 delta_t: float,
                 n_sims_per_timepoint: int,
                 n_versions: int):
        std_dev = np.sqrt(delta_t/n_sims_per_timepoint)
        if std_dev < 1e-20:
            raise ValueError("Process standard deviation is too small/negative")
        simulated = delta_t * simulator.bivariate_norm(correlation=self.rho(), no_draws=(n_timepoints * n_sims_per_timepoint, n_versions))
        return simulated
