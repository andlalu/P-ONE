from enum import Enum
from ponecpp import pone
import numpy as np

class SimulationType(Enum):
    NotSet = -1,
    Euler = 0,


class SimulatorType(Enum):
    NotSet = -1,
    CppPone = 0,
    NumPy = 1,  # to be implemented
    CuPy = 2,  # to be implemented


class Simulator:
    def _switch_rand_gen(self) -> object:
        __switcher = {
            SimulatorType.NotSet: None,
            SimulatorType.CppPone: pone.Rand(self.__seed),
            SimulatorType.NumPy: np.random.default_rng(seed=self.__seed),
            SimulatorType.CuPy: None
        }
        return __switcher[self.__type]

    def __init__(self, seed: int, simulator_type: SimulatorType):
        self.__seed = seed
        self.__type = simulator_type
        self.rand_gen = self._switch_rand_gen()

    def bivariate_norm(self, correlation: float, no_draws: (int, int)):
        if correlation > 1 or correlation < -1:
            raise ValueError("Wrong correlation")
        __switcher = {
            SimulatorType.NotSet: None,
            SimulatorType.CppPone: None,
            SimulatorType.NumPy: self.rand_gen.multivariate_normal(
                mean=[0, 0],
                cov=[[1, correlation], [correlation, 1]],
                size=no_draws,
                check_valid="ignore"),
            SimulatorType.CuPy: None
        }
        return __switcher[self.__type]



