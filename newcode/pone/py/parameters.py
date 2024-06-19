import numpy as np
from typing import List


class ParameterSet:
    def __init__(self, num_param: int, param: np.array,  num_p_param: int, p_param: np.array, num_q_param: int, q_param: np.array):
        if num_param != len(param):
            raise ValueError("Wrong param")
        self._num_param = num_param
        self._param = param
        if num_p_param != len(p_param):
            raise ValueError("Wrong P param")
        self._num_p_param = num_p_param
        self._p_param = p_param
        if num_q_param != len(q_param):
            raise  ValueError("Wrong Q param")
        self._num_q_param = num_q_param
        self._q_param = q_param

    @classmethod
    def from_list(cls, param: List[np.array]):
        __len = len(param)
        if __len < 1 or __len > 3:
            raise ValueError("Wrong param")
        if __len == 1:
            return cls(len(param[0]), param[0], 0, None, 0, None)
        if __len == 2:
            return cls(len(param[0]), param[0], len(param[1]), param[1], 0, None)
        if __len == 3:
            return cls(len(param[0]), param[0], len(param[1]), param[1], len(param[2]), param[2])

    def all_param_sets(self) -> List[np.array]:
        return [self.param, self.p_param, self.q_param]

    def param(self):
            return self._param

    def p_param(self) -> np.array:
        return self._p_param

    def q_param(self) -> np.array:
        return self._q_param

    def has_p_set(self) -> bool:
        return len(self.p_param()) > 0

    def has_q_set(self) -> bool:
        return len(self.q_param()) > 0

