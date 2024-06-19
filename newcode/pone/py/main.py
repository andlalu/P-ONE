# Sandbox to test the particle filtering
import numpy as np

from processes import Heston
from simulation import Simulator, SimulatorType

def test():
    hest = Heston(0, [np.log(100), np.sqrt(.3)], 5, 0.1, .3, -.5, 0, 0)
    simulator = Simulator(seed=1, simulator_type=SimulatorType.NumPy)
    test_res = hest.simulate(simulator, 10, 1/252, 5, 1)
    test = True

if __name__ == '__main__':
    test()

