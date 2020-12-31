# Testing the heston ccf class
from ModelCharacteristicFunctions.ccf import Heston
import numpy as np


def main():
	heston_test = Heston()
	test_argument = np.array((3, 1), dtype=complex)
	test_maturities = np.arange(.1, 1, .1)
	eta = 2
	k = 5
	vbar = .1
	sigma = .1
	rho = -.5
	eta_v = 0
	test_parameters = np.array([eta, k, vbar, sigma, rho, eta_v])
	result =  heston_test.ccf_function(test_argument, test_maturities, test_parameters)
	print(result)
	return 1

main()
