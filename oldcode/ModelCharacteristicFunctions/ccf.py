import numpy as np


class Ccf:
	pass


class Heston(Ccf):
	@staticmethod
	def ccf_function(argument: np.array, maturities: np.array, parameters: np.array) -> np.array:
		# Allocate	parameters
		# eta = parameters[0]
		k = parameters[1]
		vbar = parameters[2]
		sigmav = parameters[3]
		rho = parameters[4]
		eta_v = parameters[5]
		# ad-hoc, Heston	specific	analytic 	solution for the coefficients of the characteristic function
		vbar = k / (k + eta_v) * vbar
		k = k + eta_v
		bb = argument * sigmav * rho - k
		aa = argument * (1 - argument)
		gamma = np.sqrt(np.power(bb, 2) + aa * np.power(sigmav, 2))
		out = np.empty([argument.size, 3, maturities.size])
		test = np.outer((gamma + bb) / np.power(sigmav, 2), maturities)
		out[:][1][:] = -k * vbar * (np.outer((gamma + bb) / np.power(sigmav, 2), maturities) + 2 / np.power(sigmav, 2) * np.log(
			1 - np.multiply(np.divide((gamma + bb), (2 * gamma)), (1 - np.exp(-np.outer(gamma, maturities))))))
		out[:][2][:] = np.matlib.repmat(argument, 1, maturities.size)
		out[:][3][:] = - np.multiply(aa, np.divide((1 - np.exp(-gamma * maturities)), (
				2 * gamma - np.multiply((gamma + bb), (1 - np.exp(-gamma * maturities))))))
		return out
