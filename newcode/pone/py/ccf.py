import numpy as np


def cf_coeffs_heston(u, T_vect, etaV, k, vBar, sigmaV, rho):
    # Initialize output tensor
    num_u = len(u)
    num_T = len(T_vect)
    out_tensor = np.zeros((num_u, 3, num_T))

    # Ad-hoc, Heston specific analytic solution for the coefficients of the characteristic function
    vBar_star = k / (k - etaV) * vBar
    k_star = k - etaV
    bb = u * sigmaV * rho - k_star
    aa = u * (1 - u)
    gamma = np.sqrt(bb ** 2 + aa * sigmaV ** 2)

    out_tensor[:, 0, :] = -k_star * vBar_star * ((gamma + bb) / sigmaV ** 2 * T_vect + 2 / sigmaV ** 2 * np.log(
        1 - (gamma + bb) / (2 * gamma) * (1 - np.exp(-gamma * T_vect))))
    out_tensor[:, 1, :] = np.tile(u[:, np.newaxis], (1, num_T))
    out_tensor[:, 2, :] = -aa * (1 - np.exp(-gamma * T_vect)) / (
                2 * gamma - (gamma + bb) * (1 - np.exp(-gamma * T_vect)))

    return out_tensor
