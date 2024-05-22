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

    def cos_transform_pricer():
        a = -5.0  # integration lower bound
        b = 5.0  # integration upper bound
        c = .0  # other integration bound
        d = b  # integration bound
        n_cos = 1024  # Integration bound	for COS pricing part
        k = np.linspace(0, n_cos - 1, n_cos)  # integration grid values
        # Routine set - up:
        arg_dma = k * np.pi * (d - a) / (b - a)
        arg_cma = k * np.pi * (c - a) / (b - a)
        fact_sin = k * np.pi / (b - a)
        chi = np.cos(arg_dma) * np.exp(d) - np.cos(arg_cma) * np.exp(c)
        chi = chi + np.multiply(fact_sin, np.sin(arg_dma)) * np.exp(d)
        chi = chi - np.multiply(fact_sin, np.sin(arg_cma)) * np.exp(c)
        test = np.power((1 + np.square(fact_sin)), (-1))
        chi = np.multiply(chi, test)
        phi = np.multiply((np.sin(arg_dma) - np.sin(arg_cma)),
                          (np.append(np.multiply((b - a), np.power((k[1:] * np.pi), (-1))), 0) +
                           np.append(np.zeros((n_cos - 1, 1)), d - c)))
        v_call = 2.0 / (b - a) * (chi - phi)
        return v_call  # assume okay for now

    return out_tensor
