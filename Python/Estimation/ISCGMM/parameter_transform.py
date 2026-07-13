from __future__ import annotations

import math

import numpy as np

from Models.Heston.parameters import HestonParameters


def to_free(theta: HestonParameters) -> np.ndarray:
    theta.validate()
    return np.array(
        [theta.eta, math.log(theta.kappa), math.log(theta.vbar), math.log(theta.sigma_v), np.arctanh(theta.rho), math.log(theta.kappa_q)],
        dtype=float,
    )


def from_free(free_theta: np.ndarray, *, r: float = 0.0, q: float = 0.0) -> HestonParameters:
    free = np.asarray(free_theta, dtype=float)
    if free.shape != (6,):
        raise ValueError("free_theta must have shape (6,)")
    eta, kappa, vbar, sigma_v, rho, kappa_q = (
        float(free[0]),
        float(np.exp(free[1])),
        float(np.exp(free[2])),
        float(np.exp(free[3])),
        float(np.tanh(free[4])),
        float(np.exp(free[5])),
    )
    theta = HestonParameters(
        eta=eta,
        kappa=kappa,
        vbar=vbar,
        sigma_v=sigma_v,
        rho=rho,
        eta_v=kappa - kappa_q,
        r=float(r),
        q=float(q),
    )
    theta.validate()
    return theta


def free_parameter_bounds(natural_bounds: tuple[tuple[float, float], ...]) -> tuple[tuple[float, float], ...]:
    """Convert bounds for eta, kappa, vbar, sigma_v, rho and kappa_q."""

    if len(natural_bounds) != 6:
        raise ValueError("natural_bounds must contain six pairs")
    out: list[tuple[float, float]] = []
    for index, (lower, upper) in enumerate(natural_bounds):
        if lower >= upper:
            raise ValueError("natural parameter lower bounds must be below upper bounds")
        if index in {1, 2, 3, 5}:
            if lower <= 0.0:
                raise ValueError("positive natural parameters require strictly positive lower bounds")
            out.append((math.log(lower), math.log(upper)))
        elif index == 4:
            if not -1.0 < lower < upper < 1.0:
                raise ValueError("rho bounds must lie strictly inside (-1, 1)")
            out.append((float(np.arctanh(lower)), float(np.arctanh(upper))))
        else:
            out.append((float(lower), float(upper)))
    return tuple(out)

