"""Thin wrapper around scipy.optimize.least_squares with named, bounded, fixable parameters."""

import os

import numpy as np
from scipy.optimize import least_squares

METHODS = ("trf", "dogbox")


class Param:
    def __init__(self, value, vary=True, min=-np.inf, max=np.inf):
        self.value = value
        self.vary = vary
        self.min = min
        self.max = max


def minimize(func, params, method="trf"):
    """Minimize the sum of squares of func(values) over the varying parameters.

    func receives a dict name -> value (fixed and varying) and returns a residual array.
    Returns (values dict, report string).
    """
    names = list(params)
    free = [n for n in names if params[n].vary]
    fixed = {n: params[n].value for n in names if not params[n].vary}

    def expand(x):
        values = dict(fixed)
        values.update(zip(free, x))
        return {n: values[n] for n in names}

    def fun(x):
        return np.asarray(func(expand(x)), dtype=float)

    if not free:
        values = expand([])
        r = fun([])
        nfev, status, message = 0, 0, "No varying parameters"
        cost = 0.5 * float(np.sum(r ** 2))
        stderr = {}
    else:
        lb = np.array([params[n].min for n in free], dtype=float)
        ub = np.array([params[n].max for n in free], dtype=float)
        x0 = np.clip([params[n].value for n in free], lb, ub)
        res = least_squares(fun, x0, bounds=(lb, ub), method=method)
        values = expand(res.x)
        r, cost, nfev, status, message = res.fun, res.cost, res.nfev, res.status, res.message
        stderr = _std_errors(res, free)

    lines = [
        "[[Fit Statistics]]",
        f"    method              = {method}",
        f"    status              = {status} ({message})",
        f"    function evals      = {nfev}",
        f"    number of residuals = {len(r)}",
        f"    chi-square          = {2 * cost:.6g}",
        "[[Variables]]",
    ]
    for n in names:
        p = params[n]
        if not p.vary:
            lines.append(f"    {n}: {values[n]:.6g} (fixed)")
        else:
            err = f" +/- {stderr[n]:.3g}" if n in stderr else ""
            lines.append(f"    {n}: {values[n]:.6g}{err} (init = {p.value:.6g}, bounds = [{p.min}:{p.max}])")
    return values, os.linesep.join(lines)


def _std_errors(res, free):
    m, n = res.jac.shape
    if m <= n:
        return {}
    try:
        cov = np.linalg.inv(res.jac.T @ res.jac) * (2 * res.cost / (m - n))
        return dict(zip(free, np.sqrt(np.diag(cov))))
    except np.linalg.LinAlgError:
        return {}
