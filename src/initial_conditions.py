# -*- coding: utf-8 -*-
import numpy as np

def _prepare_input(x):
    """Ensure input is array and track if original was scalar."""
    isscalar = np.isscalar(x)
    x = np.atleast_1d(x)
    return x, isscalar

def _finalize_output(f, isscalar):
    """Return scalar if original input was scalar."""
    return f.item() if isscalar else f

def f1(x):
    """
    """
    x, isscalar = _prepare_input(x)
    f = np.zeros(len(x))
    f[x < 0] = 1.0
    xp_mask = (0 <= x) & (x <= 2)
    xp = x[xp_mask]
    f[xp_mask] = 1 - 0.125 * xp**2 * (3 - xp)
    f[x > 2] = 0.5
    return _finalize_output(f, isscalar)

def f2(x):
    """
    """
    x, isscalar = _prepare_input(x)
    f = np.zeros(len(x))
    f[x < 0] = 1.0
    xp_mask = (0 <= x) & (x <= 2)
    xp = x[xp_mask]
    f[xp_mask] = ((xp + 1) * (xp - 2)**2) * 0.25
    f[x > 2] = 0
    return _finalize_output(f, isscalar)

def f3(x):
    """
    """
    x, isscalar = _prepare_input(x)
    f = np.zeros(len(x))
    f[x < -1] = 0.5
    xp_mask = (-1 <= x) & (x <= 1)
    xp = x[xp_mask]
    f[xp_mask] = 0.5 - ((xp - 2) * (xp + 1)**2) * 0.125
    f[x > 1] = 1
    return _finalize_output(f, isscalar)

def f4(x):
    """
    """
    x, isscalar = _prepare_input(x)
    f = np.zeros(len(x))
    f[x < 0] = 1.1
    xp_mask = (0 <= x) & (x <= 2)
    xp = x[xp_mask]
    f[xp_mask] = 0.05 * xp**3 - 0.15 * xp**2 + 1.1
    f[x > 2] = 0.9
    return _finalize_output(f, isscalar)

def riemann(ul=1, ur=0, sp=0):
    """
    Riemann initial value problem
    """
    def shock(x):
        x, isscalar = _prepare_input(x)
        f = np.where(x < sp, ul, ur)
        return _finalize_output(f, isscalar)
    return shock

# From Riemann Solvers and Numerical Methods for Flid Dynamics by E. Toro
def bumping(alpha, beta):
    """
    """
    def bump(x):
        x, isscalar = _prepare_input(x)
        f = alpha*np.exp(-beta*x*x)
        return _finalize_output(f, isscalar)
    
    return bump

def heaviside(a, b, bottom=0, height=1):
    """
    """
    if (b < a):
        def step(x):
            x, isscalar = _prepare_input(x)
            f = np.zeros(len(x))
            return _finalize_output(f, isscalar)
    else:
        def step(x):
            x, isscalar = _prepare_input(x)
            f = np.zeros(len(x)) + bottom
            xp_mask = (a <= x) & (x <= b)
            f[xp_mask] = height
            return _finalize_output(f, isscalar)
    return step

