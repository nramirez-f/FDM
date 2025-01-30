# -*- coding: utf-8 -*-
import numpy as np
from advection import characteristics
from netCDFviewer import ncv

def initial_condition(x):
    """Initial condition for the Cauchy problem"""
    z = np.zeros(len(x))
    z[x < 0] = 1.0
    z[(0 <= x) & (x <= 2)] = 1 - 1 / 8 * x[(0 <= x) & (x <= 2)]**2 * (3 - x[(0 <= x) & (x <= 2)])
    z[x > 2] = 1 / 2
    return z

filepath = characteristics(-1, 3, 100, 5, 10, 1/2, initial_condition)

nc = ncv(filepath)
nc.playShape(0)
nc.close()

