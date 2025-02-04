# -*- coding: utf-8 -*-
import numpy as np

def f1(x):
    f = np.zeros(len(x))
    f[x < 0] = 1.0
    xp = x[(0 <= x) & (x <= 2)]
    f[(0 <= x) & (x <= 2)] = 1 - 0.125 * pow(xp, 2) * (3 - xp)
    f[x > 2] = 0.5
    return f

def f2(x):
    f = np.zeros(len(x))
    f[x < 0] = 1.0
    xp = x[(0 <= x) & (x <= 2)]
    f[(0 <= x) & (x <= 2)] = ((xp + 1) * pow((xp - 2), 2)) * 0.25
    f[x > 2] = 0
    return f

def f3(x):
    f = np.zeros(len(x))
    f[x < -1] = 0.5
    xp = x[(-1 <= x) & (x <= 1)]
    f[(-1 <= x) & (x <= 1)] = 0.5 - ((xp - 2) * pow((xp + 1), 2)) * 0.125
    f[x > 1] = 1
    return f

def f4(x):
    f = np.zeros(len(x))
    f[x < 0] = 1.1
    xp = x[(0 <= x) & (x <= 2)]
    f[(0 <= x) & (x <= 2)] = 0.05 * pow(xp, 3) - 0.15 * pow(xp, 2) + 1.1
    f[x > 2] = 0.9
    return f