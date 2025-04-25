from advection import *
from ncviewer import NcFile
from initial_conditions import *

x0 = 0
xf = 2
nx = 2000
T = 1
a = 1
cfl = 0.8

filepath = method_of_characteristics(x0, xf, nx, T, 125, a, heaviside(0.3, 0.7))

ncv = NcFile(filepath)
ncv.evolution(0)