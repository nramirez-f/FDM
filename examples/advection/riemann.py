from advection import *
from ncviewer import NcFile
from initial_conditions import *

x0 = -0.5
xf = 1
nx = 600
T = 1
a = 1
cfl = 0.5

filepath = cir(x0, xf, nx, T, cfl, a, riemann(1,0), sns=2)

ncv = NcFile(filepath)
ncv.evolution(0)