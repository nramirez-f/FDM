from fdm.advection import *
from ncviewer import NcFile
from initial_conditions import *

x0 = -1
xf = 10
nx = 1100
T = 10
a = 1
cfl = 0.8

filepath = lax_wendroff(x0, xf, nx, T, cfl, a, bumping(1, 8), sns=1)

ncv = NcFile(filepath)
ncv.evolution(0)