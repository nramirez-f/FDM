from fdm.advection import *
from ncviewer import NcView
from initial_conditions import *

x0 = 0
xf = 2
nx = 2000
T = 1
a = 1
cfl = 0.8

filepath = lax_friedrichs(x0, xf, nx, T, cfl, a, heaviside(0.3, 0.7), sns=10)

ncv = NcView(filepath)
ncv.evolution(0)