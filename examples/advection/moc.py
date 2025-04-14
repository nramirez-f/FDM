from advection import method_of_characteristics as moc
from ncviewer import NcFile
from initial_conditions import *

x0 = -4
xf = 4
nx = 400
T = 10
nt = 40
c = 0.5

filepath = moc(x0, xf, nx, T, nt, c, f3)

ncv = NcFile(filepath)
ncv.evolution(0)
