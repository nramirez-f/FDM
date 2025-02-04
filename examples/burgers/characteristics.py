from burgers import plot_method_of_characteristics as pmoc
from initial_conditions import *
import numpy as np

x0 = -6
xf = 6
nx = 400
T = 10
nt = 40
c = -0.5

pmoc(x0, xf, nx, T, nt, f1)
