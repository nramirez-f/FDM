import numpy as np
from sim2netCDF import ncf

def method_of_characteristics(x0, xf, nx, T, nt, c, f):
    """
    Solves the advection equation using the method of characteristics and saves the results.

    Parameters:
    -----------
    x0 : float
        The initial spatial coordinate (left boundary of the domain).
    
    xf : float
        The final spatial coordinate (right boundary of the domain).
    
    nx : int
        The number of spatial grid points.
    
    T : float
        The total simulation time.
    
    nt : int
        The number of time steps.
    
    c : float
        The wave speed (advection velocity).
    
    f : function
        The initial condition function, which defines the initial profile of the solution.

    Returns:
    --------
    str
        The file path where the simulation results are saved.
    """

    x = np.linspace(x0, xf, nx)
    dt = T / nt

    coords = {"x": x}
    vars = ["u"]
    full_path = "simulations/sim-advection-moc.nc"
    description = "Advection simulation by Method of Characteristic"
    rootgrp = ncf(full_path, coords, vars, description)

    for k in range(nt + 1):
        t = dt * k
        u = f(x - c * t)
        ncf.save(rootgrp, t, {"u": u})

    ncf.close(rootgrp)

    return full_path
