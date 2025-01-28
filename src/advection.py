import numpy as np
from sim2netCDF import ncf

def characteristics(x0, xf, nx, T, nt, c, f):
    """
    Animación de características de una ecuación diferencial.
    """
    # Configuración inicial
    x = np.linspace(x0, xf, nx)  # Dominio espacial
    dt = T / nt  # Paso de tiempo

    # Inicializar archivo NetCDF
    coords = {"x": x}
    vars = ["u"]
    full_path = "simulations/sim-advection-characteristic.nc"
    description = "Advection simulation by characteristic method"
    rootgrp = ncf(full_path, coords, vars, description)

    # Iterar sobre el tiempo y guardar datos
    for k in range(nt):
        current_time = dt * k
        u = f(x - c * current_time)  # Solución en el tiempo actual
        ncf.save(rootgrp, current_time, {"u": u})

    # Cerrar archivo NetCDF
    ncf.close(rootgrp)

    return full_path
