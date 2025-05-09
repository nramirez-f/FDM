import numpy as np
from ncfiles import NcFile
from scipy.sparse import diags

def method_of_characteristics(x0:float, xf:float, nx:int, T:float, nt:int, a:float, f, sns:int = 1, path_to_save="simulations"):
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
    
    a : float
        The wave speed (advection velocity).
    
    f : function
        The initial condition function, which defines the initial profile of the solution.

    sns : int
        snapshot step to save simulation.
    
    path_to_save : str
        path to save the simulation. Name will be advection-exact

    Returns:
    --------
    str
        The file path where the simulation results are saved.
    """

    x = np.linspace(x0, xf, nx)
    dt = T / nt

    full_path = f"advection-exact.nc"
    ncf = NcFile(full_path, title='Advection simulation by Method of Characteristics', author = 'Nramirez', institution = 'EDANYA Group (University of Malaga)', source = 'https://github.com/nramirez-f/Finite-Differences', references ='LeVeque, Randall J.: Numerical Methods for Conservation Laws 1992')
    ncf.addCoords({"x": x})
    ncf.addVars(['u'])

    for k in range(nt + 1):
        t = dt * k
        u = f(x - a * t)
        if (a > 0):
            u[0] = 0
        if (k % sns == 0):
            ncf.save(t, {"u": u})

    return full_path

def _matrix(a:float, nu, dim, method:str):

    if (method == "cir"):
        if (a > 0):
            diag = (1 - nu) * np.ones(dim)
            lower_diag = nu * np.ones(dim-1)
            A = diags([diag, lower_diag], [0, -1], shape=(dim, dim), format='csr')
        else:
            diag = (1 + nu) * np.ones(dim)
            upper_diag = (-nu) * np.ones(dim-1)
            A = diags([diag, upper_diag], [0, 1], shape=(dim, dim), format='csr')

    elif (method == "lax_friedichs"):
        lower_diag = 0.5 * (1+nu) * np.ones(dim-1)
        upper_diag = 0.5 * (1-nu) * np.ones(dim-1)
        A = diags([lower_diag, upper_diag], [-1, 1], shape=(dim, dim), format='csr')

    elif (method == "lax_wendroff"):
        lower_diag = 0.5 * nu * (nu+1) * np.ones(dim-1)
        diag = (1-nu*nu) * np.ones(dim)
        upper_diag = 0.5 * nu * (nu-1) * np.ones(dim-1)
        A = diags([lower_diag, diag, upper_diag], [-1, 0, 1], shape=(dim, dim), format='csr')

    elif (method == "beam_warming"):
        if (a > 0):
            lower_lower_diag = 0.5 * nu * (nu-1) * np.ones(dim-2)
            lower_diag = nu * (2-nu) * np.ones(dim-1)
            diag = 0.5 * (2-3*nu+nu*nu) * np.ones(dim)
            A = diags([lower_lower_diag, lower_diag, diag], [-2, -1, 0], shape=(dim, dim), format='csr')
        else:
            upper_upper_diag = 0.5 * (-nu) * ((-nu)-1) * np.ones(dim-2)
            upper_diag = (-nu) * (2+nu) * np.ones(dim-1)
            diag = 0.5 * (2+3*nu+nu*nu) * np.ones(dim)
            A = diags([diag, upper_diag, upper_upper_diag], [0, 1, 2], shape=(dim, dim), format='csr')

    elif (method == "fromm"):
        if (a > 0):
            lower_lower_diag = (-0.25) * (1-nu) * nu * np.ones(dim-2)
            lower_diag = 0.25 * (5-nu) * nu * np.ones(dim-1)
            diag = 0.25 * (1-nu) * (4+nu) * np.ones(dim)
            upper_diag = (-0.25) * (1-nu) *np.ones(dim-1)
            A = diags([lower_lower_diag, lower_diag, diag, upper_diag], [-2, -1, 0, 1], shape=(dim, dim), format='csr')

    else:
        raise RuntimeError("404 - Method Not Found")
    

    # Dirichlet Conditions
    A = A.tolil()

    A[0, :] = 0
    A[0, 0] = 1
    if (method == "beam_warming" or method == "fromm"):
        if (a > 0):
            A[1, :] = 0
            A[1, 1] = 1
        else:
            A[-2, :] = 0
            A[-2, -2] = 1

    A[-1, :] = 0
    A[-1, -1] = 1

    A = A.tocsr() 


    return A

def _boundary_conditions(u, x0, xf, f, type:str, a:float, method:str):

    # Dirichlet
    if (type == "dirichlet"):
            u[0] = f(x0)
            u[-1] = f(xf)

            if (method == "beam_warming"  or method == "fromm"):
                if (a > 0):
                    u[1] = u[0]
                else:
                    u[-2] = u[-1]
        
    return u

def _iteration(a, nu, dim, method_name, u, iteration_type="iterative"):

    if (a == 0):
        u = u
    else:
        if (iteration_type == "iterative"):
            u = _matrix(a, nu, dim, method_name) @ u 
        else:
            u

    return u

def _one_step_method(method_name:str, x0:float, xf:float, nx:int, T:float, cfl:float, a:float, f, t0:float = 0, sns:int = 1,  path_to_save="simulations"):
    """
    Solves the advection equation using the cir method and saves the results.

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

    cfl : float
        Courant-Friedichs-Levy condition of the method.
        
        Stability Condition:
            0<= cfl <= 1
    
    nt : int
        The number of time steps.
    
    a : float
        The wave speed (advection velocity).
    
    f : function
        The initial condition function, which defines the initial profile of the solution.

    t0 : float
        initial simulation time.

    sns : int
        Snapshot step to save the simulation.
    
    path_to_save : str
        path to save the simulation. Name will be advection-<method>

    Returns:
    --------
    str
        The file path where the simulation results are saved.
    """

    if (x0 >= xf):
        raise RuntimeError("Imposible Domain - xf must be greater than x0")
    
    N = nx - 2
    x = np.linspace(x0, xf, nx)
    dx = np.abs(xf - x0) / nx
    if (a == 0):
        dt = (T - t0)
    else:
        dt = (cfl * dx) / np.abs(a) 

    # Courant Number (positive)
    nu = a * dt / dx

    # Info
    info = f"Model: Advection\nMethod: {method_name}\nDimension: 1D\nMesh: [{x0}, {xf}]\ndx: {dx}\nInterval Time: [{t0}, {T}]\ndt: {dt}\nCFL: {cfl}\nCourant Number: {nu}\n"

    # Initial Condition
    u0 = f(x)

    if (not (0 <= cfl and  cfl <= 1)):
        raise RuntimeError("Unstable method - CFL condition not satisfied")

    full_path = f"advection-{method_name}.nc"
    ncf = NcFile(full_path, title=f'Advection simulation by Method {method_name}', description=info, author = 'Nramirez', institution = 'EDANYA Group (University of Malaga)', source = 'https://github.com/nramirez-f/Finite-Differences', references ='LeVeque, Randall J.: Numerical Methods for Conservation Laws 1992')
    ncf.addCoords({'x': x})
    ncf.addVars(['u'])

    # Save initial condition
    ncf.save(t0, {"u": u0})

    u = u0.copy()
    t = t0
    k = 1
    ks = 1
    while t < T:
        t = t0 + dt * k

        u = _iteration(a, nu, N+2, method_name, u)

        # Boundary conditions
        u = _boundary_conditions(u, x0, xf, f, 'dirichlet', a, method_name)
        
        # Snapshot of simulation
        if (k % sns == 0):
            ncf.save(t, {"u": u})
            ks+=1

        k+=1

    info += f"Total iterations: {k-1}\nIterations saved: {ks-1}\n"

    print(info)

    return full_path

def select_method(method_name):
    def method(x0: float, xf: float, nx: int, T: float, cfl: float, a: float, f,
               t0: float = 0, sns: int = 1, path_to_save: str = "simulations") -> str:
        """
        Solves the advection equation using the method and saves the results.

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

        cfl : float
            Courant-Friedrichs-Levy condition of the method.

            Stability Condition:
                0<= cfl <= 1

        a : float
            The wave speed (advection velocity).

        f : function
            The initial condition function, which defines the initial profile of the solution.

        t0 : float
            initial simulation time.

        sns : int
            Snapshot step to save the simulation.

        path_to_save : str
            Path to save the simulation. Name will be advection-<method>.

        Returns:
        --------
        str
            The file path where the simulation results are saved.
        """
        return _one_step_method(method_name, x0, xf, nx, T, cfl, a, f, t0, sns, path_to_save)
    
    method.__name__ = method_name
    return method

# One Step Methods
cir = select_method("cir")
lax_friedrichs = select_method("lax_friedichs")
lax_wendroff = select_method("lax_wendroff")
beam_warming = select_method("beam_warming")
fromm = select_method("fromm")
