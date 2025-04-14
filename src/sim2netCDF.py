# -*- coding: utf-8 -*-
from netCDF4 import Dataset
import time


class ncf:
    def __init__(self, full_path, coords, vars, description):
        """
        Initializes the nc object by creating a NetCDF file.
        
        Parameters:
        - full_path: Directory where the file will be saved and name of the file (with extension) Ex: simulation/sim.nc.
        - coords: Dictionary of coordinates {name: values}.
        - vars: List of variable names to store in the file.
        - description: String describing the simulation.
        """
        self.filepath = f"{full_path}"
        self.rootgrp = Dataset(self.filepath, "w", format="NETCDF4")

        # Dimensions
        self.rootgrp.createDimension("t", None)
        self.coords_names = []
        for coord_name, coord_values in coords.items():
            self.rootgrp.createDimension(coord_name, len(coord_values))
            coord_var = self.rootgrp.createVariable(coord_name, "f4", (coord_name,))
            coord_var[:] = coord_values
            coord_var.unit = "unit"  # Default unit
            self.coords_names.append(coord_name)

        # Variables
        self.rootgrp.createVariable("t", "f8", ("t",)).units = "s"
        for var_name in vars:
            var = self.rootgrp.createVariable(var_name, "f4", (*self.coords_names, "t"))
            var.units = "unit"

        # Attributes
        self.rootgrp.description = description
        self.rootgrp.history = "Created " + time.ctime(time.time())
        self.rootgrp.source = "FDM simulation"

    def save(self, current_time, vars):
        """
        Save simulation variables for the current iteration.

        Parameters:
        - current_time: Current time of the simulation (float).
        - vars: Dictionary {variable_name: numpy_array}.
        """
        time_dim = len(self.rootgrp.variables["t"])
        self.rootgrp.variables["t"][time_dim] = current_time

        for var_name, var_values in vars.items():
            if var_name in self.rootgrp.variables:
                self.rootgrp.variables[var_name][..., time_dim] = var_values
            else:
                raise ValueError(f"Variable '{var_name}' not found in the NetCDF file.")

    def close(self):
        """Close the NetCDF file."""
        self.rootgrp.close()
