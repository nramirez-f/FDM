# -*- coding: utf-8 -*-
import netCDF4 as nc
import time
import os

class NcFile:
    def __init__(self, full_path, title='', description = '', author = '', institution = '', source = '', references ='', format_file="NETCDF4"):
        """
        Initializes the NcFile object by creating a NetCDF file.

        """
        os.makedirs(os.path.dirname(full_path), exist_ok=True)

        self.filepath = full_path
        self.ffile = format_file
        ncf = nc.Dataset(self.filepath, 'w', self.ffile)

        # Attributes
        ncf.history = "Created " + time.ctime(time.time())

        if title:
            ncf.title = title
        if description:
            ncf.description = description
        if author:
            ncf.author = author
        if institution:
            ncf.institution = institution
        if source:
            ncf.source = source
        if references:
            ncf.references = references

        ncf.close()

    def addCoords(self, spatialCoords, iterName='t', iterUnit = 's'):

        ncf = nc.Dataset(self.filepath, 'a', self.ffile)
        
        ncf.createDimension(iterName, None)
        ncf.createVariable(iterName, "f8", (iterName,)).units = iterUnit

        self.coords_names = []
        for coord_name, coord_values in spatialCoords.items():
            ncf.createDimension(coord_name, len(coord_values))
            coord_var = ncf.createVariable(coord_name, "f4", (coord_name,))
            coord_var[:] = coord_values
            coord_var.unit = "unit"
            self.coords_names.append(coord_name)

        ncf.close()

    def addVars(self, vars):
        
        ncf = nc.Dataset(self.filepath, 'a', self.ffile)

        for var_name in vars:
            var = ncf.createVariable(var_name, "f4", (*self.coords_names, "t"))
            var.units = "unit"

        ncf.close()


    def save(self, current_time, vars):
        """
        Save simulation variables for the current iteration.

        Parameters:
        - current_time: Current time of the simulation (float).
        - vars: Dictionary {variable_name: numpy_array}.
        """
        ncf = nc.Dataset(self.filepath, 'a', format="NETCDF4")

        time_dim = len(ncf.variables["t"])
        ncf.variables["t"][time_dim] = current_time

        for var_name, var_values in vars.items():
            if var_name in ncf.variables:
                ncf.variables[var_name][..., time_dim] = var_values
            else:
                raise ValueError(f"Variable '{var_name}' not found in the NetCDF file.")
            
        ncf.close()
