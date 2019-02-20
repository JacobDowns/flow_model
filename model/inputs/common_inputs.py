from dolfin import *
import numpy as np
from scipy.interpolate import UnivariateSpline
from model.forward_model.support.physical_constants import *

"""
Flexible class for generating a forward model inputs file.
"""

class CommonInputs(object):

    def __init__(self, input_file_name, input_options = {}):

        # Dictionary of extra input options
        self.input_options = input_options


        ### Mesh and function spaces
        ########################################################################

        # Load mesh
        self.mesh = Mesh()
        input_file = HDF5File(self.mesh.mpi_comm(), input_file_name, "r")
        self.input_file = input_file
        input_file.read(self.mesh, "/mesh", False)
        
        # Store mesh coordinates
        self.mesh_coords = self.mesh.coordinates()[:,0]
        # Define function spaces
        self.E_cg = FiniteElement('CG', self.mesh.ufl_cell(), 1)
        self.E_dg = FiniteElement('DG', self.mesh.ufl_cell(), 0)
        self.E_r = FiniteElement('R',  self.mesh.ufl_cell(), 0)
        self.V_cg = FunctionSpace(self.mesh, self.E_cg)
        self.V_dg = FunctionSpace(self.mesh, self.E_dg)
        self.V_r = FunctionSpace(self.mesh, self.E_r)
        # A dictionary of all model inputs
        self.input_functions = {}
        # Physical constants
        self.physical_constants = pcs


        ### CG inputs
        ########################################################################

        if 'cg_fields' in input_options:
            self.cg_fields = input_options['cg_fields']
            self.interp_fields = set(self.cg_fields).intersection(set(input_options['interp_fields']))
        else :
            # Default CG fields
            self.cg_fields = ['H0_c', 'B', 'width', 'beta2']

            # Any additional fields
            if 'additional_cg_fields' in input_options:
                self.cg_fields += input_options['additional_cg_fields']

            self.interp_fields = set(['B', 'width', 'beta2', 'S_ref'])
            if 'additional_interp_fields' in input_options:
                self.interp_fields = self.interp_fields.union(set(self.cg_fields).intersection(set(input_options['additional_interp_fields'])))

        # For interpolated CG functions, we want to store the original function as well
        self.original_cg_functions = {}

        # Load all CG inputs
        for field_name in self.cg_fields:
            self.original_cg_functions[field_name] = Function(self.V_cg)
            self.input_functions[field_name] = Function(self.V_cg)

            input_file.read(self.input_functions[field_name], field_name)
            input_file.read(self.original_cg_functions[field_name], field_name)

            if field_name == 'beta2' and 'beta2' in input_options:
                self.input_functions[field_name].interpolate(Constant(input_options['beta2']))


        ### DG inputs
        ########################################################################

        if 'dg_fields' in input_options:
            self.dg_fields = input_options['dg_fields']
        else :
            # Default DG fields
            self.dg_fields = ['H0']
            # Any additional fields
            if 'additional_dg_fields' in input_options:
                self.dg_fields += self.dg_fields['additional_dg_fields']

        # Load all DG inputs
        for field_name in self.dg_fields:
            self.input_functions[field_name] = Function(self.V_dg)
            input_file.read(self.input_functions[field_name], field_name)


        ### R inputs
        ########################################################################

        if 'r_fields' in input_options:
            # Domain len is non-optional
            self.r_fields = set(['domain_len'])
            self.r_fields = self.r_fields.union(set(input_options['r_fields']))
        else :
            # Default DG fields
            self.r_fields = ['L0', 'domain_len']
            # Any additional fields
            if 'additional_r_fields' in input_options:
                self.r_fields += self.dg_fields['additional_r_fields']

        # Load all R inputs
        for field_name in self.r_fields:
            self.input_functions[field_name] = Function(self.V_r)
            input_file.read(self.input_functions[field_name], field_name)

        # Store the domain length
        self.domain_len = float(self.input_functions['domain_len'])


        ### Create interpolated CG fields
        ########################################################################

        # Dictionary of interpolated functions
        self.interp_functions = {}
        # Get the mesh coordinates
        self.mesh_coords = self.mesh.coordinates()[:,0]
        # Normalize so that coordinates go from 0 to 1
        self.mesh_coords /= self.mesh_coords.max()

        # Create interpolated functions
        for field_name in self.interp_fields:
            skip_smooth = 1
            # Number of data points to skip between interpolated points - for smoothness
            self.interp_functions[field_name] = UnivariateSpline(self.mesh_coords, project(self.input_functions[field_name]).compute_vertex_values(), k = 3, s =  0.005)


        #### Create boundary facet function
        ########################################################################
        self.boundaries = MeshFunction('size_t', self.mesh, self.mesh.topology().dim() - 1, 0)

        for f in facets(self.mesh):
            if near(f.midpoint().x(), 1):
                # Terminus
               self.boundaries[f] = 1
            if near(f.midpoint().x(), 0):
               # Divide
               self.boundaries[f] = 2


    # Update all inputs that depend on glacier length L
    def update_interp_all(self, L):
        frac = L / self.domain_len

        for field_name in self.interp_fields:
            self.input_functions[field_name].vector()[:] = \
             np.ascontiguousarray(self.interp_functions[field_name](self.mesh_coords * frac)[::-1])


    # Update only the given fields
    def update_interp_fields(self, field_names, L):
        frac = L / self.domain_len

        for field_name in field_names:
            self.input_functions[field_name].vector()[:] = \
             np.ascontiguousarray(self.interp_functions[field_name](self.mesh_coords * frac)[::-1])


    # Get value of interpolated field at a point
    def get_interp_value(self, field_name, x):
        frac = x / self.domain_len
        return self.interp_functions[field_name](frac)
