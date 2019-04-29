import numpy as np
from scipy.interpolate import UnivariateSpline
from dolfin import *
import matplotlib.pyplot as plt

class ModelWrapper(object):

    def __init__(self, input_dict):

        ### Model mesh
        ########################################################

        if 'input_file_name' in input_dict:
            self.mesh = Mesh()
            input_file_name = input_dict['input_file_name']
            self.input_file = HDF5File(self.mesh.mpi_comm(), input_dict['input_file_name'], "r")
            self.state_input_file = HDF5File(self.mesh.mpi_comm(), input_dict['state_file_name'], "r")
            self.input_file.read(self.mesh, "/mesh", False)
        else :
            # Mesh resolution
            mesh_res = 1000.
            self.mesh = IntervalMesh(int(input_dict['domain_len'] / mesh_res), 0., 1.) 

        # Get the mesh coordinates
        self.mesh_coords = self.mesh.coordinates()[:,0]
        # Normalize so that coordinates go from 0 to 1
        self.mesh_coords /= self.mesh_coords.max()
    
            
        #### Function spaces
        ########################################################

        self.E_cg = FiniteElement('CG', self.mesh.ufl_cell(), 1)
        self.E_dg = FiniteElement('DG', self.mesh.ufl_cell(), 0)
        self.E_r = FiniteElement('R',  self.mesh.ufl_cell(), 0)
        self.V_cg = FunctionSpace(self.mesh, self.E_cg)
        self.V_dg = FunctionSpace(self.mesh, self.E_dg)
        self.V_r = FunctionSpace(self.mesh, self.E_r)

        # Dictionary of function spaces
        function_spaces = {
            'CG' : self.V_cg,
            'DG' : self.V_dg,
            'R'  : self.V_r
        }

        # Flowline functions and function spaces
        flowline_fields = {
            'B' : 'CG',
            'width' : 'CG',
            'S_ref' : 'CG',
            'domain_len' : 'R'
        }

        for i in range(12):
            flowline_fields['P' + str(i)] = 'CG'
            flowline_fields['T' + str(i)] = 'CG'

        # State functions and function spaces
        state_fields = {
            'H0_c' : 'CG',
            'H0' : 'DG',
            'L0' : 'R',
            't0' : 'R'
        }

        all_fields = flowline_fields.copy()
        all_fields.update(state_fields)
        
        
        ### Create Fenics functions
        ########################################################

        self.input_functions = {}
        self.original_cg_functions = {}
            
        for field_name, space in flowline_fields.items():
            function_space = function_spaces[space]
            self.input_functions[field_name] = Function(function_space)
            if space == 'CG':
                self.original_cg_functions[field_name] = Function(function_space)

        for field_name, space in state_fields.items():
            function_space = function_spaces[space]
            self.input_functions[field_name] = Function(function_space)
            if space == 'CG':
                self.original_cg_functions[field_name] = Function(function_space)


        ### Load inputs from .h5 files
        ######################################################## 

        if 'input_file_name' in input_dict:

            for field_name, space in flowline_fields.items():
                function_space = function_spaces[space]
                f = Function(function_space)
                self.input_file.read(f, field_name)
                self.input_functions[field_name].assign(f)

            for field_name, space in state_fields.items():
                function_space = function_spaces[space]
                f = Function(function_space)
                self.state_input_file.read(f, field_name)
                self.input_functions[field_name].assign(f)

            for field_name, space in all_fields.items():
                if space == 'CG':
                    self.original_cg_functions[field_name].assign(self.input_functions[field_name])

        
        ### Interpolate all the CG functions
        ########################################################

        self.interp_functions = {}
        
        if 'input_file_name' in input_dict:
            # Along flow coordinates
            x = self.mesh_coords
            # Create interpolated functions
            for field_name in self.original_cg_functions:
                vals = project(self.input_functions[field_name]).compute_vertex_values()
                self.interp_functions[field_name] = UnivariateSpline(x, vals, k = 3, s =  0.005)
        else:
            # Length of arrays
            N = len(input_dict['B'])
            # Along flow coordinate
            x = np.linspace(0., 1., N)
            # Create interpolated functions
            for field_name in self.original_cg_functions:
                self.interp_functions[field_name] = UnivariateSpline(x, input_dict[field_name], k = 3, s =  0.005)


        ### Initial state
        ########################################################################

        if 'input_file_name' in input_dict:
            self.domain_len = float(self.input_functions['L0'])
            L0 = float(self.input_functions['L0'])
            
        else:
            # Domain length
            self.domain_len = input_dict['domain_len']
            # Get the margin position
            H_n = self.interp_functions['H0_c'](self.mesh_coords)
            first_index = np.where(abs(H_n) > 15.)[0].max()
            x = np.linspace(0., 1., N)
            L0 = x[first_index] * self.domain_len
            # Update length
            self.update_interp_all(L0)
            # Set DG initial thickness
            self.input_functions['H0'] = interpolate(self.input_functions['H0_c'], self.V_dg)            
            
        dolfin.plot(self.input_functions['P10'])
        plt.show()
            
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

        for field_name in self.original_cg_functions:
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
