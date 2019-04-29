import numpy as np
from scipy.interpolate import UnivariateSpline
from dolfin import *


class FlowModelWrapper(object):

    def __init__(self, input_dict):

        ### Model mesh
        ########################################################

        if 'input_file_name' in input_dict:
            self.mesh = Mesh()
            input_file_name = input_dict['input_file_name']
            self.input_file = HDF5File(self.mesh.mpi_comm(), input_file_name, "r")
            self.input_file.read(self.mesh, "/mesh", False)
        else :
            # Mesh resolution
            mesh_res = 1000.
            if 'mesh_res' in input_dict:
                mesh_res = input_dict['mesh_res']
            self.mesh = IntervalMesh(int(input_dict['domain_len'] / mesh_res), 0., 1.) 

        # Get the mesh coordinates
        self.mesh_coords = self.mesh.coordinates()[:,0]
        # Normalize so that coordinates go from 0 to 1
        self.mesh_coords /= self.mesh_coords.max()
        
            
        #### Function spaces
        ########################################################

        E_cg = FiniteElement('CG', mesh.ufl_cell(), 1)
        E_dg = FiniteElement('DG', mesh.ufl_cell(), 0)
        E_r = FiniteElement('R',  mesh.ufl_cell(), 0)
        self.V_cg = FunctionSpace(mesh, E_cg)
        self.V_dg = FunctionSpace(mesh, E_dg)
        self.V_r = FunctionSpace(mesh, E_r)

        # A list of CG fields
        self.cg_fields = ['B', 'width', 'S_ref', 'H0_c']
        self.cg_fields += ['P' + str(i) for i in range(12)]
        self.cg_fields += ['T' + str(i) for i in range(12)]
        # A list of DG fields
        self.dg_fields = ['H0']
        # A list of R fields
        self.r_fields = ['domain_len', 'L0']

        
        ### Create functions in the correct function space
        ########################################################

        self.input_functions = {}
        self.original_cg_functions = {}

        for field in self.cg_fields:
            self.input_functions[field] = Function(self.V_cg)
            self.original_cg_functions[field] = Function(self.V_cg)
        for field in self.dg_fields:
            self.input_functions[field] = Function(self.V_dg)
        for field in self.r_fields:
            self.input_functions[field] = Function(self.V_r)
    
        
        ### Load inputs from an .h5 file
        ######################################################## 

        if 'input_file_name' in input_dict:
            for field in self.cg_fields:
                self.input_file.read(self.input_functions[field], field)
                self.input_file.read(self.original_cg_functions[field], field)
            for field in self.dg_fields:
                self.input_file.read(self.input_functions[field], field)
            for field in self.r_fields:
                self.input_file.read(self.input_functions[field], field)
            

        ### Create interpolated CG functions
        ########################################################

        self.interp_functions = {}
        
        if 'input_file_name' in input_dict:
            for field in self.cg_fields:
                self.interp_functions[field_name] = UnivariateSpline(self.mesh_coords, project(self.input_functions[field]).compute_vertex_values(), k = 3, s =  0.005)
        else:
            # Domain length
            self.domain_len = input_dict['domain_len']
            # Length of arrays
            N = input_dict['B']
            # Along flow coordinate
            x = np.linspace(0., 1., N)
            for field in self.cg_fields:
                # Create interpolated functions                            
                self.interp_functions[field] = UnivariateSpline(x, input_dict[field], k = 3, s =  0.005)


    def __load_input_file__(self, input_file_name, state_file_name):

        input_functions = {}
        
        for field in self.cg_fields:
            input_functions_field = Function(V_cg
