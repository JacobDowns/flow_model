import numpy as np
from scipy.interpolate import UnivariateSpline
from dolfin import *

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

        # Flowline CG fields
        flowline_cg_fields = ['B', 'width', 'S_ref']
        flowline_cg_fields += ['P' + str(i) for i in range(12)]
        flowline_cg_fields += ['T' + str(i) for i in range(12)]
        # State CG fields
        state_cg_fields = ['H0_c']
        # State DG fields
        state_dg_fields = ['H0']
        # Flowline R fields
        flowline_r_fields = ['domain_len']
        # State R fields
        state_r_fields = ['L0', 't0']

        self.cg_fields = flowline_cg_fields + state_cg_fields
        self.dg_fields = state_dg_fields
        self.r_fields = flowline_cg_fields + state_cg_fields
        
        
        ### Create Fenics functions
        ########################################################

        self.input_functions = {}
        self.original_cg_functions = {}

        for field_name in self.cg_fields:
            self.input_functions[field_name] = Function(self.V_cg)
            self.original_cg_functions[field_name] = Function(self.V_cg)
        for field_name in self.dg_fields:
            self.input_functions[field_name] = Function(self.V_dg)
        for field_name in self.r_fields:
            self.input_functions[field_name] = Function(self.V_r)

    
        ### Load inputs from .h5 files
        ######################################################## 

        if 'input_file_name' in input_dict:
        
            for field_name in flowline_cg_fields:
                print(field_name)
                self.input_file.read(self.input_functions[field_name], field_name)
                #self.input_file.read(self.original_cg_functions[field_name], field_name)
                print("Done")
                print()

            """
            for field_name in state_cg_fields:
                print(field_name)
                self.state_input_file.read(self.input_functions[field_name], field_name)
                self.state_input_file.read(self.original_cg_functions[field_name], field_name)
                print("Done")
                print()

        
            for field in state_dg_fields:
                print 
                self.state_input_file.read(self.input_functions[field], field)
            for field in flowline_r_fields:
                self.input_file.read(self.input_functions[field], field)
            for field in state_r_fields:
                self.state_input_file.read(self.input_functions[field], field)
       
        #quit()        
        ### Create interpolated CG functions
        ########################################################

        self.interp_functions = {}
        
        if 'input_file_name' in input_dict:
            # Along flow coordinates
            x = self.mesh_coords
            for field in self.cg_fields:
                vals = project(self.input_functions[field]).compute_vertex_values()
                self.interp_functions[field] = UnivariateSpline(x, vals, k = 3, s =  0.005) 
        else:
            # Length of arrays
            N = input_dict['B']
            # Along flow coordinate
            x = np.linspace(0., 1., N)
            for field in self.cg_fields:                          
                self.interp_functions[field] = UnivariateSpline(x, input_dict[field], k = 3, s =  0.005)

                
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
            """
