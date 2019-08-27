import numpy as np
from scipy.interpolate import interp1d
from dolfin import *
from ice_model.ice_model import IceModel

class ModelWrapper(object):

    def __init__(self, inputs):

        self.inputs = inputs
        
        ### Model mesh
        ########################################################

        # Mesh
        self.mesh = inputs['mesh']
        # Get the mesh coordinates
        self.mesh_coords = self.mesh.coordinates()[:,0]
        # Normalize so that coordinates go from 0 to 1
        self.mesh_coords /= self.mesh_coords.max()
        # Along flow length coordinate for inputs
        self.x = inputs['x']
        # Domain length
        self.domain_len = self.x.max()
        # Initial glacier length
        self.L0 = inputs['L0']
        # Time step
        self.dt = inputs['dt']
    
            
        #### Function spaces
        ########################################################

        self.E_cg = FiniteElement('CG', self.mesh.ufl_cell(), 1)
        self.E_dg = FiniteElement('DG', self.mesh.ufl_cell(), 0)
        self.E_r = FiniteElement('R',  self.mesh.ufl_cell(), 0)
        self.V_cg = FunctionSpace(self.mesh, self.E_cg)
        self.V_dg = FunctionSpace(self.mesh, self.E_dg)
        self.V_r = FunctionSpace(self.mesh, self.E_r)
        
        self.input_functions = {}
        self.interp_functions = {}


        ### Create boundary facet function
        ########################################################################

        self.boundaries = MeshFunction('size_t', self.mesh, self.mesh.topology().dim() - 1, 0)

        for f in facets(self.mesh):
            if near(f.midpoint().x(), 1):
                # Terminus
               self.boundaries[f] = 1
            if near(f.midpoint().x(), 0):
                # Divide
               self.boundaries[f] = 2

               
        ### Ice model 
        ########################################################

        # Create ice model
        ice_params = inputs['ice_params']
        self.ice_model = IceModel(self, ice_params)


        ### Make sure everything is initialized
        ########################################################
        
        self.update_inputs(self.L0, {})
        


    # Load fields
    def load_fields(self, inputs, fields):
        for field_name in fields:
            self.input_functions[field_name] = Function(self.V_cg)
            self.interp_functions[field_name] = interp1d(self.x, inputs[field_name], kind = 'quadratic')

            
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


    # Assign model input functions
    def update_inputs(self, L, params):
 
        # Update length dependent fields
        self.update_interp_all(L)

        # Ice model parameters
        ice_params = {}
        if 'ice_params' in params:
            ice_params = params['ice_params']
                

        ### Update ice model inputs
        ########################################################
        
        # Update model inputs like sea level
        self.model.update(ice_params)
