import numpy as np
from scipy.interpolate import interp1d
from dolfin import *
from ice_model.ice_model import IceModel
from hydro_model.hydro_model import HydroModel
from smb_model.smb_model import SMBModel
import matplotlib.pyplot as plt
import copy

set_log_level(40)

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
            if near(f.midpoint().x(), 0):
                # Terminus
               self.boundaries[f] = 2
            if near(f.midpoint().x(), 1):
                # Divide
               self.boundaries[f] = 1

               
        ### Ice model 
        ########################################################

        # Create ice model
        ice_params = inputs['ice_params']
        self.ice_model = IceModel(self, ice_params)


        ### Hydrology model 
        ########################################################

        # Create ice model
        hydro_params = {}
        if 'hydro_params' in inputs:
            hydro_params = inputs['ice_params']
    
        self.hydro_model = HydroModel(self, hydro_params)


        ### SMB model 
        ########################################################

        # Create smb model
        hydro_params = {}
        if 'smb_params' in inputs:
            smb_params = inputs['smb_params']
    
        self.smb_model = SMBModel(self, smb_params)


        ### Make sure everything is initialized
        ########################################################

        self.update_inputs(self.L0, {})


    # Load fields
    def load_fields(self, inputs, fields):
        for field_name in fields:
            self.input_functions[field_name] = Function(self.V_cg)
            self.interp_functions[field_name] = interp1d(self.x / self.x.max(), inputs[field_name], kind = 'quadratic')

        self.update_interp_fields(fields, self.L0)

        
    # Update only the given fields
    def update_interp_fields(self, field_names, L):
        frac = L / self.domain_len

        for field_name in field_names:
            self.input_functions[field_name].vector()[:] = \
             np.ascontiguousarray(self.interp_functions[field_name](self.mesh_coords * frac)[::-1])

            
    # Update all fields
    def update_interp_all(self, L):
        self.update_interp_fields(self.input_functions.keys(), L)
        

    # Get value of interpolated field at a point
    def get_interp_value(self, field_name, x):
        frac = x / self.domain_len
        return self.interp_functions[field_name](frac)


    # Assign model input functions
    def update_inputs(self, L, params):

 
        # Update length dependent fields
        self.update_interp_all(L)

        ### Update ice model inputs
        ########################################################

        # Ice model parameters
        ice_params = {}
        if 'ice_params' in params:
            ice_params = params['ice_params']
        # Update model inputs like sea level
        self.ice_model.update(ice_params)

        ### Update hydrology model inputs
        ########################################################

        # Ice model parameters
        hydro_params = {}
        if 'hydro_params' in params:
            hydro_params = params['ice_params']
        self.hydro_model.update(hydro_params)

        
        ### Update SMB model inputs
        ########################################################

        # SMB model parameters
        smb_params = {}
        if 'smb_params' in params:
            smb_params = params['smb_params']
        self.smb_model.update(smb_params)


    # Step the model forward by one time step
    def step(self, params = {}, accept = True):
        # Update the model inputs
        self.update_inputs(float(self.ice_model.L0), params)
        L = self.ice_model.step(accept)
        return L


    # Get the model state so that it can be saved
    def get_state(self):

        # Pop the mesh, otherwise we can't deep copy the inputs dictionary
        mesh = self.inputs.pop('mesh', None)
        state = copy.deepcopy(self.inputs)
        self.inputs['mesh'] = mesh
        
        # Interpolate the thickness
        H_interp = interp1d(self.mesh_coords*float(self.ice_model.L0), self.ice_model.H0_c.compute_vertex_values())
        H = np.zeros_like(state['x'])
        
        # Evaluate the thickness on the original input grid
        indexes = state['x'] < float(self.ice_model.L0)
        H[indexes] = H_interp(state['x'][indexes])
        state['ice_params']['H'] = H

        # Current time
        state['t0'] = self.ice_model.t
        state['L0'] = float(self.ice_model.L0)
        state.pop('dt', None)
        
        return state
        
