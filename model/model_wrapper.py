import numpy as np
from dolfin import FiniteElement, FunctionSpace, set_log_level, \
    MeshFunction, near, Function, IntervalMesh, project, FacetNormal
import dolfin as df
from model.support.model_params import *
from model.momentum_model.momentum_model import MomentumModel
from model.mass_model.mass_model import MassModel
from model.hydro_model.hydro_model import HydroModel
import matplotlib.pyplot as plt

set_log_level(40)

class ModelWrapper(object):

    def __init__(self, inputs):

        self.inputs = inputs
        self.model_params = model_params
        self.fields = inputs['fields']
        
        ### Model mesh
        ##############################################################

        self.x = self.fields['x']
        self.mesh = IntervalMesh(len(self.x) - 1, self.x.min(), self.x.max())
        self.mesh.coordinates()[:] = self.x[:, np.newaxis]

            
        #### Function spaces
        ##############################################################

        self.E_cg = FiniteElement('CG', self.mesh.ufl_cell(), 1)
        self.E_dg = FiniteElement('DG', self.mesh.ufl_cell(), 0)
        self.V_cg = FunctionSpace(self.mesh, self.E_cg)
        self.V_dg = FunctionSpace(self.mesh, self.E_dg)


        ### CG input functions
        ##############################################################

        self.input_functions = {}
        self.state_functions = {}

               
        ### Create submodels 
        ##############################################################

        # Momentum model   
        self.momentum_model = MomentumModel(self)
        # Mass model
        self.mass_model = MassModel(self)
        # Hydrology model     
        self.hydro_model = HydroModel(self)

        self.momentum_model.update()
        self.mass_model.update()
        self.hydro_model.update()
        
        
    # Load fields
    def load_fields(self, fields):
        for field_name in fields:
            self.input_functions[field_name] = Function(self.V_cg)
            input_array = self.fields[field_name][::-1]
            self.input_functions[field_name].vector().set_local(input_array)

            
    # Step the model forward by one time step
    def solve(self, params = {}):

        # Solve for velocity
        ##############################################################
        
        self.momentum_model.solve()

        
        # Solve continuity equation and levelset equation
        ##############################################################
        
        self.mass_model.solve()
        self.momentum_model.H.assign(project(self.mass_model.H, self.V_cg))

        min_thickness = self.model_params['mass_params']['min_thickness']
        
        H_vec = self.momentum_model.H.vector().get_local()  
        indexes = H_vec <= min_thickness
        H_vec[indexes] = min_thickness
     
        phi = project(self.mass_model.phi, self.V_cg).vector().get_local()
        indexes = phi <= min_thickness
        H_vec[indexes] = min_thickness

        self.momentum_model.H.vector().set_local(H_vec)
        self.mass_model.H0.assign(project(self.momentum_model.H, self.V_dg))
        

        # Solve hydrology model 
        ##############################################################

        self.hydro_model.solve()


    def get_state(self):
        pass
