import numpy as np
from dolfin import FiniteElement, FunctionSpace, set_log_level, \
    MeshFunction, near, Function, IntervalMesh, project

from model.support.model_params import *
from model.momentum_model.momentum_model import MomentumModel
from model.mass_model.mass_model import MassModel
from model.levelset_model.levelset_model import LevelsetModel
from model.hydro_model.hydro_model import HydroModel

#set_log_level(40)

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
        # Levelset model
        self.levelset_model = LevelsetModel(self)
        

    # Load fields
    def load_fields(self, fields):
        for field_name in fields:
            self.input_functions[field_name] = Function(self.V_cg)
            input_array = self.fields[field_name][::-1]
            self.input_functions[field_name].vector().set_local(input_array)


    # Step the model forward by one time step
    def solve(self, params = {}):
        self.momentum_model.solve()
        self.mass_model.solve()
        self.momentum_model.H.assign(project(self.mass_model.H, self.V_cg))

        
