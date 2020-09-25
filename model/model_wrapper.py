import numpy as np
from dolfin import FiniteElement, FunctionSpace, set_log_level, \
    MeshFunction, near, Function, IntervalMesh, project, FacetNormal, \
    Measure
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
        
        ### Model mesh
        ##############################################################

        self.x = inputs['x']
        self.mesh = IntervalMesh(len(self.x) - 1, self.x.min(), self.x.max())
        self.mesh.coordinates()[:] = self.x[:, np.newaxis]

        self.terminus_marker = MeshFunction('size_t', self.mesh,
                                             self.mesh.topology().dim()-1, 0)
        self.ds_t = Measure('dS', domain=self.mesh,
                            subdomain_data=self.terminus_marker)

            
        #### Function spaces
        ##############################################################

        self.E_cg = FiniteElement('CG', self.mesh.ufl_cell(), 1)
        self.E_dg = FiniteElement('DG', self.mesh.ufl_cell(), 0)
        self.V_cg = FunctionSpace(self.mesh, self.E_cg)
        self.V_dg = FunctionSpace(self.mesh, self.E_dg)


        ### CG input functions
        ##############################################################

        self.cg_funcs = {}
        self.dg_funcs = {}

        for field in model_params['cg_fields']:
            self.cg_funcs[field] = Function(self.V_cg)

        for field in model_params['dg_fields']:
            self.dg_funcs[field] = Function(self.V_dg)
            
        for field in inputs['cg_inputs']:
            input_vec = inputs['cg_inputs'][field][::-1]
            self.cg_funcs[field].vector().set_local(input_vec)

        if 'dg_inputs' in inputs:
            for field in inputs['dg_inputs']:
                input_vec = inputs['dg_inputs'][field][::-1]
                self.dg_funcs[field].vector().set_local(input_vec)

            
        ### Create submodels 
        ##############################################################

        self.update()
        
        # Momentum model   
        self.momentum_model = MomentumModel(self)
        # Mass model
        self.mass_model = MassModel(self)
        # Hydrology model     
        self.hydro_model = HydroModel(self)


    def update(self):

        rho = self.model_params['physical_constants']['rho']
        rho_w = self.model_params['physical_constants']['rho_w']
        g = self.model_params['physical_constants']['g']

        ### Update ice base
        B = self.cg_funcs['B'].vector().get_local()
        H = self.cg_funcs['H'].vector().get_local()
        Bhat = np.maximum(B,-(rho/rho_w)*H)
        self.cg_funcs['Bhat'].vector().set_local(Bhat)

        ### Update basal traction on floating parts of domain
        depth = Bhat - B
        self.cg_funcs['depth'].vector().set_local(depth)
        floating = np.ones_like(depth)
        floating[depth > 1.] = 1e-12
        self.cg_funcs['floating'].vector().set_local(floating)

        ### Terminus marker
        min_thickness = self.model_params['mass_params']['min_thickness']
        index = np.where(H == min_thickness)[0].max() + 1
        terminus = np.zeros_like(H)
        terminus[index] = 1.
        self.terminus_marker.array()[:] = terminus[::-1]

        ### Update DG thickness
        self.dg_funcs['H'].assign(project(self.cg_funcs['H'], self.V_dg))
        self.dg_funcs['H0'].assign(self.dg_funcs['H'])
        
        
    # Step the model forward by one time step
    def solve(self, params = {}):

        self.update()
        
        # Solve hydrology model 
        ##############################################################

        self.hydro_model.solve()

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


    def get_state(self):
        pass


    def get_field(self, field):
        return self.cg_funcs[field].vector().get_local()[::-1]
