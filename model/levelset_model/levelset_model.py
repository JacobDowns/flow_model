from dolfin import Function, TrialFunction, parameters, project, TestFunction, solve
import dolfin as df
from dolfin import Constant, lhs, rhs
import numpy as np
import matplotlib.pyplot as plt
from model.levelset_model.support.levelset_form import *

class LevelsetModel(object):

    def __init__(self, model_wrapper):

        # Model inputs object
        self.model_wrapper = model_wrapper

        
        ### Initialize functions
        ########################################################################

        
        # Levelset function
        phi = Function(model_wrapper.V_cg)
        # Levelset at previous time step 
        phi0 = Function(model_wrapper.V_cg)
        # Time step 
        dt = Constant(model_wrapper.model_params['dt'])
        # Time derivative
        dphi_dt = (phi - phi0) / dt

        self.phi = phi
        self.phi0 = phi0
        self.dt = dt
        self.dphi_dt = dphi_dt
        self.dt = dt
        # Velocity
        self.u = model_wrapper.momentum_model.ubar
        self.adot = model_wrapper.mass_model.adot
        
        H = model_wrapper.momentum_model.H        
        min_thickness = model_wrapper.model_params['mass_params']['min_thickness']
        H_vec = H.vector().get_local()
        indexes = np.where(H_vec <= min_thickness)[0]
        index = indexes.max()
        L = model_wrapper.x[::-1][index]
        self.set_levelset(L)
        phi0.assign(phi)


        ### Variational forms
        ########################################################################

        # Test function
        v = TestFunction(model_wrapper.V_cg)
        # Trial function 
        dphi = TrialFunction(model_wrapper.V_cg)
        self.v = v
        self.dphi = dphi

        self.levelset_form = LevelsetForm(self)

        
    def set_levelset(self, L):
        x = self.model_wrapper.x[::-1]
        phi_vec = np.zeros_like(self.phi.vector().get_local())
        phi_vec = x - L
        self.phi.vector().set_local(phi_vec)
        
        
    # Step the model forward by one time step
    def solve(self):
        pass
