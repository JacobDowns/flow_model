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
        # Trial function 
        dphi = TrialFunction(model_wrapper.V_cg)
        # Levelset at previous time step 
        phi0 = Function(model_wrapper.V_cg)
        # Time step 
        dt = Constant(model_wrapper.model_params['dt'])
        # Time derivative
        dphi_dt = (dphi - phi0) / dt

        self.phi = phi
        self.dphi = dphi
        self.phi0 = phi0
        self.dt = dt
        self.dphi_dt = dphi_dt
        self.dt = dt
        # Velocity
        self.u = model_wrapper.momentum_model.ubar
        # SMB
        self.adot = model_wrapper.mass_model.adot
        
        H = model_wrapper.momentum_model.H        
        min_thickness = model_wrapper.model_params['mass_params']['min_thickness']
        H_vec = H.vector().get_local()
        indexes = np.where(H_vec <= min_thickness)[0]
        index = indexes.max()
        self.x = model_wrapper.x
    
        L = model_wrapper.x[::-1][index]
        self.set_levelset(L)
        phi0.assign(phi)


        ### Variational forms
        ########################################################################

        # Test function
        v = TestFunction(model_wrapper.V_cg)
        self.v = v
        levelset_form = LevelsetForm(self)
        self.a_levelset = levelset_form.a_levelset
        self.L_levelset = levelset_form.L_levelset

        
        
    def set_levelset(self, L):
        x = self.model_wrapper.x[::-1]
        phi_vec = np.zeros_like(self.phi.vector().get_local())
        phi_vec = (x - L) / 1e3
        self.phi.vector().set_local(phi_vec)
        
        
    # Step the model forward by one time step
    def solve(self):
        solve(self.a_levelset == self.L_levelset, self.phi, [])
        self.phi0.assign(self.phi)
        phi_vec = self.phi0.vector()[::-1]
        i = np.where(phi_vec <= 0.)[0].max()
        phi_i0 = phi_vec[i]
        phi_i1 = phi_vec[i+1]
        w = abs(phi_i0) / (phi_i1 - phi_i0)
        x_i0 = self.x[i]
        x_i1 = self.x[i+1]
        L = w*x_i0 + (1.-w)*x_i1
        self.set_levelset(L)
        print("L", L)
        #df.plot(self.phi)
        #plt.show()
