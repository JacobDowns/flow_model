from dolfin import Function, TrialFunction, parameters, project, TestFunction,\
    solve, lhs, rhs
import dolfin as df
import matplotlib.pyplot as plt
from model.mass_model.support.mass_form import *
from model.support.expressions import Abs
from ufl import max_value
from model.support.expressions import logistic

parameters['form_compiler']['cpp_optimize'] = True

class MassModel(object):

    def __init__(self, model_wrapper, params = {}):

        # Model inputs object
        self.model_wrapper = model_wrapper
        self.constants = model_wrapper.model_params['physical_constants']

        
        ### Initialize functions
        ########################################################################

        # DG thickness
        self.H = model_wrapper.dg_funcs['H']
        # Levelset function
        self.phi = model_wrapper.dg_funcs['phi']
        # Thickness at last time step
        self.H0 = model_wrapper.dg_funcs['H0']
        # Width
        self.width = model_wrapper.cg_funcs['width']
        # SMB
        self.adot =  model_wrapper.cg_funcs['adot']
        # Depth averaged velocity
        self.u = model_wrapper.cg_funcs['u']
        # Calving
        self.c = model_wrapper.cg_funcs['c']
        # Trial function 
        self.dH = TrialFunction(model_wrapper.V_dg)
        # Levelset trial function
        self.dphi = TrialFunction(model_wrapper.V_dg)
        # Time step
        self.dt = Constant(1.)
        # Time derivative
        self.dHdt = (self.dH - self.H0) / self.dt
        # Levelset time derivative
        self.dphidt = (self.dphi - self.H0) / self.dt
        # Test funciton
        self.xsi = TestFunction(model_wrapper.V_dg)

        
        ### Variational forms
        ########################################################################
        
        # Mass residual
        mass_form = MassForm(self)
        self.a_mass = mass_form.a_mass 
        self.L_mass = mass_form.L_mass
        
        # Levelset residual
        self.a_level = mass_form.a_level
        self.L_level = mass_form.L_level
        

    # Update time evolving fields
    def update(self):
        self.dt.assign(self.model_wrapper.model_params['dt'])
                       
        u = project(self.u, self.model_wrapper.V_cg).vector().get_local()
        b = self.model_wrapper.model_params['momentum_params']['b']
        dudx = project(self.u.dx(0), self.model_wrapper.V_cg).vector().get_local()
        dudx[dudx < 0.] = 0.
        sigma = np.sqrt(3.)*b*(0.5*dudx + 1e-16)**(2./3.)
        sigma_max = self.model_wrapper.model_params['mass_params']['sigma_max']
        c = np.absolute(u)*(sigma / sigma_max)
        floating = self.model_wrapper.momentum_model.floating.vector().get_local()
        c *= floating
        self.c.vector().set_local(c*0. + 45.)
        
        
    # Step the model forward by one time step
    def solve(self):
        
        self.update()

        solve(self.a_mass == self.L_mass, self.H, [])
        solve(self.a_level == self.L_level, self.phi, [])
        
        self.H0.assign(self.H)

    
