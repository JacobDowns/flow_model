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
        # Fields that need to be loaded
        self.fields = ['adot']
        # Load model fields
        model_wrapper.load_fields(self.fields)

        
        ### Initialize functions
        ########################################################################

        # DG thickness
        H = Function(model_wrapper.V_dg)
        # Levelset function
        phi = Function(model_wrapper.V_dg)
        # Thickness at last time step
        H0 = Function(model_wrapper.V_dg)
        H.assign(project(model_wrapper.momentum_model.H, model_wrapper.V_dg))
        H0.assign(H)
        # Trial function 
        dH = TrialFunction(model_wrapper.V_dg)
        # Levelset trial function
        dphi = TrialFunction(model_wrapper.V_dg)
        # Width
        width = model_wrapper.input_functions['width']
        # Velocity
        u = model_wrapper.momentum_model.ubar
        # Time step
        dt = Constant(1.)
        # Time derivative
        dHdt = (dH - H0) / dt
        # Levelset time derivative
        dphidt = (dphi - H0) / dt
        # SMB
        adot = Function(model_wrapper.V_cg)
        adot.assign(model_wrapper.input_functions['adot'])
        # Velocity
        u = model_wrapper.momentum_model.ubar
        # Calving
        c = Function(model_wrapper.V_cg)
        
        self.u = u
        self.adot = adot 
        self.H = H
        self.phi = phi
        self.H0 = H0
        self.dt = dt
        self.dHdt = dHdt
        self.dphi = dphi
        self.dphidt = dphidt
        self.c = c
        
        
        ### Variational forms
        ########################################################################

        # Test function
        xsi = TestFunction(model_wrapper.V_dg)

        self.xsi = xsi
        self.dH = dH
        
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
        self.c.vector().set_local(c*0. + 0.*75.)
        
        
    # Step the model forward by one time step
    def solve(self):
        
        self.update()

        solve(self.a_mass == self.L_mass, self.H, [])
        solve(self.a_level == self.L_level, self.phi, [])
        
        self.H0.assign(self.H)
