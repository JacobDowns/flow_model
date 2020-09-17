from dolfin import Function, TrialFunction, parameters, project, TestFunction,\
    solve, lhs, rhs
import dolfin as df
import matplotlib.pyplot as plt
from model.mass_model.support.mass_form import *

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
        H0 = Function(model_wrapper.V_dg)
        H.assign(project(model_wrapper.momentum_model.H, model_wrapper.V_dg))
        H0.assign(H)
        # Width
        width = model_wrapper.input_functions['width']
        # Velocity
        u = model_wrapper.momentum_model.ubar
        # Time step
        dt = Constant(1.)
        # Time derivative
        dHdt = (H - H0) / dt
        # SMB
        adot = Function(model_wrapper.V_cg)
        adot.assign(model_wrapper.input_functions['adot'])
        
        self.adot = adot 
        self.H = H
        self.H0 = H0
        self.dt = dt
        self.dHdt = dHdt

        
        ### Variational forms
        ########################################################################

        # Test function
        xsi = TestFunction(model_wrapper.V_dg)
        # Trial function 
        dH = TrialFunction(model_wrapper.V_dg)
        self.xsi = xsi
        self.dH = dH
        
        # Residual 
        mass_form = MassForm(self)
        self.R = mass_form.R_mass

        self.a = lhs(self.R)
        self.L = rhs(self.R)

        
    # Step the model forward by one time step
    def solve(self, dt = 1.):
        self.dt.assign(dt)
        solve(self.R == 0, self.H, [], \
              solver_parameters={"newton_solver":{"absolute_tolerance":5e-6}})
        self.H0.assign(self.H)
