import numpy as np
from dolfin import Function

"""
A subglacial hydrology model that computes effective pressure. 
This is the dumbest hydrology model ever. 
"""

class HydroModel(object):

    def __init__(self, model_wrapper):
        self.model_wrapper = model_wrapper
        self.hydro_params = model_wrapper.model_params['hydro_params']
        self.N = Function(model_wrapper.V_cg)

    def solve(self):
        # Ice density
        rho = self. model_wrapper.model_params['physical_constants']['rho']
        # Gravitational constant
        g = self.model_wrapper.model_params['physical_constants']['g']
        # Overburden fraction
        P_frac = self.hydro_params['P_frac']
        # Ice thickness
        H = self.model_wrapper.momentum_model.H.vector().get_local()
        # Overburden pressure
        P_0 = rho*g*H
        # Set the effective pressure
        self.N.vector()[:] = (1. - P_frac)*P_0
