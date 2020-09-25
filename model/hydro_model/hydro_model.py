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
        self.constants = model_wrapper.model_params['physical_constants']
        
        # Effective pressure
        self.N = model_wrapper.cg_funcs['N']
        # Overburden pressure
        self.P_0 = model_wrapper.cg_funcs['P_0']
        # Water pressure
        self.P_w = model_wrapper.cg_funcs['P_w']
        # Ice thickness
        self.H = model_wrapper.cg_funcs['H']

        
    def solve(self):
        # Ice density
        rho = self.constants['rho']
        # Gravitational constant
        g = self.constants['g']
        # Overburden fraction
        P_frac = self.hydro_params['P_frac']
        # Ice thickness
        H = self.H.vector().get_local()
        # Overburden pressure
        P_0 = rho*g*H
        self.P_0.vector().set_local(P_0)
        # Water pressure
        P_w = P_frac*P_0
        self.P_w.vector().set_local(P_w)
        # Set the effective pressure
        self.N.vector().set_local(P_0 - P_w)
