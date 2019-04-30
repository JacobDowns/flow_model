import numpy as np
from dolfin import *
from hydro_params import hydro_params

"""
A subglacial hydrology model that computes effective pressure. 
This is the dumbest hydrology model ever. 
"""

class HydroModel(object):

    def __init__(self, model_wrapper, params= {}):

        self.model_wrapper = model_wrapper
        self.hydro_params = hydro_params
        self.hydro_params.update(params)

        # Overburden pressure
        P_0 = Constant(self.model_wrapper.model.ice_constants['rho']*self.model_wrapper.model.ice_constants['g'])*self.model_wrapper.model.H_c
        # Overburden fraction
        P_frac = Constant(self.hydro_params['P_frac'])
        # Water pressure
        P_w = P_frac*P_0
        # Effective pressure
        self.model_wrapper.model.N = P_0 - P_w


    def update(self, params = {}):
        pass


