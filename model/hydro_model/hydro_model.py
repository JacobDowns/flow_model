import numpy as np
from dolfin import *
from model.hydro_model.hydro_params import hydro_params
import matplotlib.pyplot as plt

"""
A subglacial hydrology model that computes effective pressure. 
This is the dumbest hydrology model ever. 
"""

class HydroModel(object):

    def __init__(self, model_wrapper, params= {}):

        self.model_wrapper = model_wrapper
        self.hydro_params = hydro_params
        self.hydro_params.update(params)


    def update(self, params = {}):
        # Ice density
        rho = self.model_wrapper.ice_model.ice_constants['rho']
        # Gravitational constant
        g = self.model_wrapper.ice_model.ice_constants['g']
        # Overburden fraction
        P_frac = self.hydro_params['P_frac']
        # Ice thickness
        H_vec = self.model_wrapper.ice_model.H0_c.vector().get_local()
        # Overburden pressure
        P_0 = rho*g*H_vec
        # Set the effective pressure
        self.model_wrapper.ice_model.N.vector()[:] = (1. - P_frac)*P_0
