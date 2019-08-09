import numpy as np
from dolfin import *
from ...support.expressions import softplus

class LengthForm(object):
    """
    Set up the variational form for length for a marine terminating 
    glacier with a simple calving law. 
    """
    
    def __init__(self, model):

        # CG thickness
        H_c = model.H_c
        # Water depth
        D = model.D
        # Density of ice
        rho = model.ice_constants['rho']
        # Density of water
        rho_w = model.ice_constants['rho_w']
        # Min. thickness
        min_thickness = model.ice_constants['min_thickness']
        # Calving parameter
        q = model.ice_constants['q']
        # Calving thickness
        #H_calving = softplus(Constant((rho_w/rho)*(1. + q))*D,
        #                     min_thickness, alpha = 0.1)
        H_calving = Constant(30.0)
        self.H_calving = H_calving
        # Real test function
        chi = model.chi
        # Boundary measure
        ds1 = model.ds1
        # Length residual
        R_length = (H_c - H_calving)*chi*ds1(1)
        self.R_length = R_length
