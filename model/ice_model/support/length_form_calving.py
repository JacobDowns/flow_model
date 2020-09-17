import numpy as np
from dolfin import *
#from ...support.expressions import softplus
from model.support.expressions import Abs, logistic

class LengthForm(object):
    """
    Set up the variational form for length for a marine terminating 
    glacier with a simple calving law. 
    """
    
    def __init__(self, model):

        # CG thickness
        H_c = model.H_c
        # Ice base
        Bhat = model.Bhat
        # Bed
        B = model.B
        # Density of ice
        rho = model.ice_constants['rho']
        # Density of water
        rho_w = model.ice_constants['rho_w']
        # Min. thickness
        H_calving = Constant(20.)
        self.H_calving = H_calving
        # Real test function
        chi = model.chi
        # Boundary measure
        ds1 = model.ds1

        H_calving = logistic(Bhat - B, k = .01, y0 = 50.)**2 * Constant(250.)
                
        R_length = (H_c - H_calving)*chi*ds1(1)
        self.R_length = R_length
