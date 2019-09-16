import numpy as np
from dolfin import *
from model.support.expressions import softplus

class LengthForm(object):
    """
    Dirichlet condition for a land terminating glacier.
    """
    
    def __init__(self, model):

        # CG thickness
        H_c = model.H_c
        # Real test function
        chi = model.chi
        # Boundary measure
        ds1 = model.ds1
        # Length residual
        R_length = (H_c - Constant(15.))*chi*ds1(1)
        self.R_length = R_length
