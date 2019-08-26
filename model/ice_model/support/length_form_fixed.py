import numpy as np
from dolfin import *

class LengthForm(object):
    """
    Set up the variational form for fixed length. 
    """
    
    def __init__(self, model):

      
        # Glacier length
        L = model.L
        # Last glacier length
        L0 = model.L0
        # Real test function
        chi = model.chi
        # Length form
        R_length = (L - L0)*chi*dx
        self.R_length = R_length
