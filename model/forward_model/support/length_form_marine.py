import numpy as np
from dolfin import *

class LengthForm(object):
    """
    Set up the variational form for length for a marine terminating glacier. 
    """
    
    def __init__(self, model):

        # DG thickness
        H = model.H
        # CG thickness
        H_c = model.H_c
        # Bed elevation
        B = model.B
        # Sea level
        sea_level = model.sea_level
        # Density of ice
        rho = model.constants['rho']
        # Density of water
        rho_w = model.constants['rho_w']
        # Min. thickness
        min_thickness = model.min_thickness
        # Real test function
        chi = model.chi
        # Boundary measure
        ds1 = dolfin.ds(subdomain_data = model.boundaries)

        def softplus(y1,y2,alpha=1):
            # The softplus function is a differentiable approximation
            # to the ramp function.  Its derivative is the logistic function.
            # Larger alpha makes a sharper transition.
            return dolfin.Max(y1,y2) + (1./alpha)*dolfin.ln(1.+dolfin.exp(alpha*(dolfin.Min(y1,y2)-dolfin.Max(y1,y2))))

        # Grounding line thickness
        H_g = softplus(Constant(rho_w / rho)*(sea_level - B), min_thickness, alpha = 0.9)
        R_length = (H_c - H_g)*chi*ds1(1)
        self.R_length = R_length
