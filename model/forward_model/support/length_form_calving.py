import numpy as np
from dolfin import *

class LengthForm(object):
    """
    Set up the variational form for length, or more specically, the H(x=1)=0
    boundary condition used to determine length.
    """

    def __init__(self, model):

        # DG thickness
        H = model.H
        # Rate of change of H
        dHdt = model.dHdt
        # Ice sheet length
        L = model.L
        # Rate of change of L
        dLdt = model.dLdt
        # Real test function
        chi = model.chi
        # Boundary measure
        ds1 = dolfin.ds(subdomain_data = model.boundaries)
        # SMB expression
        adot_prime = model.adot_prime
        # Ice stream width
        width = model.width
        # Spatial coordinate
        x_spatial = SpatialCoordinate(model.mesh)
        # Time partial of width
        dwdt = dLdt*x_spatial[0]/L*width.dx(0)
        # Calving velocity
        U_calving = Constant(10.)*(H  - Constant(20.0))**2

        ### Length residual
        ### ====================================================================

        #R_length = (dLdt*width*H + L*dwdt*H + L*width*(dHdt - adot_prime))*chi*dx + U_calving*width*H*chi*ds1(1)
        R_length = dHdt*chi*ds1(1)
        self.R_length = R_length
