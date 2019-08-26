import numpy as np
from dolfin import *
from ...support.expressions import *

class LengthForm(object):
    """
    Set up the variational form for fixed length. 
    """
    
    def __init__(self, model):

        # Depth averaged velocity
        ubar = model.ubar
        # CG thickness
        H_c = model.H_c
        # Glacier length
        L = model.L
        # Last glacier length
        L0 = model.L0
        # Bed elevation
        B = model.B
        # Ice base
        Bhat = model.Bhat
        # Ocean depth
        D = model.D
        # Glen's n
        n = model.ice_constants['n']
        # Density of ice
        rho = model.ice_constants['rho']
        # Density of water
        rho_w = model.ice_constants['rho_w']
        # Rate factor
        A = model.ice_constants['A']
        # Gravitational constant
        g = model.ice_constants['g']
        # Min. thickness
        min_thickness = model.ice_constants['min_thickness']
        # Real test function
        chi = model.chi
        # Boundary measure
        ds1 = dolfin.ds(subdomain_data = model.boundaries)
        # Crevasse deptsh
        grounded = logistic(Bhat - B, k = 0.05, y0 = 100.)
        # Stretching rate
        R_xx = Constant(2.*1e-16**(-1./3.))*abs(ubar.dx(0) / L + Constant(1e-16))**(1./3.)
        crevasse_depth = (-R_xx + rho*g*H_c - rho_w*g*D) / ((rho - rho_w)*g)
        self.crevasse_depth = softplus(Constant(min_thickness), crevasse_depth, alpha = .007)
        # Grounding line thickness
        H_g = self.crevasse_depth

        # Calving parameter
        q = model.ice_constants['q']
        # Calving thickness
        self.H_calving = softplus(Constant((rho_w/rho)*(1. + q))*D,
                             min_thickness, alpha = 0.1)
        
        #self.H_g = H_g
        # Length form
        R_length = (L - L0)*chi*dx
        #R_length = (H_c - self.H_calving)*chi*ds1(1)
        self.R_length = R_length
