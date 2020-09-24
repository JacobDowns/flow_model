from dolfin import dx, ds, dS, avg, jump, sqrt, SpatialCoordinate, Constant, \
    lhs, rhs
from model.support.expressions import Abs
import ufl as ufl
import numpy as np

class MassForm(object):
    """
    Set up the variational form for the mass balance equation.
    """

    def __init__(self, model):

        # DG thickness trial function
        H = model.dH
        # Levelset function
        phi = model.dphi
        # Rate of change of H
        dHdt = model.dHdt
        # Rate of change of levelset function
        dphidt = model.dphidt
        # Velocity
        u = model.u
        # Calving rate
        c = model.c
        # DG test function
        xsi = model.xsi
        # SMB
        adot = model.adot
        # Ice stream width
        width = model.model_wrapper.input_functions['width']
        # Inter element flux (upwind)
        uH = avg(u)*avg(H*width) + 0.5*Abs(avg(width*u))*jump(H)
        uphi = avg(u)*avg(phi*width) + 0.5*Abs(avg(width*u))*jump(phi)
        
        ### Mass balance residual
        ###############################################################
        
        R_mass = (width*dHdt*xsi - adot*width*xsi)*dx
        R_mass += uH*jump(xsi)*dS

        self.R_mass = R_mass
        self.a_mass = lhs(R_mass)
        self.L_mass = rhs(R_mass)


        ### Levelset residual
        ###############################################################
        
        R_level = (width*dphidt*xsi - (adot - c)*width*xsi)*dx
        R_level += uphi*jump(xsi)*dS

        self.R_level = R_level
        self.a_level = lhs(R_level)
        self.L_level = rhs(R_level)
        
