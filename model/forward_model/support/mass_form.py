import numpy as np
from dolfin import *

class MassForm(object):
    """
    Set up the variational form for the mass balance equation.
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
        # Velocity
        ubar = model.ubar
        # DG test function
        xsi = model.xsi
        # Boundary measure
        ds1 = dolfin.ds(subdomain_data = model.boundaries)
        # SMB expression
        adot_prime = model.adot_prime
        # Ice stream width
        width = model.width
        # Spatial coordinate
        x_spatial = SpatialCoordinate(model.mesh)
        # Grid velocity
        v = dLdt*x_spatial[0]
        # Flux velocity
        q_vel = ubar - v
        # Flux
        q_flux = q_vel*H*width
        # Inter element flux (upwind)
        uH = avg(q_vel)*avg(H*width) + 0.5*abs(avg(width*q_vel))*jump(H)
        # Time partial of width
        dwdt = dLdt*x_spatial[0]/L*width.dx(0)


        ### Mass balance residual
        ########################################################################
        R_mass = (L*width*dHdt*xsi + L*H*dwdt*xsi + H*width*dLdt*xsi - L*width*adot_prime*xsi)*dx
        R_mass += uH*jump(xsi)*dS
        R_mass += (q_vel / sqrt(q_vel**2 + Constant(1e-10))) * q_flux*xsi*ds1(1)

        self.R_mass = R_mass
