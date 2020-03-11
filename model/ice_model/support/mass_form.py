from dolfin import dx, ds, dS, avg, jump, sqrt, SpatialCoordinate, Constant
from model.support.expressions import Abs
import ufl as ufl

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
        ds1 = ds(subdomain_data = model.boundaries)
        # SMB
        adot = model.adot
        # Ice stream width
        width = model.width / Constant(1000.)
        # Spatial coordinate
        x_spatial = SpatialCoordinate(model.mesh)
        # Grid velocity
        v = dLdt*x_spatial[0]
        # Flux velocity
        q_vel = ubar - v
        # Flux
        q_flux = q_vel*H*width
        # Inter element flux (upwind)
        uH = avg(q_vel)*avg(H*width) + 0.5*Abs(avg(width*q_vel))*jump(H)
        # Time partial of width
        dwdt = dLdt*x_spatial[0]/L*width.dx(0)

        
        ### Mass balance residual
        ########################################################################
        R_mass = (L*width*dHdt*xsi + L*H*dwdt*xsi + H*width*dLdt*xsi - L*width*adot*xsi)*dx
        R_mass += uH*jump(xsi)*dS
        R_mass += (q_vel / sqrt(q_vel**2 + Constant(1e-10))) * q_flux*xsi*ds1(1)
        
        self.R_mass = R_mass
