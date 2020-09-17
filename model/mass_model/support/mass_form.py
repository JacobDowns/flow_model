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
        # Velocity
        u = model.model_wrapper.momentum_model.ubar
        # DG test function
        xsi = model.xsi
        # SMB
        adot = model.adot
        # Ice stream width
        width = model.model_wrapper.input_functions['width']
        # Flux
        q = u*H*width
        # Inter element flux (upwind)
        uH = avg(u)*avg(H*width) + 0.5*Abs(avg(width*u))*jump(H)

        
        ### Mass balance residual
        ########################################################################
        R_mass = (width*dHdt*xsi - adot*width*xsi)*dx
        R_mass += uH*jump(xsi)*dS
        self.R_mass = R_mass
