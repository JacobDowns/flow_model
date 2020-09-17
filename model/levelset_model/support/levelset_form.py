from dolfin import dx, Constant, sqrt
from model.support.expressions import Abs
import ufl as ufl
import numpy as np

class LevelsetForm(object):
    """
    Set up the variational form for the mass balance equation.
    """

    def __init__(self, model):

        # Rate of change of phi
        dphi_dt = model.dphi_dt
        # Velocity
        u = model.u
        # SMB
        adot = model.adot
        # Ice hardness
        b = model.model_wrapper.model_params['momentum_params']['b']
        # Calving rate
        sigma = Constant(np.sqrt(3.)*b)*(0.5*u.dx(0)**2 + Constant(1e-16))**(1./3.)
        sigma_max = model.model_wrapper.model_params['levelset_params']['sigma_max']
        # Test funciton
        v = model.v
        
        
        ### Mass balance residual
        ########################################################################
        R_levelset = dphi_dt*v*dx + u*v*dx - adot*v*dx
        self.R_levelset = R_levelset
