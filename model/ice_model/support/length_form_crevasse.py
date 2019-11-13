from dolfin import ds, Constant
from ufl import max_value
from model.support.expressions import Abs, logistic

class LengthForm(object):
    """
    Set up the variational form for length using a "zero stress" 
    calving model.
    """
    
    def __init__(self, model):

        # Depth averaged velocity
        ubar = model.ubar
        # CG thickness
        H_c = model.H_c
        # Glacier length
        L = model.L
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
        # Domain length
        domain_len = model.domain_len
        # Real test function
        chi = model.chi
        # Boundary measure
        ds1 = ds(subdomain_data = model.boundaries)
        # Boundary calving (to prevent ice from flowing out of the domain)
        self.extra_calving = model.model_wrapper.input_functions['extra_calving']
        # Crevasse depth
        grounded = logistic(Bhat - B, k = .25, y0 = 25.)
        self.grounded = grounded 
        # Stretching rate
        R_xx = Constant(1e-16**(-1./3.))*abs(ubar.dx(0) / L + Constant(1e-16))**(1./3.)
        crevasse_depth = (-R_xx + rho*g*H_c - rho_w*g*D) / ((rho - rho_w)*g)*grounded
        #self.crevasse_depth = crevasse_depth
        self.crevasse_depth = max_value(Constant(min_thickness), crevasse_depth) + self.extra_calving
        # Length form
        R_length = (H_c - self.crevasse_depth)*chi*ds1(1)
        self.R_length = R_length
