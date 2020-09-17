from dolfin import Constant
from util.custom_dict import CustomDict

momentum_params = CustomDict()
# Glen's exponent
momentum_params['n'] = 3.0
# Rate factor
momentum_params['A'] = 3.5e-25
# Ice hardness
momentum_params['b'] = (momentum_params['A'] * 60**2*24*365)**(-1./momentum_params['n'])
# Basal sliding law constants:
momentum_params['mu'] = Constant(1.0)
momentum_params['A_s'] = Constant(ice_constants['rho']*ice_constants['g']*315.0/500.)
# Minimum thickness (m)
ice_constants['min_thickness'] = 15.0
# Sea level (m)
ice_constants['sea_level'] = 0.#-10000.
# Default basal traction (Pa a / m)
ice_constants['beta2'] = 1.2e-3
