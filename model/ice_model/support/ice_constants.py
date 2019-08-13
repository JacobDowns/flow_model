from dolfin import Constant

ice_constants = {}
# Seconds per year
ice_constants['spy'] = 60**2*24*365
# Ice density (kg / m^3)
ice_constants['rho'] = 917.
# Seawater density (kg / m^3)
ice_constants['rho_w'] = 1029.0
# Gravitational acceleration (m / s^2)
ice_constants['g'] = 9.81
# Glen's exponent
ice_constants['n'] = 3.0
# Sliding law exponent
ice_constants['m'] = 3.0
# Ice hardness
ice_constants['b'] = (3.5e-25 * ice_constants['spy'])**(-1./ice_constants['n'])
# Basal sliding law constants:
ice_constants['mu'] = Constant(1.0)
ice_constants['A_s'] = Constant(ice_constants['rho']*ice_constants['g']*315.0/500.)
# Minimum thickness (m)
ice_constants['min_thickness'] = 15.0
# Sea level (m)
ice_constants['sea_level'] = 0.0 #-10000.
# Default basal traction (Pa a / m)
ice_constants['beta2'] = 1.6e-3
# Fraction above flotation where calving begins
ice_constants['q'] = 0.15
