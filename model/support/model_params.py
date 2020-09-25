from dolfin import Constant
from util.custom_dict import CustomDict

model_params = CustomDict()

model_params['dt'] = 1.


### Physical constants
######################################################################

physical_constants = CustomDict()

# Seconds per year
physical_constants['spy'] = 60**2*24*365
# Ice density (kg / m^3)
physical_constants['rho'] = 917.
# Seawater density (kg / m^3)
physical_constants['rho_w'] = 1029.0
# Gravitational acceleration (m / s^2)
physical_constants['g'] = 9.81

model_params['physical_constants'] = physical_constants

# Lists of CG and DG model fields
model_params['cg_fields'] = ['B', 'H', 'width', 'beta2', 'adot', 'depth',
                             'N', 'P_0', 'P_w', 'Bhat', 'floating',
                             'u', 'c']
model_params['dg_fields'] = ['H', 'phi', 'H0']


### Momentum model parameters
######################################################################

momentum_params = CustomDict()

# Glen's exponent
momentum_params['n'] = 3.0
# Rate factor
momentum_params['A'] = 3.5e-25
# Ice hardness
momentum_params['b'] = (momentum_params['A'] * 60**2*24*365)**(-1./momentum_params['n'])


model_params['momentum_params'] = momentum_params


### Mass model parameters
######################################################################

mass_params = CustomDict()

# Minimum thickness (m)
mass_params['min_thickness'] = 10.0
mass_params['sigma_max'] = 100e3

model_params['mass_params'] = mass_params


### Hydro model parameters
######################################################################

hydro_params = CustomDict()

# Overburden pressure fraction
hydro_params['P_frac'] = 0.85

model_params['hydro_params'] = hydro_params
