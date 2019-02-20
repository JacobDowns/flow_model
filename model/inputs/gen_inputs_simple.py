from dolfin import *
import numpy as np
from common_inputs import *
import matplotlib.pyplot as plt

"""
Generate a full input file from is_flowline_partial.h5. This script just creates
an initial idealized ice sheet profile with the given length and thickness at the
divide.
"""



# Initial glacier length
L_init = 300e3
# Maximum ice thickness
H_max = 3000.


### Load bed inputs
################################################################################
input_options = {}
input_options['cg_fields'] = ['B']
input_options['interp_fields'] = ['B']
input_options['r_fields'] = []
input_options['dg_fields'] = []
inputs = CommonInputs('is_flowline_300.h5', input_options)


### Initialize ice thickness
################################################################################
# Continuous thickness
H0 = Function(inputs.V_cg)
# DG thickness
H0_c = Function(inputs.V_dg)
# Get the bed elevation at the desired terminus position
B_term = inputs.get_interp_value('B', L_init)
inputs.update_interp_all(L_init)


# Surface expression
class SExp(Expression):
    def eval(self,values,x):
        values[0] = np.sqrt((H_max + B_term)**2*(1. - x[0])) + B_term

# Surface expression
class SExp(Expression):
    def eval(self,values,x):
        values[0] = np.sqrt((H_max + B_term)**2*x[0]) + B_term

S = project(SExp(degree = 1), inputs.V_cg)
# Compute initial thickness
H0_c.assign(project(S - inputs.input_functions['B'], inputs.V_cg))
# As DG function
H0.assign(project(H0_c, inputs.V_dg))
# Initial glacier length
L0 = Function(inputs.V_r)
L0.assign(Constant(L_init))


### Width and basal tractions
################################################################################
# Ice stream width
width = Function(inputs.V_cg)
width.interpolate(Constant(1000.))
# Basal traction
beta2 = Function(inputs.V_cg)
beta2.interpolate(Constant(1e-3))


### Update hdf5 file
################################################################################
inputs.input_file.write(H0_c, 'H0_c')
inputs.input_file.write(H0, 'H0')
inputs.input_file.write(L0, 'L0')
inputs.input_file.write(width, 'width')
inputs.input_file.write(beta2, 'beta2')
inputs.input_file.close()
