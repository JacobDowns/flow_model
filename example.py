import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
from model.pdd_wrapper import PDDWrapper

inputs = {}

### Model inputs
#############################################################
inputs['model_inputs'] = {}
model_inputs = inputs['model_inputs']

# Time step
model_inputs['dt'] = 1./3.
# Length of full domain
model_inputs['domain_len'] = 450e3
# Initial glacier length 
model_inputs['L0'] = 250e3
# Along flow coordinate
x = np.linspace(0., model_inputs['domain_len'], 450)

### Bedrock elevation
######################################################################
model_inputs['B'] = 500.0 - x/5e2


### Surface
######################################################################
model_inputs['S_ref'] = np.zeros_like(x)
indexes = x <= model_inputs['L0']
model_inputs['S_ref'][indexes] = np.sqrt((3000.)**2*(1. - (x[indexes] / model_inputs['L0']))) - 300.


# Correct the surface so it's not under the bed...
index = np.where(model_inputs['S_ref'] - model_inputs['B'] <= 15.)[0].min()
model_inputs['S_ref'][index:] = model_inputs['B'][index:]


plt.plot(model_inputs['B'])
plt.plot(model_inputs['S_ref'])
plt.plot(np.zeros_like(model_inputs['B']))
plt.show()


### Ice thickness
######################################################################

model_inputs['H0_c'] = model_inputs['S_ref'] - model_inputs['B']

### Ice stream width
######################################################################

model_inputs['width'] = 1000.*np.ones_like(x)

### Temperature and precip
######################################################################
# Some monthly temp. averages (C)
T = -20.*np.ones(12)
# Some monthly precip. averages (m.w.e. / a)
P = .75*np.ones(12)

for i in range(12):
    model_inputs['T' + str(i)] = T[i]*np.ones_like(x)
    model_inputs['P' + str(i)] = P[i]*np.ones_like(x)


### Run the model
######################################################################
wrapper = PDDWrapper(inputs)

a = project(Expression("10.0*(x[0]-0.5)", degree=1), wrapper.V_cg)
b = Function(wrapper.V_cg)

#dolfin.plot(a)
#plt.show()

from model.support.expressions import *
dolfin.plot(Constant(1.0) - logistic(a-b, y0=Constant(1.)))
plt.show()


quit()

# Bed
B = wrapper.original_cg_functions['B'].vector().get_local()[::-1]
# Domain length
domain_len = wrapper.domain_len
# Mesh coordinates
x = wrapper.mesh_coords
# Run the model 
for i in range(2000):
    wrapper.step()

#dolfin.plot(wrapper.model.adot)
#plt.show()
#quit()
# Surface
S = wrapper.model.S0_c.vector().get_local()[::-1]
# Glacier length
L = float(wrapper.model.L0)
l = project(wrapper.model.l).vector().get_local()[::-1]
Bhat = project(wrapper.model.Bhat).vector().get_local()[::-1]



quit()
#plt.plot(l)


plt.plot(x * domain_len, B, 'k')
plt.plot(x * L, S, 'r')
plt.plot(x*domain_len, 0.*x, 'k--')
plt.plot(x * L, Bhat, 'r--')
plt.show()




