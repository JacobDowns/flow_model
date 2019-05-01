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
model_inputs['L0'] = 400e3
# Along flow coordinate
x = np.linspace(0., model_inputs['domain_len'], 450)

### Bedrock elevation
model_inputs['B'] = 200.*np.sin(x / 20e3)

### Surface
model_inputs['S_ref'] = np.zeros_like(x)
indexes = x <= model_inputs['L0']
model_inputs['S_ref'][indexes] = np.sqrt((3000.)**2*(1. - (x[indexes] / model_inputs['L0'])))
# Correct the surface so it's not under the bed...
index = np.where(model_inputs['S_ref'] - model_inputs['B'] <= 15.)[0].min()
model_inputs['S_ref'][index:] = model_inputs['B'][index:]

### Ice thickness
model_inputs['H0_c'] = model_inputs['S_ref'] - model_inputs['B']

### Ice stream width
model_inputs['width'] = 1000.*np.ones_like(x)

### Temperature and precip

# Some monthly temp. averages (C)
T = [-24.25, -24.87, -23.01, -15.43, -7.08, 0.08, 2.28, -0.29, -6.19, -13.70, -18.62, -22.41]
# Some monthly precip. averages (m.w.e. / a)
P = [0.300, 0.302, 0.320, 0.378, 0.367, 0.351, 0.270, 0.362, 0.430, 0.430, 0.489, 0.343]

for i in range(12):
    model_inputs['T' + str(i)] = T[i]*np.ones_like(x)
    model_inputs['P' + str(i)] = P[i]*np.ones_like(x)

wrapper = PDDWrapper(inputs)

for i in range(10):
    print(i)
    wrapper.step()

print(float(wrapper.model.L0))
dolfin.plot(wrapper.model.adot)
plt.show()




