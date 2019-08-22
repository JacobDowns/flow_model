import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
from model.pdd_wrapper import PDDWrapper
from scipy.special import expit

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
model_inputs['L0'] = 232e3
# Along flow coordinate
x = np.linspace(0., model_inputs['domain_len'], 450)

### Bedrock elevation
######################################################################
model_inputs['B'] = 1000.*(1.0 - expit((x - 250e3) / 7.5e3 )) - 980. #500.0 - x/2.5e2


### Surface
######################################################################
model_inputs['S_ref'] = np.zeros_like(x)
indexes = x <= model_inputs['L0']
model_inputs['S_ref'][indexes] = np.sqrt((2500.)**2*(1. - (x[indexes] / model_inputs['L0']))) - 300.

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
#T = -10.*np.ones(12)
# Some monthly precip. averages (m.w.e. / a)
P = .65*np.ones(12)
chi = x / x.max()
T = -10. + chi*12.

for i in range(12):
    model_inputs['T' + str(i)] = T
    model_inputs['P' + str(i)] = P[i]*np.ones_like(x)


### Run the model
######################################################################
wrapper = PDDWrapper(inputs)

#dolfin.plot(a)
#plt.show()

#from model.support.expressions import *
#dolfin.plot(Constant(1.0) - logistic(a-b, y0=Constant(1.)))
#plt.show()

# Bed
B = wrapper.original_cg_functions['B'].vector().get_local()[::-1]
# Domain length
domain_len = wrapper.domain_len
# Mesh coordinates
x = wrapper.mesh_coords
# Run the model 
for i in range(3000):
    print(i)
    wrapper.step()

### Ice thickness
######################################################################
H_c = wrapper.model.H0_c.vector().get_local()[::-1]
# Water depth
D = project(wrapper.model.D).vector().get_local()[::-1]
# Surface
S = wrapper.model.S0_c.vector().get_local()[::-1]
# Glacier length
L = float(wrapper.model.L0)
l = project(wrapper.model.l).vector().get_local()[::-1]
# Stretching
crevasse_depth = project(wrapper.model.length_form.crevasse_depth).vector().get_local()[::-1]
# Ice base
Bhat = project(wrapper.model.Bhat).vector().get_local()[::-1]
# Velocity
u0 = project(wrapper.model.ubar).vector().get_local()[::-1]
u1 = project(wrapper.model.udef).vector().get_local()[::-1]
xx, yy = np.meshgrid(np.linspace(0,1, len(u0)), np.linspace(0,1,35))
u0, ys = np.meshgrid(u0, np.linspace(0,1,35))
u1, ys = np.meshgrid(u1, np.linspace(0,1,35))
u = u0 + (u1 / 4.) * (5*ys**4 - 1.)
Bhat = project(wrapper.model.Bhat).vector().get_local()[::-1]
yy = H_c * yy + Bhat




# Piezometric head
rho_i = 917.
rho_w = 1000.
g = 9.81

"""
# Zero-stress calving
def z(H, D, e_xx):
    #B = 2.207e-3*np.exp( (3155./(T + 273.15)) - (0.16612 / (0.24 - T)**1.17))
    Rxx = 2.*(1e-16)**(-1./3.)*np.abs(e_xx + 1e-16)**(1./3.)
    z = (-Rxx + rho_i*g*H - rho_w*g*D) / ((rho_i - rho_w)*g)
    return z

#plt.plot(e_xx)
#plt.show()
"""
plt.plot(x * domain_len, B, 'k', lw = 2)
plt.plot(x * L, Bhat, 'k', lw = 2)
plt.plot(x * L, S, 'k', lw = 2)
plt.plot(x * domain_len, 0.*x, 'k--', lw = 2)
plt.contourf(xx*L, yy[::-1], u, levels = np.linspace(0., u.max(), 250), cmap = 'RdBu_r')
plt.colorbar()

plt.show()



