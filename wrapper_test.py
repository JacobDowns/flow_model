import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
from model.pdd_wrapper import PDDWrapper

inputs = {}

# Model inputs
inputs['model_inputs'] = {}
model_inputs = inputs['model_inputs']

model_inputs['input_file_name'] = '/home/jake/flow_model/output_files/inversions/cov1/flowline.h5'
model_inputs['state_file_name'] = '/home/jake/flow_model/output_files/inversions/cov1/state0.h5'
model_inputs['dt'] = 1./3.

wrapper = PDDWrapper(inputs)

sea_level = np.linspace(-75., 25, 750)
for i in range(750):
    params = {}
    params['pdd_params'] = {'monthly_dts' : -5.*np.ones(12)}
    params['ice_params'] = {'sea_level' : sea_level[i]}
    #params['step_params'] = {}
    print(sea_level[i])
    wrapper.step(step_params)
    print(float(wrapper.model.sea_level))
    print(wrapper.model.H0_c.vector().get_local().min())
        
dolfin.plot(wrapper.model.S0_c)
dolfin.plot(wrapper.model.B)
plt.show()






