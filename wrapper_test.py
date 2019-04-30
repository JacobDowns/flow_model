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


step_params = {}

step_params['pdd_params'] = {}
step_params['pdd_params']['monthly_dts'] = -30.*np.ones(12)
#step_params['pdd_params']['monthly_dps'] = np.ones(12)

wrapper.step(step_params)
dolfin.plot(wrapper.model.adot)
plt.show()





