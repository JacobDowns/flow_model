import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
from model.pdd_wrapper import PDDWrapper

inputs = {}

# Model inputs
inputs['model_inputs'] = {}
model_inputs = inputs['model_inputs']

model_inputs['input_file_name'] = '/home/jake/flow_model/output_files/inversions/covariance_comparison/flowline.h5'
model_inputs['state_file_name'] = '/home/jake/flow_model/output_files/inversions/covariance_comparison/state0.h5'
model_inputs['dt'] = 1./3.

wrapper = PDDWrapper(inputs)

temps = []
for i in range(12):
    print(wrapper.input_functions['P' + str(i)].vector().get_local().mean())
    



quit()


step_params = {}

step_params['pdd_params'] = {}
#step_params['pdd_params']['monthly_dts'] = -30.*np.ones(12)
#step_params['pdd_params']['monthly_dps'] = np.ones(12)

wrapper.step(step_params)
dolfin.plot(wrapper.model.H0_c)
plt.show()





