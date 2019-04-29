import numpy as np
import matplotlib.pyplot as plt
from model_wrapper import ModelWrapper

inputs = {}

# Domain length
domain_len = 450e3
inputs['domain_len'] = domain_len
# Along flow coordinate
x = np.linspace(0., domain_len, 450.) 
# Initial glacier length
L0 = 400e3
inputs['L0'] = L0
# Maximum ice thickness
H_max = 3000.

# Bed
inputs['B'] = np.zeros_like(x)
# Surface
#inputs['S_ref'] = np.zeros_like(x)
#indexes = x <= L0
#inputs['S_ref'][indexes] = np.sqrt((H_max)**2*(1. - (x[indexes] / L0)))
# Width
#inputs['width'] = 1000.*np.ones_like(x)
# Temperature and precip

for i in range(12):
    inputs['T' + str(i)] = -10. + x / 100e3
    inputs['P' + str(i)] = np.sin(x / 100e3)**2
    #plt.plot(inputs['P' + str(i)])

#plt.plot(S)
#plt.show()

inputs['input_file_name'] = '/home/jake/flow_model/output_files/inversions/covariance_comparison1/flowline.h5'
inputs['state_file_name'] = '/home/jake/flow_model/output_files/inversions/covariance_comparison1/state0.h5'

wrapper = ModelWrapper(inputs)

