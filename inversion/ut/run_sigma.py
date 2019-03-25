import sys
sys.path.insert(0, '/home/jake/flow_model/')
from model.transient_runner import *
import numpy as np
import copy
import os
set_log_active(False)

# Input directory
in_dir = sys.argv[1]
# Sigma index
index = int(sys.argv[2]) - 1
# Load sigma points
sigma_points = np.loadtxt(in_dir + 'prior/X.txt')
num_sigma_points = sigma_points.shape[0]
# Load sigma times
sigma_ts = np.array(np.loadtxt(in_dir + 'prior/sigma_ts.txt'), dtype = int)
# Model input file
input_file = in_dir + 'flowline.h5'
state_file = in_dir + 'state0.h5'

### Model inputs
#######################################################

# Input dictionary
inputs = {}
# Input file name
inputs['in_file'] = input_file
# Initial state
inputs['state_file_name'] = state_file
# Time step
inputs['dt'] = 1./3.
# Number of model time steps
inputs['N'] = abs(sigma_ts.max() - sigma_ts.min())*3


### Run a sigma point
#######################################################

# Check if a sigma point has been run yet
if not os.path.isfile(in_dir + 'sigmas/Y_' + str(index) + '.txt'):
    print(index)

    # Interpolated delta temp
    X_i = sigma_points[index]
    inputs['precip_param_func'] = interp1d(sigma_ts, X_i, kind = 'linear')
    #inputs['start_age'] = sigma_ts[0]
    

    ### Perform model run 
    #######################################################

    model_runner = TransientRunner(inputs)
    ages, Ls, Hs, Bs, Ps, adots = model_runner.run()

    np.savetxt(in_dir + 'sigmas/age_' + str(index) + '.txt', ages)
    np.savetxt(in_dir + 'sigmas/Y_' + str(index) + '.txt', Ls)
