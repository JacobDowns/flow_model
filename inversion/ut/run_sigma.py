import sys
sys.path.insert(0, '/home/jake/flow_model/')
from model.transient_runner import *
import numpy as np
import copy
import os
import matplotlib.pyplot as plt
set_log_active(False)

# Input directory
in_dir = sys.argv[1]
# Sigma index
index = int(sys.argv[2]) - 1
# Load sigma points
X = np.loadtxt(in_dir + 'prior/X.txt')
num_sigma_points = X.shape[0]
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
# Time step
inputs['dt'] = 1./3.
# Number of model time steps
inputs['N'] = abs(sigma_ts.max() - sigma_ts.min())*3
# Start time
#inputs['start_age'] = sigma_ts[0]


### Sensitivity run
#######################################################

x_i = X[index]
if os.path.isfile(in_dir + 'prior/param_names.txt'):
    # Load parameter names
    param_names = np.loadtxt(in_dir + 'prior/param_names.txt', dtype = str)
    # Restrict the sigma points to parameters only and get unique param. sets
    U = np.unique(X[:, -len(param_names):], axis = 0)
    # Figure out what paramter set this sigma point uses
    u_i = x_i[-len(param_names):]
    steady_index = np.linalg.norm(U - u_i, axis = 1).argmin()
    #print(np.linalg.norm(U - u_i, axis = 1))
    state_file = in_dir + 'steady/' + str(steady_index) + '.h5'
    print(state_file)
    # Set sensitivity params. 
    for j in range(len(param_names)):
        inputs[param_names[j]] = u_i[j]
    x_i = x_i[0:-len(param_names)]

inputs['state_file_name'] = state_file

### Run a sigma point
#######################################################

# Check if a sigma point has been run yet
if not os.path.isfile(in_dir + 'sigmas/Y_' + str(index) + '.txt'):
    print(index)

    # Interpolated delta temp
    inputs['precip_param_func'] = interp1d(sigma_ts, x_i, kind = 'linear')
    

    ### Perform model run 
    #######################################################

    model_runner = TransientRunner(inputs)
    ages, Ls, Hs, Bs, Ps, adots = model_runner.run()

    np.savetxt(in_dir + 'sigmas/age_' + str(index) + '.txt', ages)
    np.savetxt(in_dir + 'sigmas/Y_' + str(index) + '.txt', Ls)
