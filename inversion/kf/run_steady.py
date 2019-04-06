import sys
sys.path.insert(0, '/home/jake/flow_model/')
from model.steady_runner import *
import numpy as np
set_log_active(False)

# Input directory
in_dir = sys.argv[1]
# Sigma index
index = int(sys.argv[2]) - 1
# Load sigma points
X = np.loadtxt(in_dir + 'prior/X.txt')
# Load parameter names
param_names = np.loadtxt(in_dir + 'prior/param_names.txt', dtype = str)
# Restrict the sigma points to parameters only and get unique param. sets
X = np.unique(X[:, -len(param_names):], axis = 0)
param_vals = X[index,:]

# Load sigma times
sigma_ts = np.array(np.loadtxt(in_dir + 'prior/sigma_ts.txt'), dtype = int)
# Model input file
input_file = in_dir + 'flowline.h5'
state_file = in_dir + 'state0.h5'
steady_file = in_dir + 'steady/' + str(index)

### Model inputs
#######################################################
inputs = {}
# Time step
inputs['dt'] = 3.
# Number of model time steps
inputs['N'] = 5000
 # Start age
inputs['start_age'] = sigma_ts[0]
# Input file
inputs['in_file'] = input_file
# Initial state file
inputs['state_file_name'] = state_file
# Steady state file
inputs['steady_file_name'] = steady_file
# Desired terminus position
inputs['L_steady'] = np.loadtxt(in_dir + 'prior/L_steady.txt')
# Set sensitivity params. 
for j in range(len(param_names)):
    inputs[param_names[j]] = param_vals[j]

### Perform model run 
#######################################################
model_runner = SteadyRunner(inputs)
model_runner.run()

