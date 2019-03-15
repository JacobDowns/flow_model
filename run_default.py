import numpy as np
from model.transient_runner import *
set_log_active(False)

### Model inputs
#######################################################

# Input dictionary
inputs = {}
# Flowline inputs file
inputs['in_file'] = 'input_files/south_buizert/flowline.h5'
# Initial states input file
inputs['state_file_name'] = 'input_files/south_buizert/state0.h5'
# Time step
inputs['dt'] = 1./3.
# Return ice thickness every 10 years
inputs['snapshot_interval'] = 15*3
# Basal traction
inputs['beta2'] = 1.6e-3


### Delta P
#######################################################
# State vector times
sigma_ts = np.loadtxt( 'paleo_data/south_buizert/sigma_ts.txt')
precip_param_opt = np.loadtxt('paleo_data/south_buizert/opt_m.txt')
# Interpolated delta temp. function 
inputs['precip_param_func'] = interp1d(sigma_ts, precip_param_opt, kind = 'linear')
# Number of model time steps
inputs['N'] = 5000*3 
# Start age
inputs['start_age'] = sigma_ts[0]


### Perform the model run
#######################################################
tr = TransientRunner(inputs)
ages, Ls, Hs, Bs, Ps, adots = tr.run()

out_dir = 'output_files/south_buizert/'
np.savetxt(out_dir + 'age.txt', ages)
np.savetxt(out_dir + 'L.txt', Ls)
np.savetxt(out_dir + 'H.txt', Hs)
np.savetxt(out_dir + 'B.txt', Bs)
np.savetxt(out_dir + 'P.txt', Ps)
np.savetxt(out_dir + 'adot.txt', adots)
