import numpy as np
from model.transient_runner import *
set_log_active(False)

""" 
Example model run that starts at 12 ka BP and runs for 1000 years.
"""

### Model inputs
#######################################################

# Input dictionary
inputs = {}
# Flowline inputs file
inputs['in_file'] = 'input_files/flowline.h5'
# Initial states input file
inputs['state_file_name'] = 'input_files/state1.h5'
# Time step
inputs['dt'] = 1./3.
# Number of time steps
inputs['N'] = 500*3
# Start age
inputs['start_age'] = -11.6e3
# Return ice thickness every 50 time steps
inputs['snapshot_interval'] = 1000


### Perform the model run
#######################################################
tr = TransientRunner(inputs)
ages, Ls, Hs, Ps, adots = tr.run()

out_dir = 'output_files/'
np.savetxt(out_dir + 'age.txt', ages)
np.savetxt(out_dir + 'L.txt', Ls)
np.savetxt(out_dir + 'H.txt', Hs)
np.savetxt(out_dir + 'P.txt', Ps)
