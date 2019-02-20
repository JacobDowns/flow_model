import numpy as np
from model.transient_runner import *

""" 
Example model run that starts at 12 ka BP and runs for 1000 years.
"""

### Model inputs
#######################################################

# Input dictionary
inputs = {}
# Input file name
inputs['in_file'] = 'input_files/steady.h5'
# Time step
inputs['dt'] = 1./3.
# Number of time steps
inputs['N'] = 100*3
# Start age
inputs['start_age'] = -11.6e3
# Return ice thickness every 50 time steps
inputs['snapshot_interval'] = 100


### Perform the model run
#######################################################
tr = TransientRunner(inputs)
ages, Ls, Hs, Ps = tr.run()


### Plot
#######################################################
import matplotlib.pyplot as plt

out_dir = 'output_files/'
np.savetxt(out_dir + 'age.txt', ages)
np.savetxt(out_dir + 'L.txt', Ls)
np.savetxt(out_dir + 'H.txt', Hs)
np.savetxt(out_dir + 'P.txt', Ps)

plt.plot(ages, Ls)
plt.xlabel('Age (ka BP)')
plt.ylabel('Glacier Length (m)')
plt.show()
