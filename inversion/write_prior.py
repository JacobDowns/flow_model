from ut.prior_writer import PriorWriter
import numpy as np
import sys
from pyamg.gallery import poisson

inputs = {}

# Directory to write prior
inputs['out_dir'] = sys.argv[1]
# Sigma point times
N = 40
sigma_ts = np.linspace(-11554., 0., N)
inputs['sigma_ts'] = sigma_ts
# Prior mean 
chi = np.linspace(0., 1., N)
inputs['x'] = 0.*chi
# Prior precision matrix
delta = 7.5e3
Q = delta*np.asarray(poisson((N,)).todense())
# Prior covariance
inputs['Pxx'] = np.linalg.inv(Q)
# Sigma set type
inputs['set_type'] = 'fifth_order'
# The first weight for tuning
inputs['kappa'] = 3.

pw = PriorWriter(inputs)
