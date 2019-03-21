from ut.prior_writer import PriorWriter
import numpy as np
import sys
from pyamg.gallery import poisson

inputs = {}

# Directory to write prior
inputs['out_dir'] = sys.argv[1]
# Sigma point times
N = 44
# For shorter inversions
n = 24
sigma_ts = np.linspace(-11554., 0., N)
inputs['sigma_ts'] = sigma_ts[0:n]
# Prior mean 
inputs['x'] = np.loadtxt('paleo_data/south_buizert/opt_m.txt')[0:n]
# Prior precision matrix
delta = 7.5e3
Q = delta*np.asarray(poisson((N,)).todense())
# Prior covariance
inputs['Pxx'] = np.linalg.inv(Q)[0:n, 0:n]
# Sigma set type
inputs['set_type'] = 'fifth_order'
# The first weight for tuning
inputs['kappa'] = 1.

pw = PriorWriter(inputs)
