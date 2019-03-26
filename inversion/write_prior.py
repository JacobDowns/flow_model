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
n1 = 4
n2 = 22
sigma_ts = np.linspace(-11554., 0., N)
inputs['sigma_ts'] = sigma_ts[n1:n2]
print(inputs['sigma_ts'])
#quit()
# Prior mean 
inputs['x'] = np.loadtxt('paleo_data/south_buizert/opt_m.txt')[n1:n2]
# Prior precision matrix
delta = 7.5e3
Q = delta*np.asarray(poisson((N,)).todense())
# Prior covariance
inputs['Pxx'] = np.linalg.inv(Q)[n1:n2, n1:n2]
# Sigma set type
inputs['set_type'] = 'min'
# The first weight for tuning
inputs['kappa'] = .5


pw = PriorWriter(inputs)
