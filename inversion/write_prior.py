from ut.prior_writer import PriorWriter
import numpy as np
import sys
from pyamg.gallery import poisson

inputs = {}

# Number of points in time
N = 44
# Precip. anomaly times
sigma_ts = np.linspace(-11554., 0., N)
# Prior precision matrix
delta = 7.5e3
Q = delta*np.asarray(poisson((N,)).todense())

if sys.argv[1] == 'cov1':
    # Directory to write prior
    inputs['out_dir'] = 'output_files/inversions/cov1/'    
    # For shorter inversions
    n1 = 4
    n2 = 22
    inputs['sigma_ts'] = sigma_ts[n1:n2]
    # Prior mean 
    inputs['x'] = np.loadtxt('paleo_data/south_buizert/opt_m.txt')[n1:n2]
    # Prior covariance
    inputs['Pxx'] = np.linalg.inv(Q)[n1:n2, n1:n2]
    # Sigma set type
    inputs['set_type'] = 'fifth_order'
    # The first weight for tuning
    inputs['kappa'] = 2.9

if sys.argv[1] == 'cov2':
    # Directory to write prior
    inputs['out_dir'] = 'output_files/inversions/cov2/'    
    # For shorter inversions
    n1 = 4
    n2 = 22
    inputs['sigma_ts'] = sigma_ts[n1:n2]
    # Prior mean 
    inputs['x'] = np.loadtxt('paleo_data/south_buizert/opt_m.txt')[n1:n2]
    # Prior covariance
    inputs['Pxx'] = np.linalg.inv(Q)[n1:n2, n1:n2]
    # Sigma set type
    inputs['set_type'] = 'min'
    # The first weight for tuning
    inputs['kappa'] = .5

if sys.argv[1] == 'cov3':
    # Directory to write prior
    inputs['out_dir'] = 'output_files/inversions/cov3/'    
    # For shorter inversions
    n1 = 4
    n2 = 22
    inputs['sigma_ts'] = sigma_ts[n1:n2]
    # Prior mean 
    inputs['x'] = np.loadtxt('paleo_data/south_buizert/opt_m.txt')[n1:n2]
    # Prior covariance
    inputs['Pxx'] = np.linalg.inv(Q)[n1:n2, n1:n2]
    # Sigma set type
    inputs['set_type'] = 'gauss'
    # The first weight for tuning
    inputs['kappa'] = 20.

if sys.argv[1] == 'sensitivity':
    # Directory to write prior
    inputs['out_dir'] = 'output_files/inversions/cov3/'    
    # For shorter inversions
    n1 = 4
    n2 = 22
    inputs['sigma_ts'] = sigma_ts[n1:n2]
    # Prior mean 
    inputs['x'] = np.loadtxt('paleo_data/south_buizert/opt_m.txt')[n1:n2]
    # Prior covariance
    inputs['Pxx'] = np.linalg.inv(Q)[n1:n2, n1:n2]
    # Sigma set type
    inputs['set_type'] = 'min'
    # The first weight for tuning
    inputs['kappa'] = 0.5

    ### Sensitivity params
    ####################################################

    # Parameter names
    inputs['param_names'] = np.array(['log_A', 'beta2', 'lambda_precip', 'lambda_ice', 'lambda_snow'])
    # Parameter mean vector
    inputs['u'] = np.array([np.log(3.5e-25), 1.6e-3, 0.07, 0.008, 0.005])
    # Parameter covariance
    inputs['Puu'] = np.diag(np.array([0.25**2, 8e-5**2, 0.02**2, 0.001**2, 0.001**2]))
    
pw = PriorWriter(inputs)
