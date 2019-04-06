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
    inputs['out_dir'] = 'output_files/inversions/sensitivity/'    
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
    inputs['Puu'] = np.diag(np.array([0.1**2, 8e-5**2, 0.02**2, 0.001**2, 0.001**2]))
    # Steady state terminus position
    inputs['L_steady'] = np.array([402078.])

if sys.argv[1] == 'sensitivity1':
    # Directory to write prior
    inputs['out_dir'] = 'output_files/inversions/sensitivity1/'    
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
    inputs['kappa'] = 1.5

    ### Sensitivity params
    ####################################################

    # Parameter names
    inputs['param_names'] = np.array(['log_A', 'beta2', 'lambda_precip', 'lambda_ice', 'lambda_snow'])
    # Parameter mean vector
    inputs['u'] = np.array([np.log(3.5e-25), 1.6e-3, 0.07, 0.008, 0.005])
    # Parameter covariance
    inputs['Puu'] = np.diag(np.array([0.1**2, 8e-5**2, 0.02**2, 0.001**2, 0.001**2]))
    # Steady state terminus position
    inputs['L_steady'] = np.array([402078.])


    
if sys.argv[1] == 'sensitivity2':
    # Directory to write prior
    inputs['out_dir'] = 'output_files/inversions/sensitivity2/'    
    # For shorter inversions
    n1 = 4
    n2 = 22
    inputs['sigma_ts'] = sigma_ts[n1:n2]
    # Prior mean 
    #inputs['x'] = np.loadtxt('paleo_data/south_buizert/opt_m.txt')[n1:n2]
    # Prior mean 
    inputs['x'] = np.loadtxt('output_files/inversions/sensitivity/posterior/x_new.txt')[0:(n2 - n1)]
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
    inputs['u'] = np.array([-5.62e+01, 1.58e-03, 1.10292653e-01, 8.58e-03, 5.49e-3])
    # Parameter covariance
    inputs['Puu'] = np.diag(np.array([0.1**2, 8e-5**2, 0.02**2, 0.001**2, 0.001**2]))
    # Steady state terminus position
    inputs['L_steady'] = np.array([402078.])

if sys.argv[1] == 'sensitivity3':
    # Directory to write prior
    inputs['out_dir'] = 'output_files/inversions/sensitivity3/'    
    # For shorter inversions
    n1 = 4
    n2 = 22
    inputs['sigma_ts'] = sigma_ts[n1:n2]
    # Prior mean 
    inputs['x'] = np.loadtxt('output_files/inversions/sensitivity1/posterior/x_new.txt')[0:(n2 - n1)]
    # Prior covariance
    inputs['Pxx'] = np.linalg.inv(Q)[n1:n2, n1:n2]
    # Sigma set type
    inputs['set_type'] = 'fifth_order'
    # The first weight for tuning
    inputs['kappa'] = 1.2

    ### Sensitivity params
    ####################################################

    # Parameter names
    inputs['param_names'] = np.array(['log_A', 'beta2', 'lambda_precip', 'lambda_ice', 'lambda_snow'])
    # Parameter mean vector
    inputs['u'] = np.array([-5.64433066e+01, 1.71693691e-03, 1.42872219e-01, 5.73006099e-03, 8.13594366e-03])
    # Parameter covariance
    inputs['Puu'] = np.diag(np.array([0.1**2, 8e-5**2, 0.02**2, 0.001**2, 0.001**2]))
    # Steady state terminus position
    inputs['L_steady'] = np.array([402078.])
    
pw = PriorWriter(inputs)
