import numpy as np
import sys
sys.path.insert(0, '/home/jake/flow_model/')
from inversion.ut.kalman_update import KalmanUpdate

# Flowline
flowline = sys.argv[1]
# Input dictionary
inputs = {}
    
### Center
#############################################################

if flowline == 'cov_compare':
    inputs['in_dir'] = 'output_files/inversions/covariance_comparison'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/covariance_comparison/posterior/')

if flowline == 'cov_compare1':
    inputs['in_dir'] = 'output_files/inversions/covariance_comparison1'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/covariance_comparison1/posterior/')

if flowline == 'htm_covariance':
    inputs['in_dir'] = 'output_files/inversions/htm_covariance'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/htm_covariance/posterior/')
