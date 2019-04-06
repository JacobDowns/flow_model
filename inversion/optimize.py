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

if flowline == 'htm_covariance1':
    inputs['in_dir'] = 'output_files/inversions/htm_covariance1'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/htm_covariance1/posterior/')

if flowline == 'cov1':
    inputs['in_dir'] = 'output_files/inversions/cov1'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/cov1/posterior/')

if flowline == 'cov2':
    inputs['in_dir'] = 'output_files/inversions/cov2'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/cov2/posterior/')

if flowline == 'cov3':
    inputs['in_dir'] = 'output_files/inversions/cov3'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/cov3/posterior/')

if flowline == 'sensitivity':
    inputs['in_dir'] = 'output_files/inversions/sensitivity'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/sensitivity/posterior/')

if flowline == 'sensitivity1':
    inputs['in_dir'] = 'output_files/inversions/sensitivity1'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/sensitivity1/posterior/')

if flowline == 'sensitivity2':
    inputs['in_dir'] = 'output_files/inversions/sensitivity2'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/sensitivity2/posterior/')

if flowline == 'sensitivity3':
    inputs['in_dir'] = 'output_files/inversions/sensitivity3'
    inputs['y_ages'] = np.loadtxt('paleo_data/measurements/y_ages.txt')
    inputs['y'] = np.loadtxt('paleo_data/measurements/y_s.txt')
    inputs['Pyy'] = 1.*np.loadtxt('paleo_data/measurements/Py_s.txt')
    ku = KalmanUpdate(inputs)
    ku.optimize('output_files/inversions/sensitivity3/posterior/')
