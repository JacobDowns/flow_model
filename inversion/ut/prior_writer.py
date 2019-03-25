import numpy as np
from ut.sigma_points import SigmaPoints
import matplotlib.pyplot as plt
import os

class PriorWriter(object):

    def __init__(self, inputs):

        # Output directory
        out_dir = inputs['out_dir']
        # Sigma point times 
        sigma_ts = inputs['sigma_ts'] 
        # Prior mean
        x = inputs['x']
        # Length of state vector
        N = len(x)
        # Covariance 
        Pxx = inputs['Pxx']
        # Sigma set type
        set_type = inputs['set_type']
        # Parameter names
        if 'param_names' in inputs:
            param_names = inputs['param_names']
            # Parameter mean vector
            u = inputs['u']
            # Param covariance
            Puu = inputs['Puu']

        
        ### Plot samples from prior
        ##########################################################################

        samples = np.random.multivariate_normal(x, Pxx, 100)

        for i in range(100):
            plt.plot(samples[i])
        plt.show()
        
         
        ### Get a set of sigma points
        ##############################################

        points = SigmaPoints(x, Pxx)
        X, wm, wc = points.get_set(set_type, inputs['kappa'])

        # Plot the sigma points
        for i in range(len(X)):
            x_i = X[i]
            plt.plot(x_i, marker = 'o')

        plt.show()


        ### Write the prior info to a file
        ####################################################

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
            os.mkdir(out_dir + '/prior/')
            
        np.savetxt(out_dir + '/prior/sigma_ts.txt', sigma_ts)
        np.savetxt(out_dir + '/prior/x.txt', x)
        np.savetxt(out_dir + '/prior/Pxx.txt', Pxx)
        np.savetxt(out_dir + '/prior/w_m.txt', wm)
        np.savetxt(out_dir + '/prior/w_c.txt', wc)
        np.savetxt(out_dir + '/prior/X.txt', X)

        if 'param_names' in inputs:
            np.savetxt(out_dir + '/prior/param_names.txt', param_names)
            np.savetxt(out_dir + '/prior/u.txt', u)
            np.savetxt(out_dir + '/prior/Puu.txt', Puu)

        
