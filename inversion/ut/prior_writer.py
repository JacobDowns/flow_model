import numpy as np
from scipy.interpolate import interp1d
from ut.sigma_points import SigmaPoints
import matplotlib.pyplot as plt

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

        
        ### Plot samples from prior
        ##########################################################################

        samples = np.random.multivariate_normal(x, Pxx, 100)

        for i in range(100):
            plt.plot(samples[i])
        plt.show()
        
         
        ### Compute a minimal set of sigma points
        ##############################################

        points = SigmaPoints(x, Pxx)
        X, wm, wc = points.get_set(set_type, inputs['kappa'])


        for i in range(len(X)):
            x_i = X[i]
            plt.plot(x_i[0:3], marker = 'o')

        plt.show()


        branches = [[ np.array([v]) for v in np.unique(X[:,0]) ]]


        branches = []
    
        for i in range(1,X.shape[1]):
            x = len(np.unique(X[:, 0:i], axis = 0))
            #print(x)
            branches.append(x)


        #plt.plot(branches)
        #plt.show()
        #quit()

        branches = np.array(branches)
        print()
        print(len(branches))
        print(branches)
        print(X.shape)
        print(branches.sum() / (X.shape[0]*X.shape[1]))
        quit()
    
        #print(X)
        #quit()
        x_n = np.zeros(N)
        P_n = np.zeros((N,N))
        # Plot the sigma points and check the mean and covariance
        for i in range(len(X)):
            wm_i = wm[i]
            wc_i = wc[i]
            x_i = X[i]
            x_n += wm_i*x_i
            P_n += wc_i*np.outer(x_i - x, x_i - x)
            plt.plot(x_i)

        plt.show()


        ### Write the prior info to a file
        ####################################################
        np.savetxt(out_dir + 'prior/sigma_ts.txt', sigma_ts)
        np.savetxt(out_dir + 'prior/prior_m.txt', x)
        np.savetxt(out_dir + 'prior/prior_P.txt', Pxx)
        np.savetxt(out_dir + 'prior/w_m.txt', wm)
        np.savetxt(out_dir + 'prior/w_c.txt', wc)
        np.savetxt(out_dir + 'prior/X.txt', X)

        
