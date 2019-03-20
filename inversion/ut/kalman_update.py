import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class KalmanUpdate(object):

    def __init__(self, input_dict):

        ### Load stuff we need for the unscented transform 
        #############################################################

        # Input dictionary
        self.input_dict = input_dict
        # Input directory 
        self.in_dir = input_dict['in_dir']
        # Model time steps 
        self.model_ages = np.loadtxt(self.in_dir + '/sigmas/age_1.txt')
        # Sigma points
        self.X = np.loadtxt(self.in_dir + '/prior/X.txt')
        # Transformed sigma points
        self.Y = np.loadtxt(self.in_dir + '/sigmas/Y.txt')
        # Prior mean
        self.x = np.loadtxt(self.in_dir + '/prior/x.txt')
        # Prior covariance 
        self.Pxx = np.loadtxt(self.in_dir + '/prior/Pxx.txt')
        # Load mean weights
        self.w_m = np.loadtxt(self.in_dir + '/prior/w_m.txt')
        # Load covariance weights
        self.w_c = np.loadtxt(self.in_dir + '/prior/w_c.txt')

        #print(self.model_ages)
        #quit()

        ### Use a specially built measurement mean and covariance matrix
        #############################################################
        
        # Measurement ages
        self.y_ages = input_dict['y_ages']        
        # Measurement mean
        self.y = input_dict['y']
        dt = int(self.y_ages[1] - self.y_ages[0])
        # Measurement covariance
        self.Pyy = input_dict['Pyy']
        # Measurement indexes
        self.meas_indexes = list(range(0, len(self.model_ages), dt*3))
        # Restrict transformed sigma points
        self.Y = self.Y[:,self.meas_indexes]
        # Restrict the measurement
        self.y = self.y[0:len(self.meas_indexes)]
        self.Pyy = self.Pyy[0:len(self.meas_indexes), 0:len(self.meas_indexes)]


    ### Build the joint distribution
    def get_joint_dist(self):
        
        # Compute predicted mean
        mu = np.dot(self.w_m, self.Y)

        # Compute predicted measurement covariance
        S = np.zeros((self.Y.shape[1], self.Y.shape[1]))
        for i in range(len(self.w_c)):
            print(i)
            S += self.w_c[i]*np.outer(self.Y[i] - mu, self.Y[i] - mu)

        S += self.Pyy
        
        # Compute predicted measurement covariance
        C = np.zeros((self.X.shape[1], self.Y.shape[1]))
        for i in range(len(self.w_c)):
            print(i)
            C += self.w_c[i]*np.outer(self.X[i] - self.x, self.Y[i] - mu)

        x_full = np.block([self.x, mu])
        P_full = np.block([[self.Pxx, C], [C.T, S]])
        
        return x_full, P_full


    ### Compute a conditional probability from the joint probability
    def get_conditional_dist(self):

        ### Partition the joint distribution
        ######################################################################
        
        x_full, P_full = self.get_joint_dist()
        # Length of state variable             
        n1 = len(x_full) - len(self.y)
        # Length of measurement vector
        n2 = len(self.y)
        # State mean
        x = x_full[0:n1]
        # Measurement mean
        mu = x_full[n1:]
        # State prior
        P = P_full[0:n1, 0:n1]
        # Measurement covariance
        S = P_full[n1:, n1:]
        # Cross covariance 
        C = P_full[n1:, 0:n1]

        ### Do the Kalman update
        ######################################################################
        
        K = np.dot(C.T, np.linalg.inv(S))
        x_new = x + np.dot(K, self.y - mu)
        Pxx_new = P - np.dot(np.dot(K, S), K.T)

        #kappa = 9 - n
        #print(kappa)
        #quit()
        v = np.sqrt(np.diag(Pxx_new))
        plt.plot(x_new - 2.*v)
        plt.plot(x_new)
        plt.plot(x_new + 2.*v)
        plt.show()
        quit()

        return x_new, Pxx_new, mu


    # Do the Kalman update to incorporate the measurement and correct the prior mean
    def optimize(self, out_dir = None):

        ### Do the Kalman update
        #############################################################

        self.get_conditional_dist()
    
        
        if out_dir:
            np.savetxt(out_dir + 'mu.txt', mu)
            np.savetxt(out_dir + 'opt_m.txt', m_p)
            np.savetxt(out_dir + 'opt_P.txt', P_p)
            np.savetxt(out_dir + 'y.txt', self.meas_y)
            np.savetxt(out_dir + 'Py.txt', self.meas_Py)
            np.savetxt(out_dir + 'v.txt', v)
            np.savetxt(out_dir + 'y_v.txt', y_v)

            plt.plot(m_p)
            plt.plot(m_p + 2.0*np.sqrt(v))
            plt.plot(m_p - 2.0*np.sqrt(v))
            plt.show()

        return m_p, P_p, mu, K, self.meas_y, self.meas_Py
        
