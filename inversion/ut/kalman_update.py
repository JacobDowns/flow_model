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


        ### Use a specially built measurement mean and covariance matrix
        #############################################################
        
        # Measurement ages
        self.y_ages = input_dict['y_ages']        
        # Measurement mean
        self.y = input_dict['y']
        dt = int(self.y_ages[1] - self.y_ages[0])
        # Measurement covariance
        self.y = interp1d(input_dict['y_ages'], input_dict['y'])(self.model_ages)
        self.Pyy = np.diag(interp1d(input_dict['y_ages'], np.diag(input_dict['Pyy']))(self.model_ages))
        # Restrict
        self.y = interp1d(input_dict['y_ages'], input_dict['y'])(self.model_ages)[::3*25]
        self.Pyy = np.diag(interp1d(input_dict['y_ages'], np.diag(input_dict['Pyy']))(self.model_ages))[::3*25,::3*25]
        self.Y = self.Y[:,::3*25]

        


    ### Build the joint distribution
    def get_joint_dist(self):
        
        # Compute predicted mean
        mu = np.dot(self.w_m, self.Y)

        plt.plot(mu)
        plt.show()
        #quit()

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
      
        return x_new, Pxx_new, mu


    # Do the Kalman update to incorporate the measurement and correct the prior mean
    def optimize(self, out_dir = None):

        ### Do the Kalman update
        #############################################################

        x_new, Pxx_new, mu = self.get_conditional_dist()  
        
        if out_dir:
            np.savetxt(out_dir + 'mu.txt', mu)
            np.savetxt(out_dir + 'x_new.txt', x_new)
            np.savetxt(out_dir + 'Pxx_new.txt', Pxx_new)
            np.savetxt(out_dir + 'y.txt', self.y)
            np.savetxt(out_dir + 'Py.txt', self.Pyy)

            v = np.sqrt(np.diag(Pxx_new))

            print(x_new[18:])

            plt.plot(x_new[0:18])
            plt.plot((x_new + 2.0*v)[0:18])
            plt.plot((x_new - 2.0*v)[0:18])
            plt.show()
        
