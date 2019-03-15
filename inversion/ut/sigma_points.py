import numpy as np

"""
This object generates sigma points and weight sets given a prior
mean x and covariance matrix Pxx.
"""

class SigmaPoints(object):

    def __init__(self, x, Pxx):
         # Sample mean
         self.x = x
         # Sample covariance
         self.Pxx = Pxx
         # State dimension
         self.n = len(x)
         # Matrix square root of Pxx
         self.Pxx_sqrt = np.linalg.cholesky(Pxx)

         
    # Get a specific sigma point set
    def get_set(self, set_type, kappa):

        if set_type == 'min':
            return self.__get_set_min__(kappa)
        elif set_type == 'random':
            return self.__get_set_random__(kappa)
        elif set_type == 'fifth_order':
            return self.__get_set_fifth_order__(kappa)
        elif set_type == 'mean':
            return self.__get_mean_set__(kappa)

        return self.__get_set_gauss__(kappa)


    # Return a mean set of 2n + 1 points
    def __get_mean_set__(self, w_0):

        # Dimension
        n = self.n
        x = self.x
        Pxx_sqrt = self.Pxx_sqrt
        
        # Number os sigma pionts
        N = 2*n + 1
        # Rows are sigma points
        X = x[:,None].repeat(N, axis = 1).T
        X[1:(n+1),:] += np.sqrt(n / (1. - w_0))*Pxx_sqrt
        X[(n+1): ,:] -= np.sqrt(n / (1. - w_0))*Pxx_sqrt
        # Mean / covariance weights
        w = (np.ones(N) - w_0) / (2.*n)
        w[0] = w_0
        
        return X, w, w
        
        
    # Return a Gauss set with 2n + 1 points
    def __get_set_gauss__(self, kappa):

        n = self.n
        x = self.x
        Pxx_sqrt = self.Pxx_sqrt
        
        # Dimension
        n = len(x)
        # Spread parameter
        w_0 = 1. - (n / kappa)
        # Number of sigma pionts
        N = 2*n + 1
        # Rows are sigma points
        X = x[:,None].repeat(N, axis = 1).T
        X[1:(n+1),:] += np.sqrt(n / (1. - w_0))*Pxx_sqrt
        X[(n+1): ,:] -= np.sqrt(n / (1. - w_0))*Pxx_sqrt
        # Mean / covariance weights
        w = np.ones(N) / (2.*kappa)
        w[0] = w_0
        
        return X, w, w
        

    # Return a minimal set of n + 1 points
    def __get_set_min__(self, w0 = 0.5):
        
        n = self.n
        x = self.x
        Pxx_sqrt = self.Pxx_sqrt

        alpha = np.sqrt((1. - w0) / n)
        C = np.linalg.cholesky(np.diag(np.ones(n), 0) - (alpha**2)*np.ones((n, n)))
        W = np.diag(np.diag(np.linalg.multi_dot([w0*(alpha**2)*np.linalg.inv(C), np.ones((n,n)), np.linalg.inv(C.T)])), 0)
        W_sqrt = np.linalg.cholesky(W)
        X = np.zeros((n, n+1))
        X[:,0] = np.dot(-Pxx_sqrt, (alpha / np.sqrt(w0))*np.ones(n))
        X[:,1:] = np.linalg.multi_dot([Pxx_sqrt, C, np.linalg.inv(W_sqrt)])
        X = X.T
        X += x

        # Array of weights
        w = np.zeros(n+1)
        w[0] = w0
        w[1:] = np.diag(W, 0)

        return X, w, w


    # Return a fifth order set with 2*n^2 + 1 points
    def __get_set_fifth_order__(self, r = np.sqrt(3.)):

        # Dimension
        N = self.n
        x = self.x

        ### Generate Weights
        ########################################################

        # Coordinate for the first symmetric set
        r1 = (r*np.sqrt(N-4.))/np.sqrt(N - r**2 - 1.)
        # First symmetric set weight
        w2 = (4. - N) / (2. * r1**4)
        # Second symmetric set weight
        w3 = 1. / (4. * r**4)
        # Center point weight
        w1 = 1. - 2.*N*w2 - 2.*N*(N-1)*w3
        # Vector of weights
        w = np.block([w1, np.repeat(w2, 2*N), np.repeat(w3, 2*N*(N-1))])


        ### Generate Points
        ########################################################
        
        # First fully symmetric set
        X0 = r1*np.eye(N)
        X0_s = np.block([X0, -X0])
        
        # Second fully symmetric set
        X1 = r*np.eye(N)
        indexes_i = []
        indexes_j = []
        for i in range(1,N):
            indexes_i.append(np.repeat([i],i))
            indexes_j.append(np.arange(0,i))
        indexes_i = np.concatenate(indexes_i).ravel()
        indexes_j = np.concatenate(indexes_j).ravel()
        P1 = X1[indexes_i, :].T + X1[indexes_j, :].T
        P2 = X1[indexes_i, :].T - X1[indexes_j, :].T
        X1_s = np.block([P1, P2, -P1, -P2])

        # Full set of points (columns are points)
        X = np.block([np.zeros(N)[:,None], X0_s, X1_s])

        # Change variables
        X = x[:,None].repeat(2*N**2 + 1, axis = 1) + self.Pxx_sqrt@X

        return X.T, w, w


    # Just take random draws from the distribution (for testing)
    def __get_set_random__(self, num_draws):
        samples = np.random.multivariate_normal(self.x, self.Pxx, num_draws)
        w = np.ones(num_draws)
        return samples, w, w
        
        
