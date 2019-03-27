import numpy as np

"""
Computes sigma points for a scalar Gaussian.
"""

class SigmaPointsScalar(object):

    def __init__(self, alpha, beta, kappa):
        self.alpha = alpha
        self.beta = beta
        self.kappa = kappa
        self.lam = alpha**2 * (1. + kappa) - 1.

        wm0 = self.lam / (1. + self.lam)
        wc0 = wm0 + 1. - alpha**2 + beta
        w1 = 1. / (2. * (1. + self.lam))

        self.mean_weights = np.array([w1, wm0, w1])
        self.variance_weights = np.array([w1, wc0, w1])


    def sigma_points(self, mean, variance):
        points = np.array([mean - np.sqrt((1. + self.lam)*variance), mean, mean + np.sqrt((1. + self.lam)*variance)])
        return points
