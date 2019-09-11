import numpy as np
from scipy.special import erf, erfc

class PDDCalculator():

    def __init__(self, sigma):
        # Standard deviation for PDD model
        self.pdd_sigma = sigma
        # For unit conversion
        self.month = 365. / 12.
        # Standard deviation for computing accumulation
        self.acc_sigma = sigma - 0.5


    ### Inner iterated integral for computing pdds
    def __inner_int_pdd__(self, T_m, sigma):
        sigma2 = sigma**2
        val  = sigma2*(np.e**(-T_m**2 / (2.*sigma2)))
        val += np.sqrt(np.pi/2.)*T_m*sigma*(erf(T_m / (np.sqrt(2.)*sigma)) + 1.)
        return val


    ### Inner iterated integral for computing accumulation fraction
    def __inner_int_acc__(self, T_m, sigma):
        sigma2 = sigma**2
        c = np.sqrt(1. / sigma2)
        val = (np.sqrt(np.pi / 2.) / c) * erfc(T_m * (c / np.sqrt(2.)))
        return val


    ### Calculate PDDs given mean monthly temperature
    def __get_pdd__(self, T_m, sigma):
        sigma2 = sigma**2
        return (1./(sigma*np.sqrt(2.*np.pi))) * self.__inner_int_pdd__(T_m, sigma) * self.month


    ### Compute PDD's from mean monthly tempearture T_m (units degrees C * days)
    def get_pdd(self, T_m):
        return self.__get_pdd__(T_m, self.pdd_sigma)


    ### Get the fraction of accumulation that falls as snow
    def __get_acc_frac__(self, T_m, sigma):
        return (1./(sigma*np.sqrt(2.*np.pi))) * self.__inner_int_acc__(T_m, sigma)


    ### Get the fraction of accumulation that falls as snow
    def get_acc_frac(self, T_m):
        return self.__get_acc_frac__(T_m, self.acc_sigma)
