from model.common_runner import *
import numpy as np

class TransientRunner(CommonRunner):

    def __init__(self, input_dict):

        super(TransientRunner, self).__init__(input_dict)
        
        ### Transient inputs
        ############################################################################
  
        # Precipitation param as a function of age
        self.precip_param_func = lambda x : 0.0
        if 'precip_param_func' in input_dict:
            self.precip_param_func = input_dict['precip_param_func']
            
        # Snapshot interval : write out thickness vector periodically
        self.snapshot_interval = 1000
        if 'snapshot_interval' in input_dict:
            self.snapshot_interval = input_dict['snapshot_interval']


    # Perform a model run
    def run(self):
        # Length at each time step
        Ls = []
        # Age at each time step 
        ages = []
        # Thickness vectors through time
        Hs = []
        # Integral of accumulation  through time
        Ps = []
        # Integral of smb over time
        adots = []

        for j in range(self.N):
            # Age
            age = self.start_age + self.model.t
            ages.append(age)
            # Precip param.
            precip_param = self.precip_param_func(age)

            L = self.model.step(precip_param, accept = True)
            Ls.append(L)
            Ps.append(assemble(self.model.precip_func*dx)*L)
            adots.append(assemble(self.model.adot_prime_func*dx)*L)

            if j % self.snapshot_interval == 0:
                Hs.append(self.model.H0.vector().get_local())

        return np.array(ages), np.array(Ls), np.array(Hs), np.array(Ps), np.array(adots)
