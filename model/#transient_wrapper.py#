import numpy as np
import matplotlib.pyplot as plt
from model_wrapper import ModelWrapper
from pdd_model.pdd_model import PDDModel

class TransientWrapper(ModelWrapper):

    def __init__(self, inputs):
        
        super(TransientWrapper, self).__init__(inputs)

        ### Initialize the PDD model
        ########################################################

        pdd_params = {}
        if 'pdd_params' in inputs:
            pdd_params = inputs['pdd_params']

        self.pdd_model = PDDModel(pdd_params) 
        


      # Assign model input functions
      def update_inputs(self, L, params):

        
          
        # Update length dependent fields
        self.update_interp_all(L)

        # Update model bed elevation
        self.B.assign(self.input_functions['B'])
        # Update model basal traction
        self.beta2.assign(self.model_inputs.input_functions['beta2'])
        # Update model surface
        self.model.S0_c.assign(self.model.B + self.model.H0_c)
        # Update model width
        self.width.assign(self.input_functions['width'])
        
        # Update the surface mass balance
        self.model.adot.assign(self.pdd_model.get_adot(params))
    
        self.precip_func.assign(self.model_inputs.precip_func)
