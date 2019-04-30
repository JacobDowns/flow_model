import numpy as np
import matplotlib.pyplot as plt
from model_wrapper import ModelWrapper
from pdd_model.pdd_model import PDDModel
from hydro_model.hydro_model import HydroModel

class TransientWrapper(ModelWrapper):

    def __init__(self, inputs):
        
        super(TransientWrapper, self).__init__(inputs)

        ### Initialize the PDD model
        ########################################################

        pdd_params = {}
        if 'pdd_params' in inputs:
            pdd_params = inputs['pdd_params']

        self.pdd_model = PDDModel(pdd_params) 
        

        ### Initialize the hydrology model
        ########################################################

        hdyro_params = {}
        if 'hydro_params' in inputs:
            hydro_params = inputs['hydro_params']

        self.hydro_model = HydroModel(hydro_params)
        

    # Assign model input functions
    def update_inputs(self, L, params):
 
        # Update length dependent fields
        self.update_interp_all(L)

        # Ice model parameters
        ice_model = {}
        if 'ice_params' in params:
            ice_params = 'ice_params'
        
        # PDD model parameters
        pdd_params = {}
        if 'pdd_params' in params:
            pdd_params = params['pdd_params']

        # Hydrology model parameters
        hydro_params = {}
        if 'hydro_params' in params:
            hydro_params = params['hydro_params']
        

        ### Update ice model inputs
        ########################################################
        


        # Update model inputs like sea level
        self.model.update(ice_params)
        # Update the surface mass balance function 
        self.pdd_model.update(pdd_params)
        # Update the effective pressure
        self.hydro_model.update(hydro_params)


    # Step the model forward by one time step
    def step(self, params):
        # Update the model inputs
        self.update_inputs(params)
        self.model.step()


    
        


        
