import numpy as np
from dolfin import *
import matplotlib.pyplot as plt
from model_wrapper import ModelWrapper
from pdd_model.pdd_model import PDDModel
from hydro_model.hydro_model import HydroModel

set_log_level(40)

class PDDWrapper(ModelWrapper):

    def __init__(self, inputs):
        
        super(PDDWrapper, self).__init__(inputs)

        ### Initialize the PDD model
        ########################################################

        pdd_params = {}
        if 'pdd_params' in inputs:
            pdd_params = inputs['pdd_params']

        self.pdd_model = PDDModel(self, pdd_params) 
        

        ### Initialize the hydrology model
        ########################################################

        hydro_params = {}
        if 'hydro_params' in inputs:
            hydro_params = inputs['hydro_params']

        self.hydro_model = HydroModel(self, hydro_params)

        ### Make sure everything is initialized
        ########################################################
        
        self.update_inputs(float(self.input_functions['L0']), {})

        
    # Assign model input functions
    def update_inputs(self, L, params):
 
        # Update length dependent fields
        self.update_interp_all(L)

        # Ice model parameters
        ice_params = {}
        if 'ice_params' in params:
            ice_params = params['ice_params']
        
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
        
        # Update the surface mass balance function 
        self.pdd_model.update(pdd_params)
        # Update the effective pressure
        self.hydro_model.update(hydro_params)
        # Update model inputs like sea level
        self.model.update(ice_params)
        

    # Step the model forward by one time step
    def step(self, params = {}):
        # Update the model inputs
        self.update_inputs(float(self.model.L0), params)
        # Get step
        step_params = {}
        if 'step_params' in params:
            step_params = params['step_params']
        L = self.model.step(step_params)
        print(L)
