import numpy as np
from dolfin import *
from model.smb_model.pdd_calculator import PDDCalculator
from model.smb_model.smb_params import smb_params
import matplotlib.pyplot as plt

"""
Positive degree day model. Computes a surface mass balance function adot. 
"""

class SMBModel(object):

    def __init__(self, model_wrapper, params= {}):
        self.model_wrapper = model_wrapper
        self.smb_params = smb_params
        self.smb_params.update(params)
        self.fields = ['smb', 'submarine_melt_scale']
        model_wrapper.load_fields(self.smb_params['fields'], self.fields)
        

    """
    Update the model SMB function.
    """
    def update(self, params):
        
        # Constant offset
        scale = 1.
        if 'scale' in params:
            scale = params['scale']


        submarine_melt_rate = self.smb_params['submarine_melt_rate']

        submarine_melt_scale = self.model_wrapper.input_functions['submarine_melt_scale'].vector().get_local()[0]
        submarine_melt_rate *= submarine_melt_scale

        B = self.model_wrapper.ice_model.B.vector().get_local()
        Bhat = project(self.model_wrapper.ice_model.Bhat).vector().get_local()
        submarine_melt_rate = ((Bhat - B) > 25.)*submarine_melt_rate
        smb = scale*self.model_wrapper.input_functions['smb'].vector().get_local()

        self.model_wrapper.ice_model.adot.vector()[:] = smb + submarine_melt_rate
