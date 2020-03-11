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
        self.smb_params = params
        self.fields = ['smb']
        model_wrapper.load_fields(self.smb_params['fields'], self.fields)
        self.adot = model_wrapper.input_functions['smb']
        

    """
    Update the model SMB function.
    """
    def update(self, params):
        #print(params)
        self.model_wrapper.ice_model.adot.vector()[:] = self.adot.vector().get_local()
