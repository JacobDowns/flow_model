import numpy as np
from dolfin import *
import matplotlib.pyplot as plt
from pdd_wrapper import PDDWrapper

class TransientWrapper(PDDWrapper):

    def __init__(self, inputs):
        
        super(TransientWrapper, self).__init__(inputs)

        transient_inputs = inputs['transient_inputs']

        # The number of model steps
        self.N = transient_inputs['N']
