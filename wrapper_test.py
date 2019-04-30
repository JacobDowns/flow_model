import numpy as np
import matplotlib.pyplot as plt
from model.pdd_wrapper import PDDWrapper

inputs = {}

# Model inputs
inputs['model_inputs'] = {}
model_inputs = inputs['model_inputs']

model_inputs['input_file_name'] = '/home/jake/flow_model/output_files/inversions/cov1/flowline.h5'
model_inputs['state_file_name'] = '/home/jake/flow_model/output_files/inversions/cov1/state0.h5'
model_inputs['dt'] = 1./3.

wrapper = PDDWrapper(inputs)
wrapper.step()


