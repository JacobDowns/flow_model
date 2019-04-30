import numpy as np
import matplotlib.pyplot as plt
from model.transient_wrapper import TransientWrapper

inputs = {}

# Model inputs
inputs['model_inputs'] = {}
model_inputs = inputs['model_inputs']

model_inputs['input_file_name'] = '/home/jake/flow_model/output_files/inversions/covariance_comparison/flowline.h5'
model_inputs['state_file_name'] = '/home/jake/flow_model/output_files/inversions/covariance_comparison/state0.h5'

wrapper = TransientWrapper(inputs)

