from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from data_extractor import DataExtractor
import sys

""" 
Extract bed, surface, precip., and temp., data along a flowline specified by lat, lon coords.
"""

flowline = sys.argv[1]

if flowline == 'north_seasonal':
    flowline_data = np.loadtxt('/home/jake/flow_model/paleo_data/flowlines/north_seasonal_coords.txt')
    output_file = 'h5_files/north_seasonal.h5'
    
lons = flowline_data[:,0][::-1]
lats = flowline_data[:,1][::-1]

temp_ext = DataExtractor()
temp_ext.write_mesh(output_file, lats, lons)
