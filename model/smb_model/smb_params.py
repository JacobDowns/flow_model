import numpy as np
from util.custom_dict import CustomDict

smb_params = CustomDict()
# PDD model variance 
smb_params['pdd_var'] = 5.5
# Ablation rate for snow (m / (degree C * day))
smb_params['lambda_snow'] = 0.005
# Ablation rate for ice (m / (degree C * day))
smb_params['lambda_ice'] = 0.008
# Precipitation scaling parameter
smb_params['lambda_precip'] = 0.07
# Superimposed ice fraction
smb_params['super_ice_frac'] = 0.6
# Elevation lapse rate  (degrees C / km)
smb_params['lapse_rate'] = 5.
# Monthly temp. anomalies (C)
smb_params['monthly_dts'] = np.zeros(12)
# Monthly precip. anomalies (m.w.e. / a)
smb_params['monthly_dps'] = np.zeros(12)
