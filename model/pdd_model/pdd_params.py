import numpy as np

pdd_params = {}
# PDD model variance 
pdd_params['pdd_var'] = 5.5
# Ablation rate for snow (m / (degree C * day))
pdd_params['lambda_snow'] = 0.005
# Ablation rate for ice (m / (degree C * day))
pdd_params['lambda_ice'] = 0.008
# Precipitation scaling parameter
pdd_params['lambda_precip'] = 0.07
# Superimposed ice fraction
pdd_params['super_ice_frac'] = 0.6
# Elevation lapse rate  (degrees C / km)
pdd_params['lapse_rate'] = 5.
# Monthly temp. anomalies (C)
pdd_params['monthly_dts'] = np.zeros(12)
# Monthly precip. anomalies (m.w.e. / a)
pdd_params['monthly_dps'] = np.zeros(12)
