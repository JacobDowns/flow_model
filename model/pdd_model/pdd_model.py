import numpy as np
from dolfin import *
from pdd_calculator import PDDCalculator
from pdd_params import pdd_params
import matplotlib.pyplot as plt

"""
Positive degree day model. Computes a surface mass balance function adot. 
"""

class PDDModel(object):

    def __init__(self, model_wrapper, params= {}):

        self.model_wrapper = model_wrapper

        # PDD model parameters
        self.pdd_params = pdd_params
        self.pdd_params.update(params)
        # Surface mass blance function
        self.adot = Function(self.model_wrapper.V_cg)
        # Precipitation for the given year (for plotting) in m/a
        self.precip_func = Function(self.model_wrapper.V_cg)
        # Temperature for the given year (for plotting) in C
        self.temp = Function(self.model_wrapper.V_cg)
        # Object for calculating PDD's
        self.pdd_calc = PDDCalculator(self.pdd_params['pdd_var'])
        

    """
    Update the model SMB function.
    """
    def update(self, params):

        ### Temperature, precipitation, and other parameters to use for the PDD model
        ########################################################################

        self.pdd_params.update(params)

        # Monthly temperature anomalies (C)
        monthly_dts = self.pdd_params['monthly_dts']
        # Monthly precipitation anomalies (m.w.e. / a)
        monthly_dps = self.pdd_params['monthly_dps']
        # Elevation lapse rate (degrees C / km)
        lapse_rate = self.pdd_params['lapse_rate']
        # Ablation rate for ice (m.w.e / PDD) 
        lambda_ice = self.pdd_params['lambda_ice']
        # Ablation rate for snow (m.w.e. / (degree C * day))
        lambda_snow = self.pdd_params['lambda_snow']
        # Ablation rate for ice (m.w.e. / (degree C * day))
        lambda_precip = self.pdd_params['lambda_precip']
        # Superimposed ice fraction
        super_ice_frac = self.pdd_params['super_ice_frac']

        
        ### Compute monthly PDD's and precip.
        ########################################################################
        
        # Get the reference elevation used by climate model
        ref_elevation_vec = self.model_wrapper.input_functions['S_ref'].vector().get_local()
        # Get the modeled elevation
        modeled_elevation_vec = self.model_wrapper.model.S0_c.vector().get_local()
        # Compute the lapse rate correction in C
        lapse_correction = ((ref_elevation_vec - modeled_elevation_vec) / 1000.0) * lapse_rate
        # Total snow that has fallen for the year
        total_snowfall = np.zeros_like(self.model_wrapper.input_functions['S_ref'].vector().get_local())
        # Total number of pdds for the year
        total_pdds = np.zeros_like(self.model_wrapper.input_functions['S_ref'].vector().get_local())
        
        for i in range(12):
            # Compute the delta temp. adjusted / lapse rate corrected temp. for this month
            # Modern temp.  and precip. for a given month are computed as a 30 year modern average from Box
            modern_temp_vec = self.model_wrapper.input_functions['T' + str(i)].vector().get_local()
            temp_vec = modern_temp_vec + monthly_dts[i] + lapse_correction
            # Compute the delta temp. adjusted precip.
            modern_precip_vec = self.model_wrapper.input_functions['P' + str(i)].vector().get_local()
            # Temp. corrected precip. rate in m.w.e./a
            precip_vec = modern_precip_vec*np.e**(lambda_precip*(temp_vec - modern_temp_vec)) + monthly_dps[i]
            # Compute pdd's for this month
            total_pdds += self.pdd_calc.get_pdd(temp_vec)
            # Fraction of precip. that falls as snow
            snowfall_frac = self.pdd_calc.get_acc_frac(temp_vec)    
            # Compute snowfall for the month in m.w.e
            total_snowfall += precip_vec * (1./12.) * snowfall_frac

        # Save total snowfall for plotting
        self.precip_func.vector()[:] = total_snowfall


        ### Compute SMB from total snowfall and pdds
        ########################################################################
        
        # PDD's needed to melt given fraction of the snowpack 
        pdds_to_make_super_ice = (super_ice_frac*total_snowfall) / lambda_snow
        # PDD's that actually go to making superimposed ice
        pdds_super_ice = np.minimum(total_pdds, pdds_to_make_super_ice)
        total_pdds -= pdds_super_ice
        # Amount of superimposed ice in m.w.e
        super_ice = pdds_super_ice * lambda_snow
        # Amount of snow in m.w.e remaining after some has been converted to superimposed ice
        total_snowfall -= super_ice        
        # PDD's needed to melt all the remaining snow
        pdds_to_melt_snow = total_snowfall / lambda_snow
        # PDD's that actually go to melting snow
        pdds_melt_snow = np.minimum(total_pdds, pdds_to_melt_snow)
        total_pdds -= pdds_melt_snow
        # Amount of snow that remains in m.w.e
        total_snowfall -= pdds_melt_snow * lambda_snow
        # PDD's needed to melt the superimposed ice
        pdds_to_melt_super_ice = super_ice / lambda_ice
        # PDD's that actually go to melting superimposed ice
        pdds_melt_super_ice = np.minimum(total_pdds, pdds_to_melt_super_ice)
        total_pdds -= pdds_melt_super_ice
        # The amount of superimposed ice remaining
        super_ice -= pdds_melt_super_ice * lambda_ice
        # Compute the accumulation in m.w.e, consisting of snow and superimpsed ice
        accumulation = total_snowfall + super_ice
        # Compute the amount of ablation in m.w.e (remaining PDD's go to melting glacier ice)
        ablation = total_pdds * lambda_ice
        # Total yearly mass balance in m.i.e. assuming snowpack turns to ice at end of year
        smb = (accumulation - ablation) * (10./9.)

        self.model_wrapper.model.adot.vector()[:] = smb
