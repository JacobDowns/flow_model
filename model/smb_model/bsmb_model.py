import numpy as np
from dolfin import *
from model.smb_model.pdd_calculator import PDDCalculator
from model.smb_model.smb_params import smb_params
import matplotlib.pyplot as plt

"""
Positive degree day model. Computes a surface mass balance function adot. 
Additionally, adds submarine melt on the floating portion of the ice. Hence

"""

class SMBModel(object):

    def __init__(self, model_wrapper, params= {}):

        self.model_wrapper = model_wrapper

        # PDD model parameters
        self.smb_params = smb_params
        self.smb_params.update(params)
        # Surface mass blance function
        self.adot = Function(self.model_wrapper.V_cg)
        # Solid precip. for the given year (for plotting) in m/a
        self.snowfall_func = Function(self.model_wrapper.V_cg)
        # Liquid precip. for the given year (for plotting) in m/a
        self.rainfall_func = Function(self.model_wrapper.V_cg)
        # Temperature for the given year (for plotting) in C
        self.temp = Function(self.model_wrapper.V_cg)
        # Object for calculating PDD's
        self.pdd_calc = PDDCalculator(self.smb_params['pdd_var'])
        # Fields that need to be loaded
        self.fields = ['P' + str(i) for i in range(12)]
        self.fields += ['T' + str(i) for i in range(12)]
        self.fields.append('smb')
        self.fields.append('submarine_melt_scale')
        # Load model fields
        model_wrapper.load_fields(self.smb_params['fields'], self.fields)


        
    """
    Update the model SMB function.
    """
    def update(self, params):

        ### Temperature, precipitation, and other parameters to use for
        ### the PDD model
        ##############################################################

        
        self.smb_params.update(params)
        #self.model_wrapper.update_interp_fields(self.fields, float(self.model_wrapper.model.L0))

        # Monthly temperature anomalies (C)
        monthly_dts = self.smb_params['monthly_dts']
        # Monthly precipitation anomalies (m.w.e. / a)
        monthly_dps = self.smb_params['monthly_dps']
        # Elevation lapse rate (degrees C / km)
        lapse_rate = self.smb_params['lapse_rate']
        # Ablation rate for ice (m.w.e / PDD) 
        lambda_ice = self.smb_params['lambda_ice']
        # Ablation rate for snow (m.w.e. / (degree C * day))
        lambda_snow = self.smb_params['lambda_snow']
        # Ablation rate for ice (m.w.e. / (degree C * day))
        lambda_precip = self.smb_params['lambda_precip']
        # Superimposed ice fraction
        super_ice_frac = self.smb_params['super_ice_frac']
        # Submarine melt rate (m.i.e. / a)
        submarine_melt_rate = self.smb_params['submarine_melt_rate']
        submarine_melt_scale = self.model_wrapper.input_functions['submarine_melt_scale'].vector().get_local()[0]
        submarine_melt_rate *= submarine_melt_scale
        #print(submarine_melt_rate)

        
        ### Compute monthly PDD's and precip.
        ##############################################################
        
        # Get the reference elevation used by climate model
        ref_elevation_vec = self.model_wrapper.input_functions['S_ref'].vector().get_local()
        # Get the modeled elevation
        modeled_elevation_vec = self.model_wrapper.ice_model.S0_c.vector().get_local()        
        # Compute the lapse rate correction in C
        lapse_correction = ((ref_elevation_vec - modeled_elevation_vec) / 1000.0) * lapse_rate
        # Total snow that has fallen for the year
        total_snowfall = np.zeros_like(self.model_wrapper.input_functions['S_ref'].vector().get_local())
         # Total rain that has fallen for the year
        total_rainfall = np.zeros_like(self.model_wrapper.input_functions['S_ref'].vector().get_local())
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
            # Compute rainfall for the month in m.w.e
            total_rainfall += precip_vec * (1./12.) * (1. - snowfall_frac)


        # Save total snowfall / rainfall for plotting
        self.snowfall_func.vector()[:] = total_snowfall
        self.rainfall_func.vector()[:] = total_rainfall

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

        ### Compute SMB from total snowfall and pdds
        ########################################################################
        B = self.model_wrapper.ice_model.B.vector().get_local()
        Bhat = project(self.model_wrapper.ice_model.Bhat).vector().get_local()
        submarine_melt_rate = ((Bhat - B) > 25.)*submarine_melt_rate

        self.model_wrapper.ice_model.adot.vector()[:] = smb + submarine_melt_rate
