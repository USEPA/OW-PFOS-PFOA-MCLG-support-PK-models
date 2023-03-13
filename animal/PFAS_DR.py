# -*- coding: utf-8 -*-
"""
Created on Fri May 14 08:25:30 2021

@author: EPA
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.integrate as sci
import scipy.interpolate as scinterp
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='rpy2') # Ignore 'BayesPBPKTools' namespace warning as it's not needed for importing posterior chains
warnings.filterwarnings("ignore", category=FutureWarning) # Ignore future warnings in rpy2

import rpy2
from rpy2.robjects import r
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri as p2r
from rpy2.robjects import numpy2ri as n2r
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpack
from rpy2.robjects.packages import importr

from rpy2.robjects.conversion import localconverter

import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)

NULL = r("NULL")
data_frame = r("data.frame")
RMCSim = importr('RMCSim')

#n2r.activate()
#p2r.activate()

class PFAS_DR():
    def __init__(self, chem='PFOA', sex='Female', route='oral', strain='CD1', 
                 Qcc=8.68, dose=1, dose_i=0, sample_posterior=False, times=None, ts=0, tf=10, time_units='days',
                 param_path='../data_files', model_path='../pfoa_2compabandersenoral',
                 tm_daily_dose=None, ti_daily_dose=None, t_gest = 22, t_lact = 21, M_m_1 = 0.25, M_m_2 = 0.273,
                 t_m_1 = 1, t_m_2 = 44, n_f = 10, n_i = 10, M_i_1 = 0.0066, M_i_2 = 0.014, 
                 M_i_3 = 0.028, M_i_4 = 0.25, M_i_5 = np.nan, M_i_6 = np.nan,
                 t_i_1 = 3, t_i_2 = 10, t_i_3 = 17, t_i_4 = 75, t_i_5 = np.nan, t_i_6 = np.nan, p_milk=np.nan,
                 r_milk_0 = 0.001, r_milk_1 = 0.003, r_milk_2 = 0.0054, r_milk_3 = 0.0059, r_milk_4 = np.nan, 
                 food_dose = 0, t_m_start = 0,  t_m_end = 1, t_dose = 0.1, t_dose_i=1, t_i_start = 0, t_i_end = 1,
                 plotting=True, dev=False, sex_i='Female', t_step=0.1, t_length=1000, growth_data=pd.DataFrame(), adult_growth=False):
        """
        Any simulation-specific parameters such as time, dose, route of 
        exposure etc. are defined in the __init__. For developmental 
        simulations, t=0 represents concenption. Therefore, t<0 represents
        pre-gestational exposure while 0 < t < t_gest represents gestational
        exposure. Finally, t_gest < t < t_lact represents lactational exposure.
        
        All pk-parameters for the model are defined through the posterior fits 
        from Wambaugh et al. and defined later.

        Parameters
        ----------
        chem : str, optional
            Chemical being simulated (PFOS or PFOA). The default is 'PFOA'.
        sex : str, optional
            Sex of animal (Male of Female). The default is 'Female'.
        route : str, optional
            Route of exposure (IV or oral). The default is 'oral'.
        strain : str, optional
            Strain or species of animal. If mouse, specify 'CD1' or 'C57B6'. 
            If running a rat or monkey simulation, specify strain='rat' or
            strain='monkey', respectively. The default is 'CD1'.
        Qcc : float, optional
            Cardiac output for animal (L/h/kg**0.74). The default is 8.68.
        dose : float, optional
            mg/kg dose given to adult animal. The default is 1.
        dose_i : float, optional
            mg/kg dose given to juvenile animal. The default is 0.
        sample_posterior : boolean, optional
            Run an adult simulation with draws from posterior distribution. 
            The default is False.
        times : ndarray, optional
            Time for simulation output. If None, a default np.linspace is
            generated based on start and stop times. The default is None.
        ts : float, optional
            Simulation start time. The default is 0.
        tf : float, optional
            Simulation end time. The default is 10.
        time_units : float, optional
            Units for time-scale in simulation. The default is 'days'.
        param_path : float, optional
            Path to the directory where RData objects from Wambaugh et al. 
            are located. The default is '../data_files'.
        model_path : float, optional
            Path to compiled RMCSim model. The default is '../pfoa_2compabandersenoral'.
        tm_daily_dose : TYPE, optional
            Dosing schedule for adult dosing. 
            None: Single dose at begining of dosing simulation
            [t_1, t_2, ...]: Daily dosing at each hour specified in list.
            For example, [0] represents once daily dosing at hour 0
            'weekly': Once weekly dosing
            The default is None.
        ti_daily_dose : TYPE, optional
            Dosing schedule for juvenile dosing during post-weaning.
            None: Single dose at begining of dosing simulation
            [t_1, t_2, ...]: Daily dosing at each hour specified in list.
            For example, [0] represents once daily dosing at hour 0
            'weekly': Once weekly dosing
            The default is None.
        t_gest : int, optional
            Length (in days) of gestation. The default is 22.
        t_lact : int, optional
            Length (in days) of lactation. The default is 21.
        M_m_1 : float, optional
            Adult body weight (kg) at t_m_1. The default is 0.25.
        M_m_2 : float, optional
            Adult body weight (kg) at t_m_2. The default is 0.273.
        t_m_1 : float, optional
            Time corresponding to M_m_1. The default is 1.
        t_m_2 : TYPE, optional
            Time corresponding to M_m_1. The default is 44.
        n_f : int, optional
            Number of fetuses in developmental simulation. The default is 10.
        n_i : int, optional
            Number of live infants born requring lactional trasnfer. The default is 10.
        M_i_1 : float, optional
            Mass (kg) of infant at time t_i_1. The default is 0.0066.
        M_i_2 : TYPE, optional
            Mass (kg) of infant at time t_i_2. The default is 0.014.
        M_i_3 : TYPE, optional
            Mass (kg) of infant at time t_i_3. The default is 0.028.
        M_i_4 : TYPE, optional
            Mass (kg) of infant at time t_i_4. The default is 0.25.
        M_i_5 : TYPE, optional
            Mass (kg) of infant at time t_i_5. The default is np.nan.
        M_i_6 : TYPE, optional
            Mass (kg) of infant at time t_i_6. The default is np.nan.
        t_i_1 : TYPE, optional
            Time (days) corresponding to infant mass M_i_1. The default is 3.
        t_i_2 : TYPE, optional
            Time (days) corresponding to infant mass M_i_2. The default is 10.
        t_i_3 : TYPE, optional
            Time (days) corresponding to infant mass M_i_3. The default is 17.
        t_i_4 : TYPE, optional
            Time (days) corresponding to infant mass M_i_4. The default is 75.
        t_i_5 : TYPE, optional
            Time (days) corresponding to infant mass M_i_5. The default is np.nan.
        t_i_6 : TYPE, optional
            Time (days) corresponding to infant mass M_i_6. The default is np.nan.
        p_milk : float, optional
            Per body mass rate of milk ingestion (kg/kg/d). 
            If np.nan, a kg/kg/d ingestion rate is calculated from r_milk_*
            The default is np.nan.
        r_milk_0 : TYPE, optional
            Initial milk consumption rate (kg_milk/day). The default is 0.001.
        r_milk_1 : TYPE, optional
            Week 1 milk consumption rate (kg_milk/day). The default is 0.003.
        r_milk_2 : TYPE, optional
            Week 2 milk consumption rate (kg_milk/day). The default is 0.0054.
        r_milk_3 : TYPE, optional
            Week 3 milk consumption rate (kg_milk/day). The default is 0.0059.
        r_milk_4 : TYPE, optional
            Week 4 milk consumption rate (kg_milk/day). The default is np.nan.
        food_dose : int, optional
            food_dose flag for feed consumption in RMCSim model. The default is 0.
        t_m_start : int, optional
            Start day for adult (maternal) dosing. The default is 0.
        t_m_end : int, optional
            End day for adult (maternal) dosing. The default is 1.
        t_dose : float, optional
            Length of dose for adult (dam). For dietary dose, t_dose = 1 for 
            dose spread out over whole day. The default is 0.1.
        t_dose_i : float, optional
            Length of dose for juvenile. For dietary dose, t_dose = 1 for 
            dose spread out over whole day. The default is 1.
        t_i_start : int, optional
            Start day for juvenile dosing. The default is 0.
        t_i_end : int, optional
            End day for juvenile dosing. The default is 1.
        plotting : bool, optional
            If plots are needed set True. If False (e.g. the model needs to be compiled), no axes generated. The default is True.
        dev : bool, optional
            True for a developmental simulation and False for a strictly 
            adult simulation. The default is False.
        sex_i : str, optional
            Sex of infant/juvenile. The default is 'Female'.
        t_step : float, optional
            Time step for ODE solver. The default is 0.1.
        t_length : int, optional
            Number of times to simulate. The default is 1000.
        growth_data : pandas dataframe, optional
            If there is measured growth, pass this information with a pandas 
            dataframe. The default is an empty pd.DataFrame(), i.e. no growth data.
        adult_growth : bool, optional
            True for adult growth or False for constant animal body weight. 
            If True, a linear interpolation with t_m_2 and M_m_2 is used for 
            the adult simulation
            The default is False.

        Returns
        -------
        None.

        """
        
        # Hard coded attributes
        self.rdata_dict = {'PFOS': 'PFOS-analysis-2016-02-01.RData',
                           'PFOA': 'PFOA-analysis-2015-09-17.RData'
                           }
        
        # Molecular weights from Wambaugh et al., 2013
        self.MW_dict = {'PFOS': 538.22,
                        'PFOA': 414.07}
        
        # Default dev parameters for PFOS/PFOA infant model
        self.dev_params = {'PFOS': {'rat': {'Vcc_i': 0.284, # Infant volume of dsitribution
                                            'P_milk': 0.13, # Partition coefficeint into milk
                                            'r_f_m': 0.833, # fetal:maternal ratio at birth
                                            'half_life': 38.8 # Infant half-life
                                            },
                                    'CD1': {'Vcc_i': 0.268, 
                                            'P_milk': 0.32, # SAME AS MOUSE PFOA
                                            'r_f_m': 0.41,
                                            'half_life': 36.87
                                        }
                                    },
                           'PFOA': {'rat': {'Vcc_i': 0.18,
                                            'P_milk': 0.11,
                                            'r_f_m': 0.42,
                                            'half_life': 2.23
                                            },
                                    'CD1': {'Vcc_i': 0.18,
                                            'P_milk': 0.32, 
                                            'r_f_m': 0.25,
                                            'half_life': 18.65
                                        }
                                    }
                           }
        
        
        

    
        # User-defined attributes for class
        self.chem=chem 
        self.sex=sex
        self.route=route
        self.strain=strain 
        self.Qcc=Qcc
        self.dose=dose
        self.dose_i = dose_i
        self.sample_posterior=sample_posterior
        self.times=times
        self.ts=ts # Start of simulation
        self.tf=tf # End of simulation
        self.time_units = time_units
        self.param_path=param_path
        self.model_path=model_path
        self.model_name = model_path.split('/')[-1]
        self.tm_daily_dose=tm_daily_dose # If None, then only one dose
        self.ti_daily_dose=ti_daily_dose
        self.t_gest = t_gest
        self.t_lact = t_lact
        self.M_m_1 = M_m_1
        self.M_m_2 = M_m_2
        self.t_m_1 = t_m_1
        self.t_m_2 = t_m_2
        self.n_f = n_f
        self.n_i = n_i
        self.M_i_1 = M_i_1
        self.M_i_2 = M_i_2
        self.M_i_3 = M_i_3
        self.M_i_4 = M_i_4
        self.M_i_5 = M_i_5
        self.M_i_6 = M_i_6
        self.t_i_1 = t_i_1
        self.t_i_2 = t_i_2
        self.t_i_3 = t_i_3
        self.t_i_4 = t_i_4
        self.t_i_5 = t_i_5
        self.t_i_6 = t_i_6
        self.r_milk_0 =r_milk_0
        self.r_milk_1 = r_milk_1
        self.r_milk_2 = r_milk_2
        self.r_milk_3 = r_milk_3
        self.r_milk_4 = r_milk_4
        self.p_milk = p_milk # Per body mass rate of milk ingestion (kg/kg/d).
        self.food_dose = food_dose
        self.t_m_start =t_m_start # Start of adult dosing
        self.t_m_end = t_m_end # End of adult dosing
        self.t_dose = t_dose
        self.t_dose_i = t_dose_i
        self.t_i_start = t_i_start
        self.t_i_end=t_i_end
        self.plotting=plotting
        self.dev = dev
        self.sex_i = sex_i
        self.dose_orig = dose
        self.dose_i_orig = dose_i
        self.t_step = t_step #simulation steps
        self.t_length = t_length # Time vector length
        self.t_adult = self.t_gest + self.t_lact # Time when infant model switches to adult (wambaugh)
        self.growth_data = growth_data
        self.growth_data_orig = growth_data.copy()
        self.adult_growth = adult_growth
        
        # Private functions to init
        self._define_R_code()
        
        # Set attr for dev parameters
        if strain != 'monkey':
            for k, v in self.dev_params[chem][strain].items():
                setattr(self, k, v)
        
        if plotting:
            if dev:
                self.fig, self.ax = plt.subplots(2,1)
            else:
                self.fig, self.ax = plt.subplots(1,1)
                self.ax = [self.ax]
    def update_p(self, **kwargs):
        """Update simulation paramaters after defining the PFAS object"""
        
        for k,v in kwargs.items():
            setattr(self, k, v)
            if k == 'dose':
                self.dose_orig = self.dose
            if k == 'dose_i':
                self.dose_i_orig = self.dose_i
            if k == 'model_path':
                self.model_name = v.split('/')[-1]
    
    def _load_posterior(self):
        time_cols = ["Qcc", "ka", "k12","Tmc"] # Params for converting hr --> day
        
        rdata_name = self.rdata_dict[self.chem]
        parm_path = self.param_path
        chem = self.chem
        strain = self.strain
        sex = self.sex
        
        params = p2r.rpy2py(self.pk_run.load_model_params(rdata_name, parm_path, chem, strain, sex))

        params['Qcc'] = self.Qcc
        
        if self.dev and self.model_name == 'pfoa_2compabandersenoral':
            sex_i = self.sex_i
            params_i = p2r.rpy2py(self.pk_run.load_model_params(rdata_name, parm_path, chem, strain, sex_i))
            params_i.columns = [x+'_i' for x in params_i.columns]
            params = pd.concat([params, params_i], axis=1)
            time_cols.extend(["ka_i", "k12_i", "Tmc_i"])
        elif self.dev and self.model_name == 'pfoa_2compabandersenoral_1cmptDev':
            # Dev-specific parameters
            params['half_life'] = self.half_life
            params['Vcc_i'] = self.Vcc_i
            params['P_milk'] = self.P_milk
            params['r_f_m'] = self.r_f_m
            params['n_i'] = self.n_i
            params['food_dose'] = self.food_dose
            
            sex_i = self.sex_i
            params_i = p2r.rpy2py(self.pk_run.load_model_params(rdata_name, parm_path, chem, strain, sex_i))
            params_i['Qcc'] = self.Qcc
            
        
        if self.sample_posterior == False:
            params = params.median(axis=0)
            if self.dev:
                params_i = params_i.median(axis=0)
        

        
        self.original_params = params.copy()
        if self.dev:
            self.original_params_i = params_i.copy()
        if self.time_units == 'days' and self.sample_posterior:
            params.loc[:, time_cols] = params.loc[:, time_cols]*24
            if self.dev:
                params_i.loc[:, time_cols] = params_i.loc[:, time_cols]*24
        elif self.time_units == 'days' and not self.sample_posterior:
            params.loc[time_cols] = params.loc[time_cols]*24
            if self.dev:
                params_i.loc[time_cols] = params_i.loc[time_cols]*24
        self.transformed_params = params.copy()
        if self.dev:
            self.transformed_params_i = params_i.copy()
        
        self.new_params = params
        if self.dev:
            self.new_params_i = params_i
    
    def _define_R_code(self):
        # **************************************************
        # USER: Make sure path to Rtools is appropriate here.
        # **************************************************
        
        R_str = """
library(RMCSim) # If not using RMCSim as library, will need to update
Sys.setenv(PATH = paste("C:/Rtools40/usr/bin", Sys.getenv("PATH"), sep=";"))

compile_pfas_model <- function(model.path) {
    #Sys.setenv(PATH = paste("C:/Rtools40/usr/bin", Sys.getenv("PATH"), sep=";"))
    #compile_model("pfoa_2compabandersenoral")
    compile_model(model.path)
}

load_model_params <- function(rdata_name, parm.path, chem, strain, sex) {
    suppressWarnings(load(paste(parm.path, rdata_name, sep='/')))
    parm.cols <- c("Rv2v1", "ka","Vcc","k12","Tmc", "KT","free","Qfilc", "Vfilc")
    
    if (chem == 'PFOS' & strain == 'rat' & sex == 'Male') {
        mcmc.parms <- PFOS.malerat.sim.mcmc
    } else if (chem == 'PFOS' & strain == 'rat' & sex == 'Female') {
        mcmc.parms <- PFOS.femalerat.sim.mcmc
    } else if (chem == 'PFOS' & strain == 'CD1' & sex == 'Male') {
        mcmc.parms <- PFOS.CD1.male.sim.mcmc
    } else if (chem == 'PFOS' & strain == 'CD1' & sex == 'Female') {
        mcmc.parms <-PFOS.CD1.female.sim.mcmc
    } else if (chem == 'PFOS' & strain == 'monkey') {
        mcmc.parms <-PFOS.monkey.sim.mcmc
    } else if (chem == 'PFOA' & strain == 'CD1' & sex == 'Female') {
        mcmc.parms <-PFOA.CD1.sim.mcmc
    } else if (chem == 'PFOA' & strain == 'rat' & sex == 'Male') {
        mcmc.parms <-PFOA.malerat.sim.mcmc
    } else if (chem == 'PFOA' & strain == 'rat' & sex == 'Female') {
        mcmc.parms <-PFOA.femalerat.sim.mcmc
    } else if (chem == 'PFOA' & strain == 'monkey') {
        mcmc.parms <-PFOA.monkey.sim.mcmc
    } else if (chem == 'PFOA' & strain == 'C57B6') {
        mcmc.parms <-PFOA.C57B6.sim.mcmc
    }
        

    new.params <- data.frame(exp(mcmc.parms[,parm.cols]))
    return(new.params)
    
    }

            """
        self.pk_run = rpack.STAP(R_str, "pk_run")
            
    
    def _create_dosing_df(self):
        if self.route == 'oral':
            adult_sv = 'd_m'
            pup_sv = 'd_i'
        elif self.route == 'IV':
            adult_sv = 'd_m_iv'
            pup_sv = 'd_i_iv'
        if self.times is None:
            sim_days_adult = int(np.ceil(np.abs(self.t_m_start-self.t_m_end)))
            sim_days_pup = int(np.ceil(np.abs(self.t_i_start-self.t_i_end)))
        else:
            sim_days_adult = int(np.max(self.times))
            sim_days_pup = int(np.max(self.times))
    
        # Set up adult dosing
        #--------------------
        if not self.tm_daily_dose:
            # Single dose at t_m_start
            adult_dose_on = pd.DataFrame({'var': [adult_sv], 'time': [self.t_m_start],
                                          'value': [self.dose_model/self.t_dose],
                                          'method': 'replace'})
            
        elif self.tm_daily_dose == 'weekly':
            # Assume once weekly dosing
            tz = np.arange(sim_days_adult) + self.t_m_start
            tz = tz[::7]
            adult_dose_on = pd.DataFrame({'var': [adult_sv]*len(tz),
                                          'time': tz,
                                          'value': [self.dose_model/self.t_dose]*len(tz),
                                          'method': ['replace']*len(tz)})
        else:
            tz = np.arange(sim_days_adult) + self.t_m_start
            adult_dose_on = pd.DataFrame({'var': [adult_sv]*len(tz),
                                          'time': tz,
                                          'value': [self.dose_model/self.t_dose]*len(tz),
                                          'method': ['replace']*len(tz)})
        # Set up pup dosing
        #------------------
        if not self.ti_daily_dose:
            pup_dose_on = pd.DataFrame({'var': [pup_sv], 'time': [self.t_i_start],
                                          'value': [self.dose_i_model/self.t_dose_i],
                                          'method': ['replace']})
        else:
            tz_pup = np.arange(sim_days_pup) + self.t_i_start
            pup_dose_on = pd.DataFrame({'var': [pup_sv]*len(tz_pup),
                                        'time': tz_pup,
                                        'value': [self.dose_i_model/self.t_dose_i]*len(tz_pup),
                                        'method': ['replace']*len(tz_pup)})
        adult_dose_off = adult_dose_on.copy()
        adult_dose_off.loc[:, 'value'] = 0
        adult_dose_off.loc[:, 'time'] = adult_dose_off.loc[:, 'time'] + self.t_dose
        
        pup_dose_off = pup_dose_on.copy()
        pup_dose_off.loc[:, 'value'] = 0
        pup_dose_off.loc[:, 'time'] = pup_dose_off.loc[:, 'time'] + self.t_dose_i
        
        df_dose = pd.concat([adult_dose_on, adult_dose_off, pup_dose_on, pup_dose_off], axis=0)
        df_dose.sort_values(by='time', inplace=True)
        df_dose.reset_index(inplace=True, drop=True)
        return(df_dose)
    
    def _correct_AUC(self, df):
        # Zero out infant concentrations and AUC before conception
        new_df = df.copy()
        AUC_0 = new_df.loc[new_df.time == 0, 'AUC_i'].values[0]
        new_df.loc[new_df.time >= 0, 'AUC_i'] = new_df.loc[new_df.time >= 0, 'AUC_i']-AUC_0
        new_df.loc[new_df.time < 0, ['C_i', 'AUC_i']] = 0
        
        return new_df
    
    def dev_sim(self, times, df_dose):
        """Run a developmental simulation where t < 0 indicates dosing before concenption"""
        r.load_model(self.model_path)
        
        new_params, new_params_i = self.new_params.copy(), self.new_params_i.copy()

        
        times = np.round(times, 8)
        MW = self.MW_dict[self.chem]
        conv = 1000*MW/10**6 # Convert mg --> umol
        
        M_pups = self.n_f * self.M_i_1
        if times[0] <= self.t_gest: # Simulation starts before birth
            # Generate points to interpolate for dam mass up to moment of birth    

            if self.t_m_2 <= self.t_gest: # CASE 1: t_m_2 <= t_gest
                if self.M_m_2 < (self.M_m_1 + M_pups): #If M_m_2 is not large enough to account for the mass of the pups at birth, calculate a third maternal mass that is large enough.
                    M_m_birth = self.M_m_1 + M_pups # Maternal mass just before 
                else:
                    M_m_birth = self.M_m_2
                t_m_vec = [self.t_m_1, self.t_m_2, self.t_gest]
                M_m_vec = [self.M_m_1, self.M_m_2, M_m_birth]
            
            elif self.t_m_1 <= self.t_gest: # CASE 2: t_m_1 <= t_gest
                M_m_birth = self.M_m_2 + M_pups
                t_m_vec = [self.t_m_1, self.t_gest]
                M_m_vec = [self.M_m_1, M_m_birth]
            
            else:
                raise ValueError("Time t_m_1 must be before time t_gest.")
            
            M_mf_in = pd.DataFrame({'times': t_m_vec, 'mass': M_m_vec})
            M_i_in = pd.DataFrame({'times': [0,1], 'mass': [0,0]}) # No infant mass with t < t_gest
            R_milk_in = pd.DataFrame({'times': [0,1], 'rate': [0,0]}) # No milk rate with t < t_gest
    
            M_m_b = M_m_birth - M_pups
            
            # Ratio of concentrations in mother and mother + fetus(es). At conception,
            # this ratio is 1. At birth this ratio can be calculated as shown below.
            # assume a linear transition in the values of this ratio between conception
            # and birth.
            r_m_mf_0 = 1.0
            r_m_mf_b = (M_m_b + M_pups) / (M_m_b + self.r_f_m * M_pups)
            r_m_mf_in = pd.DataFrame({'times': [0, self.t_gest], 'ratio': [r_m_mf_0, r_m_mf_b]})
            
        
            # Initialize the parameters
            parms = r.initParms(p2r.py2rpy(new_params))
            # Initialize the state variables
            Y0 = r.initStates(p2r.py2rpy(parms))
            times1 = times[times<=self.t_gest]
            df_dose1 = df_dose[df_dose.time <= np.max(times1)]
            df_dose1_index = df_dose1.index
            times1 = np.sort(np.unique(np.round(np.append(times1, df_dose1.time.values),8))) # Make sure dosing times are integrated

            
            
            # Convert the dataframes
            with localconverter(ro.default_converter + p2r.converter):
                M_mf_in = ro.conversion.py2rpy(M_mf_in)
                M_i_in = ro.conversion.py2rpy(M_i_in)
                R_milk_in = ro.conversion.py2rpy(R_milk_in)
                r_m_mf_in = ro.conversion.py2rpy(r_m_mf_in)
                df_dose1 = ro.conversion.py2rpy(df_dose1)


            out1 = p2r.rpy2py(data_frame(r.run_model(self.model_name, n2r.numpy2rpy(times1), Y0, parms,
                                    forcing=r.list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in),
                                    fcontrol=r.list(method="linear", rule=2, ties="ordered"),
                                    event_list=r.list(data=df_dose1), #method='lsodes',
                                    rtol=1e-8, atol=1e-8
                                    )))
            
            
            
            self.out1 = out1
            

            if np.max(times) <= self.t_gest:
                if times[0] < 0:
                    out1 = self._correct_AUC(out1)
                return out1
        # Continue with a simulation of the lactation/nursing period if necessary...
        if times[0] <= self.t_gest:
            # Apportion the amount of substance in the mother at the moment of birth to
            # the mother and the infant(s) according to body mass. Use these values to
            # set the initial conditions for the post-birth simulation.
            day_idx = int(1/self.t_step) # Number of index on tail for one day
            M_m_birth = float(out1.iloc[-1]["M_mf"] - M_pups)
            M_i_birth = self.M_i_1
            
            #C_m_birth = float(out1.iloc[-1]["C_m_kg"])
            #C_i_birth = float(out1.iloc[-1]["C_i_kg"]) # Original fetal concentration
            
            #C_i_birth = np.mean(out1.iloc[-day_idx:]["C_i_kg"]) # Average over last day
            C_i_birth = np.mean(out1.iloc[-day_idx:]["C_i"]) # Average over last day
            
            Adeep_birth = float(out1.iloc[-1]["Adeep"])
            Afil_birth = float(out1.iloc[-1]["Afil"])
            Agut_birth = float(out1.iloc[-1]["Agut"])

            #Y0.rx["A_mf"] = C_m_birth * M_m_birth
            #A_i = C_i_birth * M_i_birth
            A_i = C_i_birth  * self.Vcc_i * M_i_birth
            
            #Y0.rx["A_mf"] = float(out1.iloc[-1]["A_mf"]) 
            Y0.rx["A_mf"] = float(np.mean(out1.iloc[-day_idx:]["A_mf"]) - self.n_f*A_i) # Average over final day of gestation for transport
            
            #print('A_mf: %0.3f'%(C_m_birth * M_m_birth))
            Y0.rx["A_i"] = float(A_i)
            Y0.rx["Adeep"] = float(Adeep_birth) #- self.n_f*A_i
            Y0.rx["Afil"] = float(Afil_birth)
            Y0.rx["Agut"] = float(Agut_birth)
            
            Y0.rx["AUC_m"] = float(out1.iloc[-1]["AUC_m"])
            Y0.rx["AUC_i"] = float(out1.iloc[-1]["AUC_i"])
            Y0.rx["d_m"] = float(out1.iloc[-1]["d_m"])
            Y0.rx["d_m_iv"] = float(out1.iloc[-1]["d_m_iv"])
            Y0.rx["d_i"] = float(out1.iloc[-1]["d_i"])
            Y0.rx["T_in"] = float(out1.iloc[-1]["T_in"])
            Y0.rx["T_out"] = float(out1.iloc[-1]["T_out"] + A_i * (self.n_f - self.n_i))
       
        else:
             # Initialize the parameters
            parms = r.initParms(p2r.py2rpy(new_params))
            parms.rx['P_milk'] = self.P_milk
            parms.rx['r_f_m'] = self.r_f_m
            parms.rx['n_i'] = self.n_i
            parms.rx['food_dose'] = self.food_dose
            parms = r.initParms(parms)
            
            # Initialize the state variables
            Y0 = r.initStates(p2r.py2rpy(parms))
        
        r_m_mf_in = pd.DataFrame({'times': [0,1], 'ratio': [0,0]})
        
        # Generate points to interpolate for calculations of dam mass after moment of birth.
        if self.t_m_2 <= self.t_gest:
            if self.M_m_2 < (self.M_m_1 + M_pups):
                M_m_birth = self.M_m_1
            else:
                M_m_birth = self.M_m_2 - M_pups
            t_m_vec = [self.t_gest, self.t_gest+1]
            M_m_vec = [M_m_birth, M_m_birth]
        elif self.t_m_1 <= self.t_gest:
            M_m_birth = self.M_m_2
            t_m_vec = [self.t_gest, self.t_gest+1]
            M_m_vec = [M_m_birth, M_m_birth]
        else:
            raise ValueError("Time t_m_1 must be before time t_gest.")
        
        M_mf_in = pd.DataFrame({'times': t_m_vec, 'mass': M_m_vec})
        
        # Rat infant mass interpolation points.
        if np.isnan(self.M_i_5):
            t_i_vec = self.t_gest + np.array([-1e-6, self.t_i_1, self.t_i_2, self.t_i_3, self.t_i_4])
            M_i_vec = np.array([self.M_i_1, self.M_i_1, self.M_i_2, self.M_i_3, self.M_i_4])
        else:
            t_i_vec = self.t_gest + np.array([-1e-6, self.t_i_1, self.t_i_2, self.t_i_3, self.t_i_4, self.t_i_5, self.t_i_6])
            M_i_vec = np.array([self.M_i_1, self.M_i_1, self.M_i_2, self.M_i_3, self.M_i_4, self.M_i_5, self.M_i_6])
        M_i_in = pd.DataFrame({'times': t_i_vec, 'mass': M_i_vec})
        
        # Rat breast milk ingestion rate interpolation points.
        if np.isnan(self.p_milk):
            t_milk_0 = self.t_gest
            t_milk_1e = t_milk_0 + 7
            t_milk_2e = t_milk_0 + 14
            t_milk_3e = t_milk_0 + 21
            r_milk_0 = self.r_milk_0
            r_milk_1e = r_milk_0 + 7.0 * (self.r_milk_1 - self.r_milk_0) / 3.5
            r_milk_2e = r_milk_1e + 7.0 * (self.r_milk_2 - r_milk_1e) / 3.5
            r_milk_3e = r_milk_2e + 7.0 * (self.r_milk_3 - r_milk_2e) / 3.5
            
            if np.isnan(self.r_milk_4):
                t_milk_4e = t_milk_0 + self.t_lact
                r_milk_4e = 0.0
                t_milk_vec = np.array([t_milk_0, t_milk_1e, t_milk_2e, t_milk_3e, t_milk_4e])
                r_milk_vec = np.array([r_milk_0, r_milk_1e, r_milk_2e, r_milk_3e, r_milk_4e])
            else:
                t_milk_4e = t_milk_0 + 28
                r_milk_4e = r_milk_3e + 7.0 * (self.r_milk_4 - r_milk_3e) / 3.5
                t_milk_5e = t_milk_0 + self.t_lact
                r_milk_5e = 0.0
                t_milk_vec = np.array([t_milk_0, t_milk_1e, t_milk_2e, t_milk_3e, t_milk_4e, t_milk_5e])
                r_milk_vec = np.array([r_milk_0, r_milk_1e, r_milk_2e, r_milk_3e, r_milk_4e, r_milk_5e])
        else: # base milk consumption on infant mass if given a p_milk
            t_wean = self.t_gest + self.t_lact
            M_fn = scinterp.interp1d(t_i_vec, M_i_vec)
            M_wean = M_fn(t_wean)
            t_milk_vec = np.append(t_i_vec, t_wean)
            r_milk_vec = self.p_milk * np.append(M_i_vec, M_wean)
            r_milk_vec = r_milk_vec[t_milk_vec <= t_wean]
            t_milk_vec = t_milk_vec[t_milk_vec <= t_wean]
            t_milk_vec = np.append(t_milk_vec, t_wean)
            r_milk_vec = np.append(r_milk_vec, 0.0)
        
        R_milk_in = pd.DataFrame({'times': t_milk_vec, 'rate': r_milk_vec})
        
        #R_milk_in = pd.DataFrame({'times': t_milk_vec, 'rate': r_milk_vec/conv*1000})
        #display(R_milk_in)
        
        times2 = times[(times >= self.t_gest) & (times < self.t_adult)]
        df_dose2 = df_dose[(df_dose.time >= self.t_gest) & (df_dose.time < self.t_adult)]
        
        #df_dose2 = df_dose.loc[~df_dose.index.isin(df_dose1_index)]
        times2 = np.sort(np.unique(np.round(np.append(times2, df_dose2.time.values), 8))) # Make sure dosing times are integrated

        # Convert the dataframes
        
        
        with localconverter(ro.default_converter + p2r.converter):
            M_mf_in = ro.conversion.py2rpy(M_mf_in)
            M_i_in = ro.conversion.py2rpy(M_i_in)
            R_milk_in = ro.conversion.py2rpy(R_milk_in)
            r_m_mf_in = ro.conversion.py2rpy(r_m_mf_in)
            df_dose2 = ro.conversion.py2rpy(df_dose2)
        out2 = p2r.rpy2py(data_frame(r.run_model(self.model_name, n2r.numpy2rpy(times2), Y0, parms,
                                    forcing=r.list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in),
                                    fcontrol=r.list(method="linear", rule=2, ties="ordered"),
                                    event_list=r.list(data=df_dose2), #method='lsodes', 
                                    rtol=1e-8, atol=1e-8
                                    )))
        self.out2 = out2
        
        
        if times[0] <= self.t_gest:
            out = pd.concat([out1[:-1], out2], axis=0)
            #out = pd.concat([out1, out2], axis=0)
            out.reset_index(drop=True, inplace=True)
        else:
            out = out2
        
        if times[0] < 0:
            out = self._correct_AUC(out)
   
        if np.max(times) < self.t_adult:
            return out
        
        # If times > t_adult, need to switch to wambaugh model and continue growth
        times3 = times[times >= self.t_adult]
        df_dose3 = df_dose[(df_dose.time >= self.t_adult)].copy()
        df_dose3.replace({'d_i': 'd_m'}, inplace=True) # Pup is now the adult
        df_dose3_index = df_dose3.index
        times3 = np.sort(np.unique(np.round(np.append(times3, df_dose3.time.values),8)))

        # Generate points to interpolate for calculations of mass after switch to adult model
        # This uses the data already extracted from the literature
        if self.strain == 'rat' and self.growth_data.empty:
            growth_data = pd.read_csv(self.param_path+'/rat_growth.csv', index_col=0, header=None).T
            growth_data = growth_data[(growth_data.strain == 'sprague-dawley') & (growth_data.sex == self.sex_i)]
            growth_data['age'] = growth_data['age'].astype(float)
            growth_data['BW_mean'] = growth_data['BW_mean'].astype(float)
            
        elif not self.growth_data.empty:
            growth_data = self.growth_data
        
        growth_data['age'] = growth_data['age']+ self.t_gest # Adjust for gestation
        growth_data['BW_mean'] = growth_data['BW_mean']/1000 # Convert to kg
        
        
        if times[0] <= self.t_adult:
            current_age = out.iloc[-1]['time']
            current_weight = out.iloc[-1]['M_i']
        else:
            current_age = self.t_i_4
            current_weight = self.M_i_4

        growth_data = growth_data[(growth_data.age > current_age) & (growth_data.BW_mean > current_weight)]

        t_m_vec = [current_age]
        t_m_vec.extend(list(growth_data.age))
        M_m_vec = [current_weight]
        M_m_vec.extend(list(growth_data.BW_mean))
        
        M_mf_in = pd.DataFrame({'times': t_m_vec, 'mass': M_m_vec})
        M_i_in = pd.DataFrame({'times':[0, 1], 'mass':[0, 0]})
        R_milk_in = pd.DataFrame({'times':[0, 1], 'rate':[0, 0]})
        r_m_mf_in = pd.DataFrame({'times':[0, 1], 'ratio':[1, 1]})
        
        
        # Transfer infant amounts to main model
        # Initialize the parameters
        
        parms = r.initParms(p2r.py2rpy(new_params_i)) # Original
        
        # Initialize the state variables
        Y0 = r.initStates(p2r.py2rpy(parms))
        if times[0] <= self.t_adult:
            A_m_init = float(out.iloc[-1]["A_i"]) # New initial condition is infant amount
        else:
            A_m_init = 0
        Y0.rx["A_mf"] = A_m_init
        Y0.rx['T_in'] = Y0.rx["A_mf"]
        
        with localconverter(ro.default_converter + p2r.converter):
            M_mf_in = ro.conversion.py2rpy(M_mf_in)
            M_i_in = ro.conversion.py2rpy(M_i_in)
            R_milk_in = ro.conversion.py2rpy(R_milk_in)
            r_m_mf_in = ro.conversion.py2rpy(r_m_mf_in)
            df_dose3 = ro.conversion.py2rpy(df_dose3)
        
        out3 = p2r.rpy2py(data_frame(r.run_model(self.model_name, n2r.numpy2rpy(times3), Y0, parms,
                                    forcing=r.list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in),
                                    fcontrol=r.list(method="linear", rule=2, ties="ordered"),
                                    event_list=r.list(data=df_dose3)#, method='lsodes'
                                    )))
        # Zero out the "fetal" concentrations for grwoing pup
        out3['C_i'] = 0
        out3['C_i_kg'] = 0
        out3_mapper = {'C_m': 'C_i',
                       'C_i': 'C_m',
                       'C_m_kg': 'C_i_kg',
                       'C_i_kg': 'C_m_kg',
                       'M_mf': 'M_i',
                       'M_i': 'M_mf'
                       
            }
        out3.rename(columns = out3_mapper, inplace=True)
        self.out3 = out3
        out = pd.concat([out, out3], axis=0)
        out.reset_index(drop=True, inplace=True)
        

        
        return out

        
    def adult_sim(self, times, df_dose):
        """Run an adult simulation where t is never less than 0"""
        r.load_model(self.model_path)
        #assert np.min(times) >= 0, "Times for adult simulation must be greater than zero"
        new_params = self.new_params.copy()
        M_mf = self.M_m_1 # Constant M_mf mass at M_m_1
        if self.adult_growth:
            M_mf_in = pd.DataFrame({'times': [self.t_m_1,self.t_m_2], 'mass':[self.M_m_1, self.M_m_2]})
        else:
            M_mf_in = pd.DataFrame({'times': [0,1], 'mass':[M_mf, M_mf]})
        M_i_in = pd.DataFrame({'times':[0, 1], 'mass':[0, 0]})
        R_milk_in = pd.DataFrame({'times':[0, 1], 'rate':[0, 0]})
        r_m_mf_in = pd.DataFrame({'times':[0, 1], 'ratio':[1, 1]})
    
        # Initialize the parameters
        parms = r.initParms(p2r.py2rpy(new_params))
        Y0 = r.initStates(p2r.py2rpy(parms))
        
        times = np.sort(np.unique(np.append(times, df_dose.time.values)))
        with localconverter(ro.default_converter + p2r.converter):
            M_mf_in = ro.conversion.py2rpy(M_mf_in)
            M_i_in = ro.conversion.py2rpy(M_i_in)
            R_milk_in = ro.conversion.py2rpy(R_milk_in)
            r_m_mf_in = ro.conversion.py2rpy(r_m_mf_in)
            df_dose = ro.conversion.py2rpy(df_dose)

        out = p2r.rpy2py(data_frame(r.run_model(self.model_name, n2r.numpy2rpy(times), Y0, parms, 
                                                forcing=r.list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in),
                                                fcontrol=r.list(method="linear", rule=2, ties="ordered"),
                                                event_list=r.list(data=df_dose),
                                                method='lsodes', atol=1e-10, rtol=1e-10
                                                )))
        return out
    
    def adult_sim_sample(self, times, df_dose):
        """ Run an adult simulation where t is never less than 0
            and using the posterior distributions"""
        r.load_model(self.model_path)
        new_params = self.new_params
        
        M_mf = self.M_m_1 # Constant M_mf mass at M_m_1
        M_mf_in = pd.DataFrame({'times': [0,1], 'mass':[M_mf, M_mf]})
        M_i_in = pd.DataFrame({'times':[0, 1], 'mass':[0, 0]})
        R_milk_in = pd.DataFrame({'times':[0, 1], 'rate':[0, 0]})
        r_m_mf_in = pd.DataFrame({'times':[0, 1], 'ratio':[1, 1]})
        
        times = np.sort(np.unique(np.append(times, df_dose.time.values)))
        with localconverter(ro.default_converter + p2r.converter):
            M_mf_in = ro.conversion.py2rpy(M_mf_in)
            M_i_in = ro.conversion.py2rpy(M_i_in)
            R_milk_in = ro.conversion.py2rpy(R_milk_in)
            r_m_mf_in = ro.conversion.py2rpy(r_m_mf_in)
            df_dose = ro.conversion.py2rpy(df_dose)
            
        all_sims = []
        for ind, row in new_params.iterrows():
            # Initialize parameters for this chain
            parms = r.initParms(p2r.py2rpy(row))
            iter_out = p2r.rpy2py(data_frame(r.run_model(self.model_name, n2r.numpy2rpy(times), parms=parms, 
                                                forcing=r.list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in),
                                                fcontrol=r.list(method="linear", rule=2, ties="ordered"),
                                                event_list=r.list(data=df_dose),
                                                method='lsodes', atol=1e-10, rtol=1e-10
                                                )))
            iter_out.loc[:, 'niter'] = ind # Save the iteration number for later
            all_sims.append(iter_out)
            #if int(ind) >100:
            #    break
        
        out = pd.concat(all_sims)
        return out
    
    
    def run_model(self, dose_units='mg/kg', sens=False, pk_params = pd.Series()):
        """
        Parameters
        ----------
        dose_units : str, optional
            Units on applied dose. The default is 'mg/kg'.
        sens : bool, optional
            True if calling tun_model for a sensitivity analysis. The default is False.
        pk_params : pd.Series, optional
            Manually update selected Wambaugh et al., 2013 posterior pk 
            parameters. The default is an empty pd.Series(), i.e. no change 
            to default pk parameters.

        """
        
        MW = self.MW_dict[self.chem]
        conv = 1000*MW/10**6 # Convert mg --> umol
        #print('current dose', self.dose)
        #print('current infant dose', self.dose_i)
        if dose_units == 'mg/kg': # convert mg dose to umol
            self.dose_model = self.dose/conv
            self.dose_i_model = self.dose_i/conv

             
        if not sens:
            # If this is not a senstivity analysis, load the default params
            self._load_posterior()
            self.new_params.update(pk_params) # Manually change pk parmaters
        df_dose = self._create_dosing_df()
        self.df_dose = df_dose
        if self.times is None:
            times = np.arange(self.ts, self.tf, self.t_step)
            #times = np.linspace(self.ts, self.tf, self.t_length)
        else:
            times = self.times
        
        if self.dev:
            print('Running developmental model')
            print('current adult dose', self.dose)
            print('current infant dose', self.dose_i)
            self.output = self.dev_sim(times, df_dose)
        
        elif self.sample_posterior:
            print('Running adult model with posterior distributions')
            self.output = self.adult_sim_sample(times, df_dose)
        else:
            print('Running adult model for %s with %s mg/kg %s' % (self.strain, self.dose, self.chem))
            self.output = self.adult_sim(times, df_dose)
            self.output = self.output[self.output.time.isin(times)]
        
        

        if dose_units == 'mg/kg':
            # convert umol/L concentration back to mg/L
            do_not_convert = ['time', 'ind', 'k_i', 'r_m_f', 'M_mf', 'M_i', 'R_milk', 'r_m_mf']
            self.output.loc[:, ~self.output.columns.isin(do_not_convert)] = self.output.loc[:, ~self.output.columns.isin(do_not_convert)]*conv
    def plot_time_course(self, label_style='dose', label=None, disp_mean=False, **kwargs):
        """
        Helper function to plot simulations from animal studies

        Parameters
        ----------
        label_style : str, optional
            Label style for each simulation. The default is 'dose'.
        label : TYPE, optional
            If label_style=None, manually define label. The default is None.
        disp_mean : bool, optional
            Plot the average concentration over the last week. The default is False.
        **kwargs : TYPE
            matplotlib kwargs.

        Returns
        -------
        None.

        """
        
        if label_style == 'dose':
            label = '%s mg/kg/day' % self.dose_orig
        
        if self.sample_posterior:
            CI = self.output.groupby('time')['C_m'].quantile([0.05, 0.5, 0.95])
            time = CI.index.unique(level='time')
            CI = CI.unstack().reset_index().drop('time', axis=1)
            self.ax[0].plot(time, CI.loc[:, 0.5], label=label, **kwargs)
            kwargs['linestyle'] = '--'
            self.ax[0].plot(time, CI.loc[:, 0.05], **kwargs)
            self.ax[0].plot(time, CI.loc[:, 0.95], **kwargs)
        
        else:
            if self.dev:
                
                
                self.ax[0].plot(self.output.loc[self.output['time'] <self.t_adult, 'time'], self.output.loc[self.output['time'] <self.t_adult, 'C_m'], **kwargs, label=label)
                
                self.ax[1].plot(self.output['time'], self.output['C_i'], **kwargs)
                self.ax[0].set_xlim(*self.ax[1].get_xlim())
            else:
                self.ax[0].plot(self.output['time'], self.output['C_m'], **kwargs, label=label)
                
        
        if disp_mean:
            if 'alpha' in kwargs.keys():
                kwargs.pop('alpha')
            #mean_conc = [self.output.loc[self.output['time'] > (self.output['time'].max()-self.tf), 'C_m'].mean()]*len(self.output['time'])
            mean_conc = [self.output.loc[self.output['time'] > (self.output['time'].max()-7), 'C_m'].mean()]*len(self.output['time'])
            self.ax[0].plot(self.output['time'], mean_conc, linestyle='--', alpha=1, **kwargs)
        
        
        if label:
            #pass
            self.ax[0].legend(loc='upper left', bbox_to_anchor=(1, 1))
        
        if self.dev:
            self.ax[0].set_ylabel('Dam\n%s conc. [mg/L]'%self.chem)
            self.ax[1].set_xlabel('Time [days]')
            self.ax[1].set_ylabel('Pup\n%s conc. [mg/L]'%self.chem)
        else:
            self.ax[0].set_ylabel('%s serum conc. (mg/L)'%self.chem)
            self.ax[0].set_xlabel('Time (days)')
    def add_features(self, disp_birth=False, disp_lact=False):
        """
        Helper function to add lines for birth and lactational periods

        Parameters
        ----------
        disp_birth : TYPE, optional
            Display birth on ploe. The default is False.
        disp_lact : TYPE, optional
            Display end of lactation on plot. The default is False.

        Returns
        -------
        None.

        """
        ylim0 = self.ax[0].get_ylim()
        if self.dev:
            ylim1 = self.ax[1].get_ylim()
        if disp_birth:
            self.ax[0].vlines(self.t_gest, *ylim0, color='k', linestyle='--', label='birth')
            self.ax[0].set_ylim(ylim0)
            if self.dev:
                self.ax[1].vlines(self.t_gest, *ylim1, color='k', linestyle='--', label='birth')
                self.ax[1].set_ylim(ylim1)
        
        if disp_lact:
            self.ax[0].vlines(self.t_gest+self.t_lact, *ylim0, color='k', linestyle='-.', label='end lact.')
            self.ax[0].set_ylim(ylim0)
            if self.dev:
                self.ax[1].vlines(self.t_gest+self.t_lact, *ylim1, color='k', linestyle='-.', label='end lact.')
                self.ax[1].set_ylim(ylim1)
        self.ax[1].legend(loc='best')
        self.fig.tight_layout()

    def plot_data(self, time, conc, lifestage='adult', label=None, label_style=None, **kwargs):
        """
        Helper function to add measured data to simulation plots

        Parameters
        ----------
        time : list
            List of measured times.
        conc : list
            List of measured concentrations.
        lifestage : TYPE, optional
            Defines which plot the data are added to. The default is 'adult'.
        label : TYPE, optional
            Label for data. The default is None.
        label_style : TYPE, optional
            Use pre-defined label styles of 'sex/dose/route' or 'dose/route' if desired. The default is None.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if lifestage == 'adult':
            idx = 0
        else:
            idx = 1
        
        if label_style == 'sex/dose/route' and not label:
            label = '%s: %s mg/kg, %s' % (self.sex, self.dose_orig, self.route)
        elif label_style == 'dose/route' and not label:
            label = '%s mg/kg, %s' % (self.dose_orig, self.route)
        
        self.ax[idx].scatter(time, conc, label=label, **kwargs)
        if label_style:
            self.ax[idx].legend(loc='best')
    
    def run_sensivity(self, params, dose_metrics, update_params = {}, rel_step=1.0e-2, plot=False):
        """
        

        Parameters
        ----------
        params : list
            list of parameters in PK model to vary for sensitivity analysis.
        dose_metrics : list
            List of dose metrics for comparison.
        update_params : dict, optional
            Model parameters to update from defaults. The default is {}.
        rel_step: float, optional
            Relative step for changing parameter. The default is 1e-2

        Returns
        -------
        None.

        """
        self.run_model(sens=False) # Run model and initialize paramter values
        metrics = self.get_internal_dose_metrics(metrics=dose_metrics, name=0)
        self.metrics = metrics # Baseline metrics for default parameters

        
        # Save the default parameters for the sensitivity analysis
        default_p = self.new_params[params].copy()
        
        sens_list = []
        #print('Initial params')
        #display(self.new_params)
        for i, p in enumerate(params):
            new_p = default_p.copy()
            tmp_p = default_p[p]
            
            # Adjust the parameters value
            if np.abs(tmp_p) >= rel_step/(1. + rel_step):
                new_val = tmp_p*(1 + rel_step)
            elif tmp_p >= 0:
                new_val = rel_step
            else:
                new_val = -rel_step
            self.new_params[p] = new_val
            
            self.update_p(**{'growth_data': self.growth_data_orig.copy()}) # Have to redefine growth data each time because of conversion

            
            self.run_model(sens=True) # Run model with updated parameters
            tmp_metrics = self.get_internal_dose_metrics(metrics=dose_metrics, name=p)
            self.tmp_metrics = tmp_metrics
            tmp_sens = (tmp_metrics - metrics)*default_p[p] / ((new_val - default_p[p])*metrics)
            tmp_sens.name = p
            sens_list.append(tmp_sens)
            self.new_params[default_p.index] = default_p # Reset parameters back to default

        self.sens_list = sens_list
        sens_df = pd.concat(sens_list, axis=1)
        self.sens_df = sens_df
        
        if plot:
            fig, ax = plt.subplots(1,1,figsize=(5,len(params)))
            test = pd.melt(sens_df.reset_index(), id_vars='index', value_vars = sens_df.columns).sort_values(by='value').fillna(0)
            test['index'] = test['index'].str.replace('AU','')
            
            self.test = test
            g = sns.barplot(x='value',y='variable',hue='index', data=test, palette='viridis', ax=ax, hue_order=[x.replace('AU','') for x in dose_metrics])
            legend_labels, _= ax.get_legend_handles_labels()
            g.legend(loc='center left', bbox_to_anchor=(1, 0.5))


    
    def calc_internal_dose(self, dose_metric, days=7):
        """ Available dose metrics:
            
            Adult (non-devlopmental)
            #-----------------------
            Cmax: Maximum adult concentration from non-developmental study
            Css: Average adult concentration over n days of non-developmental study
            Cavg: Average adult concentration over entire non-developmental study
            AUC: Area under the curve for adult in non-developmental study
            
            Dam (developmental)
            #------------------
            Cmax_dam: Maximum dam concentration during gestation
        """
        
        t_lact_end = self.t_gest+self.t_lact-self.t_step
        
        # Adult
        #------
        if dose_metric == 'C7avg':
            metric = self.output.loc[self.output['time'] > (self.output['time'].max()-days), 'C_m'].mean()
        if dose_metric == 'Cfil7avg':
            metric = self.output.loc[self.output['time'] > (self.output['time'].max()-days), 'Cfil'].mean()
        elif dose_metric == 'AUC': # USE WHATS IN THE MODEL
            #AUC = sci.trapz(self.output['C_m'], x=self.output['time'])
            metric = self.output.iloc[-1]['AUC_m']
        elif dose_metric == 'Cavg':
            metric = self.output['C_m'].mean()
        elif dose_metric == 'Cmax':
            metric = self.output['C_m'].max()
        elif dose_metric == 'AUCavg':
            metric = self.output.iloc[-1]['AUC_m']/self.tf
        
        # Dam Development
        #----------------
        elif dose_metric == 'Cmax_pre':
            metric = self.output.loc[(self.output.time < 0), 'C_m'].max()
        elif dose_metric == 'AUCavg_pre':
            metric = self.output.loc[self.output['time']==0-self.t_dose, 'AUC_m'].values[0]/np.abs(self.ts)
        
        elif dose_metric == 'Cmax_dam':
            metric = self.output.loc[(self.output.time >= 0) & (self.output.time <= self.t_gest), 'C_m'].max()
        elif dose_metric == 'Cavg_dam':
            metric = self.output.loc[(self.output.time >= 0) & (self.output.time <= self.t_gest), 'C_m'].mean()
        elif dose_metric == 'AUCavg_dam_gest':
            metric = self.output.loc[self.output['time']==self.t_gest, 'AUC_m'].values[0]/self.t_gest
        elif dose_metric == 'AUCavg_dam_lact':
            AUC_gest = self.output.loc[self.output['time']==self.t_gest, 'AUC_m'].values[0] # AUC at end of gestation
            AUC_lact = self.output.loc[self.output['time']==t_lact_end, 'AUC_m'].values[0]  # AUC at end of lactation
            metric = (AUC_lact-AUC_gest)/self.t_lact # Divide by total days of lactation
        elif dose_metric == 'AUCavg_dam_gest_lact':
            metric = self.output.loc[self.output['time']==t_lact_end, 'AUC_m'].values[0]/t_lact_end
        
        # Pup development
        elif dose_metric == 'Cavg_pup_gest':
            metric = self.output.loc[(self.output.time >= 0) & (self.output.time < self.t_gest), 'C_i'].mean()
        elif dose_metric == 'Cmax_pup_gest':
            metric = self.output.loc[(self.output.time >= 0) & (self.output.time < self.t_gest), 'C_i'].max()
        elif dose_metric == 'Cavg_pup_lact':
            metric = self.output.loc[(self.output.time >= self.t_gest) & (self.output.time <= t_lact_end), 'C_i'].mean()
        elif dose_metric == 'Cmax_pup_lact':
            metric = self.output.loc[(self.output.time >= self.t_gest) & (self.output.time <= t_lact_end), 'C_i'].max()
        elif dose_metric == 'Cavg_pup_gest_lact':
            metric = self.output.loc[(self.output.time >= 0) & (self.output.time <= t_lact_end), 'C_i'].mean()
        elif dose_metric == 'Cmax_pup_wean':
            metric = self.output.loc[(self.output.time >= t_lact_end), 'C_i'].max()
        elif dose_metric == 'Css_pup_wean':
            metric = self.output.loc[self.output['time'] > (self.output['time'].max()-days), 'C_i'].mean()
        elif dose_metric == 'AUCavg_pup_gest':
            metric = self.output.loc[self.output['time']==self.t_gest, 'AUC_i'].values[0]/self.t_gest
        elif dose_metric == 'AUCavg_pup_lact':
            AUC_gest = self.output.loc[self.output['time']==self.t_gest, 'AUC_i'].values[0] # AUC at end of gestation
            AUC_lact = self.output.loc[self.output['time']==t_lact_end, 'AUC_i'].values[0]  # AUC at end of lactation
            metric = (AUC_lact-AUC_gest)/self.t_lact # Divide by total days of lactation
        elif dose_metric == 'AUCavg_pup_gest_lact':
            metric = self.output.loc[self.output['time']==t_lact_end, 'AUC_i'].values[0]/t_lact_end
        elif dose_metric == 'AUCavg_pup_diet':
            t_diet = self.tf - (self.t_gest+self.t_lact) # total time for dietary exposure
            metric = self.output.iloc[-1]['AUC_i']/t_diet
        elif dose_metric == 'AUCavg_pup_total':
            # Need to add AUC from gestational/lactaional to dietary AUC because of model switch
            AUC_gest_lact = self.output.loc[self.output['time']==t_lact_end, 'AUC_i'].values[0]
            AUC_diet = self.output.iloc[-1]['AUC_i']
            metric = (AUC_gest_lact+AUC_diet)/self.tf
        return metric

    
    def get_internal_dose_metrics(self, metrics=['mean_conc', 'AUC', 'Cmax'],name=None):
        hero_metrics = {}
        for metric in metrics:
            hero_metrics[metric] = self.calc_internal_dose(metric)
        df = pd.Series(hero_metrics)
        if name is not None:
            df.name = name
        else:
            df.name = '%s mg/kg' % self.dose_orig
        return df
    
    def plot_internal_dose_response(self, response, yerr=None):
        fig, ax  = plt.subplots(1,1)
        if yerr:
            ax.errorbar(self.internal_dose, response, yerr=yerr, fmt='o', linestyle='')
        else:
            ax.plot(self.internal_dose, response, 'o', linestyle='')
        return fig, ax
    def calc_pk_params_sample(self):
        all_iters = self.output.niter.unique()
        halft_dist, Vd_beta_dist, CLC_dist = [], [], []
        for niter in all_iters:
            beta, Vd_beta, CLC = self.calc_pk_params(phase='beta', disp=False, plot_phase=False, niter=niter)
            halft_dist.append(beta)
            Vd_beta_dist.append(Vd_beta)
            CLC_dist.append(CLC)
        low_halft, med_halft, high_halft = np.quantile(halft_dist, [0.05, 0.5, 0.95])
        low_Vd_beta, med_Vd_beta, high_Vd_beta = np.quantile(Vd_beta_dist, [0.05, 0.5, 0.95])
        low_CLC, med_CLC, high_CLC = np.quantile(CLC_dist, [0.05, 0.5, 0.95])
        print("half_t (days): %0.5f (%0.5f - %0.5f)" % (med_halft, low_halft, high_halft))
        print("Vd_beta (L/kg): %0.5f (%0.5f - %0.5f)" % (med_Vd_beta, low_Vd_beta, high_Vd_beta))
        print("CLC (L/kg/day): %0.5f (%0.5f - %0.5f)" % (med_CLC, low_CLC, high_CLC))        
    
    def calc_pk_params(self, phase='beta', disp=True, plot_phase=True, niter=None):
        #assert self.output, "Simulation must be run first (run_model)"
        
        if niter:
            data = self.output[self.output.niter == niter].copy()
        else:
            data = self.output.copy()
        n_time, _ = data.shape
        idx = int(n_time*0.15)
        if phase == 'beta':
            coeffs = np.polyfit(data.time[-idx:], np.log(data.C_m[-idx:]), deg=1)# Slope is first term
        elif phase == 'alpha':
            coeffs = np.polyfit(data.time[1:idx], np.log(data.C_m[1:idx]), deg=1)# Slope is first term
        ke = -coeffs[0]
        
        if plot_phase:
            self.ax[0].plot(data.time, np.exp(np.polyval(coeffs, data.time)), 'r')
        
        beta = np.log(2)/ke
        dose_mg = self.dose_orig*self.M_m_1 # mg dose applied
        AUC = sci.trapezoid(data.C_m, x=data.time)
        AUC_extra = data.C_m.iloc[-1]/ke
        AUC_inf = AUC+AUC_extra
        # CLC = D/AUC_inf
        CLC = self.dose_orig/(AUC_inf)
        Vd_beta = CLC/ke
        if disp:
            print("half_t: %0.3f, Vd_%s: %0.5f, CLC: %0.4f, ke: %0.5f"%(beta, phase, Vd_beta, CLC, ke))
        return beta, Vd_beta, CLC
        
        
    