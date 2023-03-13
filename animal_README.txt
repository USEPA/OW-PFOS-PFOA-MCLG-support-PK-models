============
ANIMAL FILES
============

===========
Data Files:
===========
These are auxillary files for model parameters and growth curves that are called during the anlaysis.

animal/data_files/NTP_growth.xlsx								BW Growth file for NTP study
animal/data_files/rat_growth.csv								BW Growth file for generic rat
animal/data_files/PFOA-analysis-2015-09-17.RData						Posterior distributions for PFOA from Wambaugh et al., 2013
animal/data_files/PFOS-analysis-2016-02-01.RData						Posterior distributions for PFOS from Wambaugh et al., 2013
animal/HECD 4-16 PFOA PFOS PBPK Data Needs - Update with Butenhoff Study_09.28.2021.xlsx	List of studies to loop through for dose-response simulations

==========
Notebooks:
==========
For the animal model, each analysis is conducted using a jupyter notebook to lay out which objects are created for a given simulation and demonstrate which methods of PFAS_DR.py are used.

dose-response notebooks
-----------------------
These notebooks generate the internal dose metrics for the studies of interest

animal/dose_response/calc_internal_metrics.ipynb		Dose response for all adult studies across all species
animal/dose_response/rat_dev_DR.ipynb				Dose response simulations for all RAT developmental studies
animal/dose_response/mouse_dev_DR.ipynb				Dose response simulations for all MOUSE developmental studies

model-validation notebooks
--------------------------
These notebooks compare the model in vivo datasets for validation

animal/model_compare/wambaugh_2013_reconstruction_update.ipynb	Compare ported model predictions to those published in Wambaugh et al., 2013
animal/model_compare/OutOfSample_compare.ipynb					Compare ported model to out-of-sample adult PK data in rats
animal/model_compare/sensitivity_analysis.ipynb					Uses sensitivity analysis method for adult and developmental sensitivity analysis
animal/model_compare/model_test_summary.ipynb					Unity line plots for comparing predicted concentrations against all in vivo data (Wanbuagh training data, test data, and developmental)

=========
PK model:
=========
The actual files that run the model.

animal/PFAS_DR.py									This is the class that runs everything
animal/pfoa_2compabandersenoral_1cmptDev.model		This model is used for both adult and dev PK animal studies.

==============
Model outputs:
==============
animal/dose_response/internal_dose_metrics_nonDev_FINAL.xlsx					Internal dose metrics for acute/chronic/subchronic endpoints
animal/dose_response/internal_dose_metrics_RAT_DevRepro_PFOA_FINAL.xlsx			PFOA reproductive and developmental dose-response in RAT
animal/dose_response/internal_dose_metrics_RAT_DevRepro_PFOS_FINAL.xlsx			PFOS reproductive and developmental dose-response in RAT
animal/dose_response/internal_dose_metrics_MOUSE_DevRepro_PFOA_FINAL.xlsx		PFOA reproductive and developmental dose-response in MOUSE
animal/dose_response/internal_dose_metrics_MOUSE_DevRepro_PFOS_FINAL.xlsx		PFOS reproductive and developmental dose-response in MOUSE

