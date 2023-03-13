==========================
Files in the human folder:
==========================
pfas_tk.model –  Model file that contains the core ODEs of the model, along with various state, output, input, and parameter definitions. Parameters values are changed by R scripts that call the model. The model is written in MCSim language, which is compiled to C and called from the R scripts. The installation of MCSim is required to compile the model. The R library deSolve is required to run the model. The R library xlsx is also loaded to read a spreadsheet of PODs.

sim_functions.R – Original file housing main functions to run the model.

human_functions.R – new code for human model, new parameters are in p_PFOA and p_PFOS, houses all the main functions to run the model, though for developmental scenarios the root function remains in sim_functions. For the adult simulations, all the code to run the model is in this file.

example.R – Compiles model and loads functions, example calls for forward dosimetry, forward dosimetry with multiple doses, reverse dosimetry, dose metric calculation

HED_calc.R – File to calculate human equivalent doses for PODs for the
Proposed Approaches to the Derivation of a Draft Maximum Contaminant
Level Goal for PFOA and PFOS submitted for SAB review.
HED_calc_update.R is an updated version of this file.

HED_calc_update.R – Updated HED calculation file for new PODs for the
Proposed Maximum Contaminant Level Goal for PFOA/PFOS in Drinking Water.
Reads PODs and study information from "For
PODHED_BMD Modeling Summary.xlsx". Outputs HED (in mg/kg/d) in
"HEDs.xlsx". The "xlsx" package, which is used to read and write excel
files, may be incompatible with a new version of RStudio. The
incompatibility does not occur for the basic RGui .

validation.R – File to compare output of human developmental model to selected datasets. For validation, PFOA/S exposure levels were set to match maternal serum levels at delivery. Maternal levels in the Mogensen 2015 study were back-calculated from the levels they reported in children at delivery, which they estimated from maternal levels. Model predictions are compared to measured values at 6 and 19 months (Fromme 2010), 36 mo. (Granum 2013), and at 11, 19, and 60 mo. (Mogensen 2015).

LSA.R – File to perform local sensitivity analysis on selected parameters

verner_rep.R – Code to replicate results of Verner model, which were obtained from running the original acslX model and saving the output in csv form.

human_functions_vbw.R – Functions to replicate Verner model. Uses different parameters, different milk intake table, different r_m_mf table, and different BW curve. Human_child_v is similar to human_child in human_functions.R and calls human_mother_child_v which is similar to human_mother_child in sim_functions.R

Parameters:
===========

Constants:
--------- 
Vd
	170 ml/kg for PFOA, 230 ml/kg for PFOS from (Thompson 2010 2919278), in Methods paragraph end of page 391 and beginning of 392.
Cl
	Calculated from the t1/2 values and the volume of distribution using:
	Cl = Vd * ln(2) / t1/2 ,
	See Cl calc.xlsx. Clearance values in ml/kg/d are divided by 1000 to convert to l/kg/d when clearance is defined in the code. T1/2 values are from 2.7 (Li 2017, Technical report: link, 1st sentence discussion (2.67 for Ronneby), 3rd paragraph discussion (2.72 for C8 population)) for PFOA and 3.4 for PFOS (Li, 2018, 4238434; Table 4).
r_f_m
	See table in main document.
P_milk
	See table in main document.

Age-dependent parameters (input from csv files and implemented as model ‘inputs’ with linear interpolation):

Male BW curve [human_male_mass_0_to_20.csv and human_male_mass_20_to_75.csv]

Quote from Kapraun, DF et al. 2022, Tox. Sci. 189(2), p. 155-174:
    * it is zero before the instant of birth;
    * it is piecewise linear and continuous for ages greater than zero;
    * it reproduces values from Figure 2 (sampled at 3-month intervals starting with age 0 months and ending with age 36 months) of {Kuczmarski, 2002, 3490881@@author-year}, which depicts human female body mass versus age for ages 0 to 36 months;
    * it reproduces values from Figure 10 (sampled a 1-year intervals starting with age 4 years and ending with age 20 years) of {Kuczmarski, 2002, 3490881@@author-year}, which depicts human female body mass versus age for ages 2 to 20 years; and
    * it has values of 67.9, 70.2, 72.7, 73.6, 73.9 and 69.0 kg at ages 25, 35, 45, 55, 65, and 75 years in agreement with Table 8-5 of Chapter 8 of {U.S. EPA, 2011, 786546@@author-year}.
Data was extracted from the two above sources for male children and adults.

Milk consumption
----------------
95th percentile for milk intake (in L/kg*day) is in r_milk_bw95.csv. These values are from Paul’s analysis in r_milk.xlsx {S1:S5}. The linear interpolation of this set of points is design to match the average values from the EPA exposure factors handbook table 15-1, upper percentile estimate {Q2:Q5}. For validation, a mean milk intake function was developed. The result is in r_milk_bw.csv. The raw data from table 15-1 is in cells {P2:P5}. The midpoints for interpolation are in {R1:R5}. The value in {R1} was chosen (not calculated) to result in a smooth function (see graph). Several modifications of these files were generated for specific validation scenarios, which are covered below. r_milk_bw_95_lsa.csv represents a 1% increase in milk intake compared to r_milk_bw_95.csv to examine the effect of a small increase in milk consumption of the serum concentration in children.

Drinking water consumption
--------------------------
95th percentile for drinking water consumption (in ml/kg*d) is in water_95.csv. Mean drinking water consumption is in water.csv. These files are loaded for individuals who are not breastfeeding or pregnant during the simulations, i.e., adults in non-developmental scenarios or the child in a formula-fed scenario. The source for this data was Goeden et al. (2018) model (MDH), which is sourced from the data in the EPA exposure factors handbook Table 3-1 (2011 version, not the 2019 update). These values are in the sheet ModelParameters cells {C16:C27} and {E16:E27} in the May 2018 MDH model spreadsheet. Ages at the ‘turning-points’ for the interpolation were confirmed to match (within the error of a single time-step) with the turning-points in column E. 

water_95_bf.csv and water_bf.csv are the 95th percentile and mean water intakes for the breastfed infant. These are modifications of the respective formula-fed files discussed above and differ in that they initially start with 0 water consumption, then at weaning (1 yr) water consumption is brought up to the value in the formula-fed scenario over a time span of 0.000001 yr (31.5 s).

water_95_lact.csv and water_lact.csv are the 95th percentile and mean water intakes for the pregnant woman. The 95th percentile lactation value matches the number from the MDH model in cell {O6} of the breastfed scenario sheet. The value for the mean is from the EPA exposure factors hand book (2019 update) Table 3-3, rounded to 23 from 22.9 ml/kg*d.

Several modifications of these files were generated for specific validation scenarios, which are covered below.
---------------------------------------------------------------------------------------------------------------
Volume of distribution adjustment
	vd_child_con.csv is used in the scripts and reflects a constant volume of distribution, i.e. Vd is multiplied by 1 at all ages. Vd_child.csv was used to test the effect of greater Vd in children and was sourced from the MDH model. In that model the linear interpolation is hardcoded, and the values at the ‘turning-points’ were taken from the values in columns H (Vd multiplier) and columns A (time in years) of the 2018 MDH model spreadsheet. Functions in human_functions.R and sim_functions.R need to be modified to change the Vd multiplier, it is not passed as a function input.
Validation Scenarios:
	Separate files were generated for validation scenario to account for different lengths of breastfeeding. The three cohorts examined are documented in Fromme 2010, Mogensen 2015, and Granum 2013. The Granum 2013 paper is on a subset of the cohort described in Magnus 2006.  
	In the Fromme 2010 study, breastfeeding information is only available 6 months after birth. At this point 37 of 50 participants with breast feeding status were exclusively breastfed, 6 predominantly breastfed, 6 partially breastfed, and 1 received no breast milk (See 1st paragraph of page 3, just below Figure 1). In the simulation, we treat this as exclusive breastfeeding from birth to 6 months of age. This matches the decision made by Verner, 2015 (see Table 1, “German study”). Data files with the ‘fromme’ suffix account for weaning at 6 months. The calculation for water consumption at 6 months is in water_weaning_calc.xlsx sheet water_mean cell C4.
	The Granum 2013 study had an average breast-feeding duration of 12.8 yr (Table 1). I believe Verner 2013 incorrectly read the value from the following row. Verner cites 13.2 months for breastfeeding length in their Table 1 (“Norwegian Study”), whereas this is the average start of day-care in Granum 2013. Because the model only has breastfeeding parameters up to 1 year (and generation of breastmilk intake for older infants would require the implementation of new data), the simulation of this study used the default scenario (with mean breastmilk consumption) with 1 year of breastfeeding. Having a slightly shorter breastfeeding duration also compensates for the lack of a description of the weaning process, where breastmilk intake is likely decreasing in children breastfeeding after 1 year.
	In the Mogensen 2015 study, the median length of exclusive breastfeeding was 4.5 months, and the median length of partial breastfeeding was 4.0 months (Table 1, survey results at age 11 mo. and older). Two scenarios were developed for this study, milk and water intake values with weaning at 4.5 months and at 8.5 months. Files with a suffix of “mogensen” represent a 8.5 month breastfeeding period, and files with a suffix of “mogensen45” represent a 4.5 month breastfeeding period. The default scenario in the model has milk intake for 1 year, followed by exposure to other PFAS sources at weaning. This timespan is more typical of total (exclusive and partial) breastfeeding, as opposed to exclusive breastfeeding which typically lasts up to around 6 mo. of age. Because of this, 8.5 months was chosen as the breastfeeding duration for the simulation of this study.



