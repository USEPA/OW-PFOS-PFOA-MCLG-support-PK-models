##------------------------------------------------------------------------------
## example.R
##
## Compiles model and loads functions, example calls for forward dosimetry,
## forward dosimetry with multiple doses, reverse dosimetry, dose metric calculation
##
##------------------------------------------------------------------------------

## Load required functions and libraries.
source("compile_pfas_tk.R")
source("sim_functions.R")
source("human_functions.R")

## Reverse dosimetry to determine human intake to acheive X mg/L plasma concentration
## This is for a child, same maternal and child intake (after weaning)
target_conc <- 1 # mg/L
human_dose <- child_rev(target_conc,age=3,sex="Female",chem="PFOA") # in mg/kg/d

## This is the same developmental scenario, but with reverse dosimetry based on the maternal plasma level. (Age is still for the child)
human_dose <- child_rev(target_conc,age=3,sex="Female",chem="PFOA",maternal=T) # in mg/kg/d

## This if for an adult starting with no PFAS body burden
target_conc <- 1 # mg/L
human_dose <- adult_rev(target_conc,age=25,sex="Female",chem="PFOA") # in mg/kg/d

## Function to perform forward dosimetry (internal plasma concentration from intake)
## Child function
human_dose <- 1 # mg/kg/d
human_conc <- human_child_conc(child_age=3,sex="Female",chem="PFOA",dose=human_dose) # in mg/L

## Developmental scenario, maternal plasma concentration
human_dose <- 1 # mg/kg/d
human_conc <- human_child_conc(child_age=3,sex="Female",chem="PFOA",dose=human_dose,maternal=T) # in mg/L

## Adult function
human_dose <- 1 # mg/kg/d
human_conc <- human_adult_conc(age=25,sex="Female",chem="PFOA",dose=human_dose) # in mg/L

## Adult function with multiple dose levels
human_dose <- 1 # mg/kg/d
dose_change  <-  data.frame(d_m=0.5,time=10)# dose in mg/kg/d, time in years from start of simulation)
human_conc <- human_adult_conc(age=25,sex="Female",chem="PFOA",dose=human_dose,dose_change=dose_change) # in mg/L

human_dose <- 0.2 # mg/kg/d
dose_change  <-  data.frame(d_m=c(0.1,0.05),time=c(18,22))## dose in mg/kg/d, time in years from start of sim.)
human_conc <- human_adult_conc(age=25,sex="Female",chem="PFOA",dose=human_dose,dose_change=dose_change) # in mg/L

## Function to generate time course data
## Adult function
human_data <- human_adult(age=50,sex="Female",chem="PFOA",dose=1)
colnames(human_data)# Variables available
plot(human_data[,"years"],human_data[,"C_m"],type='l',
     xlab="Age (yr)",ylab="Plasma Concentration (mg/L)")

## Adult function with changing dose
human_data_change <- human_adult(age=50,sex="Female",chem="PFOA",dose=1,dose_change=data.frame(d_m=0.5,time=30))
plot(human_data[,"years"],human_data[,"C_m"],type='l',
     xlab="Age (yr)",ylab="Plasma Concentration (mg/L)")
lines(human_data_change[,"years"],human_data_change[,"C_m"],lty=2)

## Multiple dose changes
human_data <- human_adult(age=50,sex="Female",chem="PFOA",dose=1,dose_change=data.frame(d_m=c(0.5,0.25),time=c(20,30)))
plot(human_data[,"years"],human_data[,"C_m"],type='l',
     xlab="Age (yr)",ylab="Plasma Concentration (mg/L)")

## Child function
human_data <- human_child(child_age=5,sex="Female",chem="PFOA")
colnames(human_data)# Variables available
plot(human_data[,"years"],human_data[,"C_i"],type='l',col=2,lty=2,
     xlab="Child Age (yr)",ylab="Plasma Concentration (mg/L)")
lines(human_data[,"years"],human_data[,"C_m"])
legend(x="topleft",legend=c("Maternal","Filial"),lty=c(1,2),col=c(1,2))

## Exposure based on drinking water consumption:
human_data <- human_child(child_age=5,sex="Female",chem="PFOA",DW=T)
plot(human_data[,"years"],human_data[,"C_i"],type='l',col=2,lty=2,
     xlab="Child Age (yr)",ylab="Plasma Concentration (mg/L)")
lines(human_data[,"years"],human_data[,"C_m"])
legend(x="topleft",legend=c("Maternal","Filial"),lty=c(1,2),col=c(1,2))

## Fillial dose metrics:
## Peak concentration during gestation/nursing
calc_internal_dose(human_data[,"time"],"Cmax_pup",C_i=human_data[,"C_i"],AUC_i=human_data[,"AUC_i"],
                   t_conc=(24.25*365),t_birth=(25*365),t_wean=(26*365))
## Average concentration during gestation/nursing
calc_internal_dose(human_data[,"time"],"Cavg_pup_gest_lact",C_i=human_data[,"C_i"],AUC_i=human_data[,"AUC_i"],
                   t_conc=(24.25*365),t_birth=(25*365),t_wean=(26*365))
## Average concentration during gestation
calc_internal_dose(human_data[,"time"],"Cavg_pup_gest",C_i=human_data[,"C_i"],AUC_i=human_data[,"AUC_i"],
                   t_conc=(24.25*365),t_birth=(25*365),t_wean=(26*365))
## Average concentration during nursing
calc_internal_dose(human_data[,"time"],"Cavg_pup_lact",C_i=human_data[,"C_i"],AUC_i=human_data[,"AUC_i"],
                   t_conc=(24.25*365),t_birth=(25*365),t_wean=(26*365))

