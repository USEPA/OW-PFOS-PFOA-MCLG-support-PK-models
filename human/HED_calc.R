##------------------------------------------------------------------------------
## HED_calc.R
##
## File to calculate human equivalent doses for PODs for the SAB draft.
## HED_calc_update.R is an updated version of this file.
##
##------------------------------------------------------------------------------

source("compile_pfas_tk.R")# Only necessary first time or after .model change
source("sim_functions.R")
source("human_functions.R")
require("xlsx")

## Read spreadsheet that contains point of departure (POD) internal dose metrics
## derived from animal studies and store PODs and related information in a data
## frame called "animal_pods".
animal_pods <- read.xlsx("./PODHED_BMD Modeling Summary-211105.xlsx",sheetIndex=1,colIndex=2:39,rowIndex=5:65)# animal pods in mg/L (based on TZ output)

## Define a mapping between dose metric codes and descriptions of the dose
## metrics.
dm_list <- list("C7avg","AUCavg_pup_gest_lact","AUCavg_pup_total","AUCavg_pup_gest","AUCavg_pup_lact",
                "AUCavg_dam_gest","AUC","Cavg_pup_gest")
names(dm_list) <- c("Average concentration over final week of study (C7,avg)",
                    "AUC normalized per day during gestation/lactation (AUCavg,pup,gest,lact)",
                    "AUC normalized per day over entire study(AUCavg,pup,total)",
                    "AUC normalized per day during gestation (AUCavg,pup,gest)",
                    "AUC normalized per day during lactation (AUCavg,pup,lact)",
                    "AUC normalized per day during gestation (AUCavg,dam,gest)",
                    "AUC for duration of study (AUC)",
                    "Average pup concentration during gestation (Cpup,avg,gest)")

## Define a function that uses information from one row of the "animal_pods"
## data frame to calculate a human equivalent dose (HED).
calc_animal_hed <- function(input_row){
    dm_name_text <- input_row["Dose.Metric.Selected"][[1]]
    if(is.na(dm_name_text)){
        return(NA)
    }
    dm_name <- dm_list[[dm_name_text]]
    dm_val_text <- input_row[input_row["POD"][[1]]][[1]]
    dm_val <- as.numeric(tail(strsplit(dm_val_text," ")[[1]],1))
    chem <- input_row["PFOA.PFOS"][[1]]
    if(dm_name%in%c("C7avg","AUCavg_pup_total")){
        if(input_row["Run.Name"][[1]]=="Zhang et al., 2020, 6505878, Length of Diestrus, C7avg, 1 sd"){
            # 28-d study, not at steady-state
            ## hed <- rev_internal_dose_adult_28d(dm_name,dm_val,age=25,sex="Female",chem=chem)
            hed <- rev_css(dm_val,chem)
            return(hed)
        } else {
            hed <- rev_css(dm_val,chem)
            return(hed)
        }
    } else if (dm_name%in%c("AUCavg_dam_gest","AUCavg_pup_gest","AUCavg_pup_lact","AUCavg_pup_gest_lact",
                            "Cavg_pup_gest")){
        if (grepl("Male",input_row["Animal.Group.Name"][[1]],fixed=T)){
            hed_m <- rev_internal_dose_child(dm_name,dm_val,age=1,sex="Male",chem=chem)
        } else {hed_m <- NULL}
        if (grepl("Female",input_row["Animal.Group.Name"][[1]],fixed=T)){
            hed_f <- rev_internal_dose_child(dm_name,dm_val,age=1,sex="Female",chem=chem)
        } else {hed_f <- NULL}
        hed <- min(hed_m,hed_f)
        return(hed)
    } else if (input_row["Run.Name"][[1]]%in%c("Butenhoff et al., 2012, 2919192, Leydig Cell Adenomas in the Testes, AUC (mg*d/L), 10% rd","Butenhoff et al., 2012, 2919192, Leydig Cell Adenomas in the Testes, AUC (mg*d/L), 4% rd")) {
        hed <- rev_css(dm_val/(365*2),chem)# Convert AUC to average C, then compare to human Css
        return(hed)
    } else {
        print(input_row["Run.Name"][[1]])
        return(NA)
    }
}

## Calculate one HED (mg/kd/g) for each row of the data frame "animal_pods" and
## store the results in a vector.
animal_HEDs <- apply(animal_pods,1,calc_animal_hed)# mg/kg/d

## Read spreadsheet that contains point of departure (POD) internal dose metrics
## derived from epidemiological studies and store PODs and related information in a data
## frame called "human_pods".
human_pods <- read.xlsx("./PODHED_BMD Modeling Summary-211105.xlsx",sheetIndex=2,colIndex=2:17,rowIndex=5:22)## human pods in ng/mL or ug/L

## Define a function that uses information from one row of the "human_pods"
## data frame to calculate a human equivalent dose (HED).
calc_human_hed <- function(input_row){
    chem <- input_row["PFOA.PFOS"][[1]]
    if(input_row["Short.Citation"][[1]]=="Shearer et al. (2021) 7161466"){
        hed <- as.numeric(input_row["Cancer.Slope.Factor."][[1]])/rev_css(1/1000,chem)
        ## dose in mg/kg/d that results in a 1 ng/mL increase in Css
        return(hed)
    }
    dm_text <- input_row[input_row["POD"][[1]]][[1]]
    dm_val <- as.numeric(dm_text)/1000# Conversion from mcg/L(equal to ng/mL) to mg/L.
    if(input_row["Short.Citation"][[1]]=="Budtz-Jorgensen et al. (2018) 5083631"){
        hed_f <- child_rev(dm_val,age=5,sex="Female",chem=chem)#mg/kg/d
        hed_m <- child_rev(dm_val,age=5,sex="Male",chem=chem)#mg/kg/d
        hed <- min(hed_f,hed_m)
        return(hed)
    } else if(input_row["Short.Citation"][[1]]=="Chu et al. (2020) 6315711") {# Maternal blood ad delivery
        hed_f <- child_rev(dm_val,age=0,sex="Female",chem=chem,maternal=T)#mg/kg/d
        hed_m <- child_rev(dm_val,age=0,sex="Male",chem=chem,maternal=T)#mg/kg/d
        hed <- min(hed_f,hed_m)
        return(hed)
    } else if(input_row["Short.Citation"][[1]]=="Govarts et al. (2016) 3230364") {# Cord blood
        hed_f <- child_rev(dm_val,age=0,sex="Female",chem=chem,maternal=F)#mg/kg/d
        hed_m <- child_rev(dm_val,age=0,sex="Male",chem=chem,maternal=F)#mg/kg/d
        hed <- min(hed_f,hed_m)
        return(hed)
    } else if(input_row["Short.Citation"][[1]]=="Sagiv et al. (2018) 4238410") {# Median GA 9 weeks
        age = -0.75*(1-(9/39))# Length of gestation (yr) * Fraction of GA remaining (0.77)
        hed_f <- child_rev(dm_val,age=age,sex="Female",chem=chem,maternal=T)#mg/kg/d
        hed_m <- child_rev(dm_val,age=age,sex="Male",chem=chem,maternal=T)#mg/kg/d
        hed <- min(hed_f,hed_m)
        return(hed)
    } else if(input_row["Short.Citation"][[1]]=="Starling et al. (2017) 3858473") {# Median GA 27 weeks
        age = -0.75*(1-(27/39))# Length of gestation (yr) * Fraction of GA remaining (0.74)
        hed_f <- child_rev(dm_val,age=age,sex="Female",chem=chem,maternal=T)#mg/kg/d
        hed_m <- child_rev(dm_val,age=age,sex="Male",chem=chem,maternal=T)#mg/kg/d
        hed <- min(hed_f,hed_m)
        return(hed)
    } else if(input_row["Short.Citation"][[1]]=="Wikström et al. (2020) 6311677") {# Median GA 10 weeks
        age = -0.75*(1-(10/39))# Length of gestation (yr) * Fraction of GA remaining (0.74)
        hed_f <- child_rev(dm_val,age=age,sex="Female",chem=chem,maternal=T)#mg/kg/d
        hed_m <- child_rev(dm_val,age=age,sex="Male",chem=chem,maternal=T)#mg/kg/d
        hed <- min(hed_f,hed_m)
        return(hed)
    } else {
        hed <- rev_css(dm_val,chem)
        return(hed)
    }
    return("HED not calculated")
}

## Calculate one HED (mg/kd/g) for each row of the data frame "human_pods" and
## store the results in a vector.
human_HEDs <- apply(human_pods,1,calc_human_hed)# mg/kg/d

## Save animal and human HED vectors to 2 sheets in an excel file
write.xlsx(animal_HEDs,file="./HEDs.xlsx",col.names=FALSE,row.names=FALSE,sheetName="Animal HEDs")# units of mg/kg/d
write.xlsx(human_HEDs,file="./HEDs.xlsx",col.names=FALSE,row.names=FALSE,sheetName="Human HEDs",append=T)# units of mg/kg/d
