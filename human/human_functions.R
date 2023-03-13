##------------------------------------------------------------------------------
## human_functions.R - new code for human model, new parameters are in p_PFOA and p_PFOS,
## houses all the main functions to run the model, though for developmental scenarios the
## root function remains in sim_functions.
## For the adult simulations, all the code to run the model is in this file.
##
##------------------------------------------------------------------------------

## Function to model a human developmental scenario. Output is a dataframe of model output. Also houses main parameter definitions for PFOA and PFOS which are passed to human_mother_child, which is the main function to run the model (located in human_functions.R)
human_child <- function(child_age=3,sex=NULL,chem=NULL,dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv"){
    if(is.null(sex)){
        print("Sex of child not specified.")
        return(NULL)
    }
    if (chem=="PFOA"){
        ## p_PFOA=c("CL"=0.085/1000,"Vd"=0.170, "CLr"=0.085/1000, "P_milk"=0.0577, "r_f_m"=0.783, "pow_milk"=1) # PFOA parameters
        ## Vd = 170 ml/kg (Thompson 2010)
        ## t1/2 = 2.7 (Li 2017)
        ## Cl = Vd * ln(2) / t1/2 = 170 mL/kg * ln(2) / (2.7 * 365 d/y) = 0.120 ml/kg/d /1000 = L/kg/d
        ## See "Milk_Plasma_Partition.xlsx" for P_milk
        ## See "Cord_Mother_Ratios.xlsx" for r_f_m
        p_PFOA=c("CL"=0.120/1000,"Vd"=0.170, "CLr"=0.120/1000, "P_milk"=0.049, "r_f_m"=0.83, "pow_milk"=1) # PFOA parameters
        p=c(p_PFOA,"d_m"=dose)# 1 ug/kg/d
    } else if (chem=="PFOS"){
        ## p_PFOS=c("CL"=0.081/1000,"Vd"=0.230, "CLr"=0.081/1000, "P_milk"=0.0138, "r_f_m"=0.454, "pow_milk"=1) # PFOS parameters
        ## Vd = 230 ml/kg (Thompson 2010)
        ## t1/2 = 3.4 (Li 2018)
        ## Cl = Vd * ln(2) / t1/2 = 0.230 mL/kg * ln(2) / (3.4 y * 365 d/y) = 0.128 ml/kg/d
        ## See "Milk_Plasma_Partition.xlsx" for P_milk
        ## See "Cord_Mother_Ratios.xlsx" for r_f_m
        p_PFOS=c("CL"=0.128/1000,"Vd"=0.230, "CLr"=0.128/1000, "P_milk"=0.0160, "r_f_m"=0.40, "pow_milk"=1) # PFOS parameters
        p=c(p_PFOS,"d_m"=dose)# 1 ug/kg/d
    } else {
        print("Specify PFOA or PFOS.")
        return(NULL)
    }
    if(child_age<=0){
        run_child_age=1
    } else {
        run_child_age=child_age
    }
    out <- human_mother_child(age_end=run_child_age+25,p=p,sex=sex,mName="pfas_tk",milkin=milkin,waterin_m=waterin_m,waterin_i=waterin_i,DW=DW)
    out <- cbind(out,(out[,"time"]/365)-25)
    colnames(out)[ncol(out)] <- "years"
    if(child_age<=0){
        out <- out[out[,"years"]<=child_age,]
    }
    return(out)
}

## Function to predict a maternal or child concentration for a developmental scenario. Calls human_child and outputs the last value for child (C_i) or maternal (C_m) serum concentration.
human_child_conc <- function(child_age=3,sex=NULL,chem=NULL,dose=NULL,maternal=F,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv"){
    out <- human_child(child_age,sex=sex,chem=chem,dose=dose,DW=DW,waterin_m=waterin_m,waterin_i=waterin_i,milkin=milkin)
    out_C_i <- tail(out[,"C_i"],1)
    out_C_m <- tail(out[,"C_m"],1)
    ##print(sprintf('Child Plasma concentration: %.3f mg/L\nMaternal Plasma concentration: %.3f mg/L',out_C_i,out_C_m))
    if(maternal){
        return(out_C_m)
    }
    return(out_C_i)
}

## Adult equivalent of human_child and human_mother_child. Sets up the simulation (including parameters) for an adults, sets up the model input, runs the model, and outputs a dataframe of model output.
human_adult <- function(age=30,sex=NULL,chem=NULL,dose=0.001,dose_change=NULL,DW=FALSE){# 0 initial PFAS concentration
    if(is.null(sex)){
        print("Sex not specified.")
        return(NULL)
    }
    if (chem=="PFOA"){
        ## p_PFOA=c("CL"=0.085/1000,"Vd"=0.170, "CLr"=0.085/1000, "P_milk"=0.0577, "r_f_m"=0.783, "pow_milk"=1) # PFOA parameters
        p_PFOA=c("CL"=0.120/1000,"Vd"=0.170, "CLr"=0.120/1000, "P_milk"=0.049, "r_f_m"=0.83, "pow_milk"=1) # PFOA parameters
        p=c(p_PFOA,"d_m"=dose)# 1 ug/kg/d
    } else if (chem=="PFOS"){
        ## p_PFOS=c("CL"=0.081/1000,"Vd"=0.230, "CLr"=0.081/1000, "P_milk"=0.0138, "r_f_m"=0.454, "pow_milk"=1) # PFOS parameters
        p_PFOS=c("CL"=0.128/1000,"Vd"=0.230, "CLr"=0.128/1000, "P_milk"=0.0160, "r_f_m"=0.40, "pow_milk"=1) # PFOS parameters
        p=c(p_PFOS,"d_m"=dose)# 1 ug/kg/d
    } else {
        print("Specify PFOA or PFOS.")
        return(NULL)
    }
    if(is.null(sex)){
        print("Sex not specified.")
        return(NULL)
    }
    t_start = 0
    t_end = age*365 # End of simulation (days)
    times = seq(from=t_start, to=t_end, by=0.25)  # Sequence for model results
                                        # Set parameters and initial values of state variables.
    p=init_human_p(p)
    if(DW){
        p["water_dose"]=TRUE
    }
    parms = initParms(newParms = c(p[c("CL","Vd","F_abs","P_milk","r_f_m","t_con",
                                       "water_dose","pow_dose","pow_milk")],"n_i"=1))
    Y0 = initStates(parms)

    Y0["d_m"] = p[["d_m"]] # Dose (mg/kg/d)
    Y0["A_mf"] = p[["A_m_init"]]

    if (sex=="Female"){
        M_m_df1 = read.csv("human_female_mass_0_to_20.csv", header=FALSE,
                           col.names=c("Age (y)", "Mass (kg)"),
                           check.names=FALSE)

        M_m_df2 = read.csv("human_female_mass_20_to_75.csv", header=FALSE,
                           col.names=c("Age (y)", "Mass (kg)"),
                           check.names=FALSE)
    } else if (sex=="Male"){
        M_m_df1 = read.csv("human_male_mass_0_to_20.csv", header=FALSE,
                           col.names=c("Age (y)", "Mass (kg)"),
                           check.names=FALSE)
        M_m_df2 = read.csv("human_male_mass_20_to_75.csv", header=FALSE,
                           col.names=c("Age (y)", "Mass (kg)"),
                           check.names=FALSE)
    } else {
        stop("Failed to assign BW.")
    }

    M_m_df = cbind(times=c(M_m_df1[["Age (y)"]],M_m_df2[["Age (y)"]])*365, mass=c(M_m_df1[["Mass (kg)"]],M_m_df2[["Mass (kg)"]]))

    M_mf_in = M_m_df


                                        # Infant mass interpolation points (zero weight).
    t_m_vec = M_m_df[,"times"] * 365
    tparms = c(min(t_m_vec), max(t_m_vec))
    M_i_in = cbind(times=tparms, mass=c(0, 0))

    R_milk_in = cbind(times=tparms,c(0,0))

                                        # For the following arrays, since the value is constant (up to child-birth),
                                        #the exact time range doesn't matter, as the interpolation method will
                                        # extrapolate the same value to all other times.

                                        # Human breast milk ingestion rate = 0, interpolation points to birth.
    R_milk_in = cbind(times=c(0,1), rate=c(0, 0))

                                        # Infant dose = 0, interpolation points to birth.
    d_i_in = cbind(times=tparms, dose=c(0, 0))

                                        # Ratio of concentration in mother to mother + fetus, assumed = 1
                                        # But giving 3 time-points/values, since departure from 1 would only occur
                                        # after conception.
    r_m_mf_in = cbind(times=c(0,24.25*365,t_end), ratio=c(1,1,1))

                                        # Clearance for menstruating women: menstruation assumed to start at age 12.4 y
                                        # and to cease from conception (24.25 y) until the end of breastfeeding (26 y),
                                        # and to cease with menopause at age 50 y. p[["CLr"]] is the average clearance
                                        # rate when a woman is menstruating.
                                        # Set up data frame for control of CLr in model:
    df_CLr = data.frame(var=rep("CLr",5), # CLrv has 5 entries
                        time  = c(0,       12.4, 24.25,         26, 50)*365,
                        value = c(0, p[["CLr"]],     0, p[["CLr"]],  0),
                        method=rep("replace",5))

    Vdaf_input = read.csv("vd_child_con.csv", header=FALSE,# Constant vd=adult in children
    ## Vdaf_input = read.csv("vd_child.csv", header=FALSE,
                       col.names=c("Age (y)", "Vdaf"),
                       check.names=FALSE)


    Vdaf_m_in = cbind(times = Vdaf_input[,1]*365,Vdaf_m_in = Vdaf_input[,2])

    Vdaf_i_in = cbind(time=c(0,t_end), Vdaf_i_in = c(1,1))

    DW_input = read.csv("water_95.csv", header=FALSE,
                       col.names=c("Age (y)", "DW (ml/kg/d)"),
                       check.names=FALSE)

    DW_m_in = cbind(times = DW_input[,1]*365,DW_m_in = DW_input[,2])

    DW_i_in = cbind(time=c(0,t_end), DW_i_in = c(0,0))

    df_events = rbind(df_CLr)

    if(is.null(dose_change)){
        ## Run model simulation
        out = run_model("pfas_tk", times, Y0, parms, rtol=1e-6, atol=1e-6,
                        list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in, d_i_in,
                             Vdaf_m_in, Vdaf_i_in, DW_m_in, DW_i_in),
                        list(method="linear", rule=2, ties="ordered"),
                        list(data=df_events))
    } else {
        dose_events <- data.frame(var=rep("d_m",nrow(dose_change)),time=dose_change$time*365,value=dose_change$d_m,method=rep("replace",nrow(dose_change)))
        event_list <- rbind(df_events,dose_events)
        event_list_sorted <- event_list[order(event_list$time),]
        event_list_sorted <- event_list_sorted[event_list_sorted$time<=t_end,]
        times <- cleanEventTimes(times,event_list_sorted$time)
        times <- sort(c(times,event_list_sorted$time))
        out = run_model("pfas_tk", times, Y0, parms,rtol=1e-6, atol=1e-6,
                        list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in, d_i_in,
                             Vdaf_m_in, Vdaf_i_in, DW_m_in, DW_i_in),
                        list(method="linear", rule=2, ties="ordered"),
                        list(data=event_list_sorted))
    }

    out <- cbind(out,(out[,"time"]/365))
    colnames(out)[ncol(out)] <- "years"
    return(out)
}

## Outputs the final serum concentration for a non-developmental scenario. Calls human_adult.
human_adult_conc <- function(age=25,sex=NULL,chem=NULL,dose=NULL,dose_change=NULL,DW=FALSE){
    out <- human_adult(age,sex=sex,chem=chem,dose=dose,dose_change=dose_change,DW=DW)
    out_C_m <- tail(out[,"C_m"],1)
    ##print(sprintf('Adult Plasma concentration: %.3f mg/L',out_C_m))
    return(out_C_m)
}

## Reverse dosimitry for a serum concentration for a non-developmental scenario. Target_conc is in mg/L. Since the model is linear with dose the output dose is calculated by multiplying the target_conc by the ratio of a standard dose and standard concentration. The standard dose is 1 mg/kg/d, so it is implicit when dose is set.
adult_rev <- function(target_conc,age=25,sex=NULL,chem=NULL,DW=FALSE){
    standard_conc <- human_adult_conc(age=age,sex=sex,chem=chem,dose=1,DW=DW)
    dose <- target_conc/standard_conc
    return(dose)
}

## Reverse dosimitry for a serum concentration for a developmental scenario. Can be based on a maternal or child concentration (set with the maternal flag). Target_conc is in mg/L. Since the model is linear with dose the output dose is calculated by multiplying the target_conc by the ratio of a standard dose and standard concentration. The standard dose is 0.001 mg/kg/d, and it is explicit when dose is set. 0.001 mg/kg/d is used because 1 mg/kg/d is a very large dose for PFOA/S, especially for developmental scenarios and resulted in numerical instability (failure in model convergance).
child_rev <- function(target_conc,age=5,sex=NULL,chem=NULL,maternal=F,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv"){
    dose <- 0.001
    standard_conc <- human_child_conc(child_age=age,sex=sex,chem=chem,dose=dose,maternal=maternal,DW=DW,waterin_m=waterin_m,waterin_i=waterin_i,milkin=milkin)
    dose <- target_conc/standard_conc*dose
    return(dose)
}

## Calculates a steady state serum concentration (Css) using a clearance value (CL). Css = dose/CL. Dose in mg/kg/d.
calc_css <- function(dose=0.001,chem="PFOA"){
    if (chem=="PFOA"){
        CL=0.120/1000
    } else if (chem=="PFOS"){
        CL=0.128/1000
    } else {
        print("Specify PFOA or PFOS.")
        return(NULL)
    }
    Css <- dose/CL
    return(Css)# mg/L
}

## Reverse dosimetry for steady-state serum concentration (Css) using a clearance value (CL) dose = Css*Cl. Dose is in mg/kg/d, Css is in mg/L.
rev_css <- function(Css,chem="PFOA"){
    if (chem=="PFOA"){
        CL=0.120/1000
    } else if (chem=="PFOS"){
        CL=0.128/1000
    } else {
        print("Specify PFOA or PFOS.")
        return(NULL)
    }
    dose <- Css*CL
    return(dose)# mg/L
}

## Function to convert an excel column value to a numerical value (A-Z are 1-26, AA-AZ are 27 to 52, BA-BZ are 53-78, and so on). This was used to identify which columns in HED_calc.R and HED_calc_update.R to read, but isn't called by any script or other function.
string_index <- function(string_ind){
    az_ind <- sapply(strsplit(string_ind,""),function(x){print(x);match(x,LETTERS)})
    num_ind <- rep(0,length(az_ind))
    for(i in 1:length(az_ind)){
            num_ind[i] <- rev(az_ind)[i]*(26^(i-1))
    }
    return(sum(num_ind))
}

## Reverse dosimetry for a developmental scenario. Accepts a dose metric name (dose_metric) and value (dose_metric_value), along with the age, sex of the child, and chemical (PFOA or PFOS). The dose_metric (see calc_internal_dose) accepts vectors from the model output from human_child. As in the case of reverse dosimetry, the dose-metric values is multiplied by the ratio of a standard dose (1 mg/kg/d) and a standard dose metric value (DM).
rev_internal_dose_child <- function(dose_metric, dose_metric_value, age, sex, chem, DW=FALSE){
    human_data <- human_child(child_age=age,sex=sex,chem=chem,dose=1,DW=DW)
    DM <- calc_internal_dose(human_data[,"time"],dose_metric,
                             C_m=human_data[,"C_m"],C_i=human_data[,"C_i"],
                             AUC_m=human_data[,"AUC_m"],AUC_i=human_data[,"AUC_i"],
                             t_conc=(24.25*365),t_birth=(25*365),t_wean=(26*365))
    dose <- as.numeric(dose_metric_value)/DM
    return(dose)
}

## Reverse dosimetry for a non-developmental scenario. Accepts a dose metric name (dose_metric) and value (dose_metric_value), along with the age, sex of the individual, and chemical (PFOA or PFOS). The dose_metric (see calc_internal_dose) accepts vectors from the model output from human_adult. As in the case of reverse dosimetry, the dose-metric values is multiplied by the ratio of a standard dose (1 mg/kg/d) and a standard dose metric value (DM).
rev_internal_dose_adult <- function(dose_metric, dose_metric_value, age, sex, chem, DW=FALSE){
    human_data <- human_adult(age=age,sex=sex,chem=chem,dose=1,DW=DW)
    DM <- calc_internal_dose(human_data[,"time"],dose_metric,C_m=human_data[,"C_m"],AUC_m=human_data[,"AUC_m"])
    dose <- as.numeric(dose_metric_value)/DM
    return(dose)
}

## No longer used
rev_internal_dose_adult_28d <- function(dose_metric, dose_metric_value, age, sex, chem,DW=FALSE){
    human_data <- human_adult(age=age,sex=sex,chem=chem,dose=0,dose_change=data.frame(d_m=1,time=(age-(28/365))),DW=DW)
    DM <- calc_internal_dose(human_data[,"time"],dose_metric,C_m=human_data[,"C_m"],AUC_m=human_data[,"AUC_m"])
    dose <- as.numeric(dose_metric_value)/DM
    return(dose)
}

## Function that outputs a dose metric value (metric) and accepts as input, vectors from model output from a developmental or adult scenario.
calc_internal_dose <- function(times, dose_metric, C_m=NULL, C_i=NULL, AUC_m=NULL, AUC_i=NULL,
                               t_conc=NULL, t_birth=NULL, t_wean=NULL, days=7){
### times in days, dose_metric of interest, maternal concentration vector, infant concentration vector,
### time of conception, time of birth, time of weaning, days for C7avg DM.
    if (dose_metric == 'C7avg'){ ### Adult
        end = length(times)
        start = which(times==(times[end]-7))
        metric = (AUC_m[end]-AUC_m[start])/(times[end]-times[start])
    } else if (dose_metric == 'AUC'){
        metric = tail(AUC_m,1)
    } else if (dose_metric == 'Cavg'){
        metric = tail(AUC_m,1)/tail(times,1)
    } else if (dose_metric == 'Cmax'){
        metric = max(C_m)
    } else if (dose_metric == 'AUCavg'){ ### Dam Developmental
        metric = (AUC_m[times==t_birth]-AUC_m[times==t_conc])/(t_birth-t_conc)
    } else if (dose_metric == 'Cmax_pre'){
        metric = max(C_m[times<t_conc])
    } else if (dose_metric == 'AUCavg_pre'){
        metric = AUC_m[times==t_conc]/times[times==t_conc]
    } else if (dose_metric == 'Cmax_dam'){
        metric = max(C_m[times>=t_conc&times<=t_birth])
    } else if (dose_metric == 'Cavg_dam'){
        metric = (AUC_m[times==t_birth]-AUC_m[times==t_conc])/(t_birth-t_conc)
    } else if (dose_metric == 'AUCavg_dam_gest'){
        metric = (AUC_m[times==t_birth]-AUC_m[times==t_conc])/(t_birth-t_conc)
    } else if (dose_metric == 'AUCavg_dam_lact'){
        metric = (AUC_m[times==t_wean]-AUC_m[times==t_birth])/(t_wean-t_birth)
    } else if (dose_metric == 'AUCavg_dam_gest_lact'){
        metric = (AUC_m[times==t_wean]-AUC_m[times==t_conc])/(t_wean-t_conc)
    } else if (dose_metric == 'Cmax_pup'){ ### Pup Developmental
        metric = max(C_i)
    } else if (dose_metric == 'Cavg_pup_gest'){ ### Pup Developmental
        metric = (AUC_i[times==t_birth]-AUC_i[times==t_conc])/(t_birth-t_conc)
    } else if (dose_metric == 'Cmax_pup_gest'){
        metric = max(C_i[times>=t_conc&times<=t_birth])
    } else if (dose_metric == 'Cavg_pup_lact'){
        metric = (AUC_i[times==t_wean]-AUC_i[times==t_birth])/(t_wean-t_birth)
    } else if (dose_metric == 'Cmax_pup_lact'){
        metric = max(C_i[times>=t_birth&times<=t_wean])
    } else if (dose_metric == 'Cavg_pup_gest_lact'){
        metric = (AUC_i[times==t_wean]-AUC_i[times==t_conc])/(t_wean-t_conc)
    } else if (dose_metric == 'Cmax_pup_wean'){
        metric = max(C_i[times>=t_wean])
    } else if (dose_metric == 'AUCavg_pup_gest'){
        metric = AUC_i[times==t_birth]/(t_birth-t_conc)
    } else if (dose_metric == 'AUCavg_pup_lact'){
        AUC_gest = AUC_i[times==t_birth]
        AUC_lact = AUC_i[times==t_wean]
        metric = (AUC_lact-AUC_gest)/(t_wean-t_birth) # Divide by total days of lactation
    } else if (dose_metric == 'AUCavg_pup_gest_lact'){
        metric = AUC_i[times==t_wean]/(t_wean-t_conc)
    } else if (dose_metric == 'AUCavg_pup_diet'){
        #t_diet = tail(times,1)
        metric = tail(AUC_i,1)/tail(times,1)
    } else if (dose_metric == 'AUCavg_pup_total'){
        metric = tail(AUC_i,1)/tail(times,1)
    }
    return(metric)
}
