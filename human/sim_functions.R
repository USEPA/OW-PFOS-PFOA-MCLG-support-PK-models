#------------------------------------------------------------------------------
# sim_functions.R
#
# Original file housing main functions to run the model.
# Minimal changes from lipophilic_tk model (Kapraun, 2022).
# Includes some legacy functions that are no longer used.
#
# References:
# -- Kapraun, D.F. et al. (2022). A Generic Pharmacokinetic Model for
#    Quantifying Mother-to-Offspring Transfer of Lipophilic Persistent
#    Environmental Chemicals.
#
#------------------------------------------------------------------------------

## Legacy function
## bw_int<-function(tdoses=1, doses=c(0,10), BWs=rbind(c(250,250),c(250,250))) {
##   # number of doses should be # of columns in BWs
##   bw_int=c()
##   bwd=BWs[,1]
##   for (td in tdoses) {
##     for (i in 1:length(bwd)) {bwd[i]=approx(doses,BWs[i,],td)$y}
##     bw_int=cbind(bw_int,bwd)
##   }
##   return(bw_int)
## }

## Legacy function
## cavg_nonpreg<-function(rparms=NULL, times=-28:0, doses=1,
##                        BWs=rbind(c(0,250),c(28,250))) {
##   tparms=c(t_m_start = min(times), # Time (d) since conception "maternal" dose starts.
##            t_m_end = max(times))  # Time (d) since conception "maternal" dose ends.

##   p=init_rat_p(rparms)
##   parms = initParms(p[c("CL","Vd","F_abs","P_milk","r_f_m","t_con",
##                         "food_dose","pow_dose","pow_milk","n_i")])
##                        # Update calculated parameters.

##   t_m_vec = BWs[,1]+min(times)
##   nt = length(t_m_vec)

##   # Rat infant mass interpolation points (zero weight).
##   M_i_in = cbind(times=tparms, mass=c(0, 0))

##   # Rat breast milk ingestion rate interpolation points (zero ingestion).
##   R_milk_in = cbind(times=tparms, rate=c(0, 0))

##   # Infant dose = 0, interpolation points to birth.
##   d_i_in = cbind(times=c(0,1), dose=c(0, 0))

##   # Clearance for menstruating women, not used for rat.
##   # Set up data frame for control of CLr in model:
##   df_CLr = data.frame(var=c("CLr","CLr"), # CLr has 2 entries
##                       time  = c(-t_exp, 0),
##                       value = c(0,  0),
##                       method=rep("replace",2))

##   # Nominal ratio of concentrations in mother and mother + fetus(es) is 1.0 for simulation.
##   r_m_mf_in = cbind(times=tparms, ratio=c(1,1))

##   # Get initial values of state variables.
##   Y0 = initStates(parms)
##   cavg=NULL
##   for(n in 1:length(doses)) {
##     # Rat mass interpolation points.
##     M_mf_in = cbind(times=BWs[,1], mass=BWs[,n+1]/1000)
##       # Body weights in column n+1, g => kg

##     Y0[["d_m"]] = doses[n]

##     # Run model simulation to birth
##     out1 = run_model("pfas_tk", times, Y0, parms, rtol=1e-6, atol=1e-6,
##                      list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in, d_i_in),
##                      list(method="linear", rule=2, ties="ordered"),
##                      list(data=df_CLr))
##     cavg = c(cavg, mean(out1[,"C_m"]))
##   }
##   plot(out1[,"time"],out1[,"C_m"],type="l")
##   plot(out1[,"time"],out1[,"d_m"],type="l")
##   cavg_nonpreg=cbind(doses,cavg)
## }

## Legacy function
## cavg_preg<-function(rparms=NULL, dstart=6, dend=15, times=0:18, doses=1,
##                     BWs=rbind(c(0,250),c(28,250))) {
##   tparms=c(t_m_start = min(times), # Time (d) since conception simulation starts.
##            t_m_end = max(times))  # Time (d) since conception simulation ends.

##   p=init_rat_p(rparms)
##   parms = initParms(p[c("CL","Vd","F_abs","P_milk","r_f_m","t_con",
##                         "food_dose","pow_dose","pow_milk","n_i")])
##                         # Update calculated parameters.

##   t_m_vec = min(BWs[,1]):max(BWs[,1]) # Time vector of BWs

##   nt = length(t_m_vec)

##   # Rat infant mass interpolation points (zero weight).
##   M_i_in = cbind(times=tparms, mass=c(0, 0))

##   # Rat breast milk ingestion rate interpolation points (zero ingestion).
##   R_milk_in = cbind(times=tparms, rate=c(0, 0))

##   # Infant dose = 0, interpolation points to birth.
##   d_i_in = cbind(times=c(0,1), dose=c(0, 0))

##   # Clearance for menstruating women, not used for rat.
##   # Set up data frame for control of CLr in model:
##   df_CLr = data.frame(var=rep("CLr",2), # CLr has 2 entries
##                       time  = c(t_m_start, t_m_end),
##                       value = c(0,  0),
##                       method=rep("replace",2))

##   # Nominal ratio of concentrations in mother and mother + fetus(es) is 1.0 for simulation.
##   r_m_mf_in = cbind(times=tparms, ratio=c(1,1))

##   # Get initial values of state variables.
##   Y0 = initStates(parms)
##   cavg=NULL
##   for(n in 1:length(doses)) {
##     # Rat mass interpolation points.
##     M_mf_in = cbind(times=BWs[,1], mass=BWs[,n+1]/1000)
##     # Body weights in column n+1, g => kg

##     Y0[["d_m"]] = 0 # Dosing begins at zero
##     # Dosing events. Turn dosing on and off at specified times.
##     df_d_m = data.frame(var=rep("d_m",2), # 2 changes, on and off
##                         time  = c(dstart, dend),
##                         value = c(doses[n],  0),
##                         method=rep("replace",2))
##     # Combined events
##     df_c = rbind(df_CLr,df_d_m)
##     df_c = df_c[order(df_c["time"]), ]

##     # Run simulation
##     out1 = run_model("pfas_tk", times, Y0, parms, rtol=1e-6, atol=1e-6,
##                      list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in, d_i_in),
##                      list(method="linear", rule=2, ties="ordered"),
##                      list(data=df_c))
##     cavg = c(cavg, mean(out1[(times>=6),"C_m"]))
##   }
##   plot(out1[,"time"],out1[,"C_m"],type="l")
##   plot(out1[,"time"],out1[,"d_m"],type="l")
##   plot(out1[,"time"],out1[,"M_mf"],type="l")
##   cavg_nonpreg=cbind(doses,cavg)
## }

init_human_p <- function(new_p = NULL) {
  # Default parameter values for human simulations.
  p = c(t_con = 24.25 * 365, # Age (d) of conception.
        t_gest = 0.75 * 365, # Duration (d) of pregnancy/gestation.
        t_lact = 1.0 * 365,  # Duration (d) of lactation/nursing.
        CL = 0.066,   # Clearance (mL/kg/d) of substance in this species.
        CLr = 0.095,  # Clearance (mL/kg/d) for reproductive-age women
        Vd = 380,     # Volume of distribution (mL/kg) of substance in this species.
        P_milk = 0.03, # Milk/maternal blood concentration partition coefficient.
        r_f_m = 1.0,   # Ratio of fetal and maternal concentrations.
        F_abs = 1.0,   # Fraction of nominal dose that is absorbed.
                                        #food_dose = 0,  # A boolean flag. Dose is delivered in food if nonzero.
        water_dose = 0,  # A boolean flag. Dose is delivered in food if nonzero.
        d_m = 0,        # Dose rate (mg/kg/d or mg/kg) to mother.
        pow_dose = 1.0, # When 1, d_m is mg/kg/d, when 0, d_m is mg/d
        pow_milk = 0.0, # when 1, milk transfer rate is kg/kg/d, when 0, is kg/d
        A_m_init = 0)   # Initial amount (mg) in mother.

  # If "new_p" is not NULL, replace default parameter values with values from
  # that vector.
  if (!is.null(new_p)) {
    if (!all(names(new_p) %in% c(names(p)))) {
      stop("Illegal parameter name in 'new_p'.")
    } else {
      p[names(new_p)] <- new_p
    }
  }
  return(p)
}

human_mother_child <- function(age_start=0, a_start=0, age_end=30, p=NULL, sex="Female",
                               mName="pfas_tk", milkin="r_milk_bw_95.csv",
                               waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",DW=FALSE) {
  t_start = age_start*365
  t_end = 25*365 # 25 years of age, end of pregnancy, converted to days
  times = seq(from=t_start, to=t_end, by=0.25)  # Sequence for model results
  times = sort(unique(c(times,(25-(0.75*(seq(1,38,1)/39)))*365)))# For weeks of pregnancy
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


  M_m_df1 = read.csv("human_female_mass_0_to_25_preg.csv", header=FALSE,
                    col.names=c("Age", "Mass (kg)"),
                    check.names=FALSE)

  M_m_df1[,1] = M_m_df1[,1]*365 # Convert age from years to days

  M_m_df2 =read.csv("human_female_mass_25_to_35_postnatal.csv",
                  header=FALSE, col.names=c("Age", "Mass (kg)"),
                  check.names=FALSE)
  M_m_df2[,1] = M_m_df2[,1]*365 # Convert age from years to days

  if (sex=="Female"){
      M_i_df1 = read.csv("human_female_mass_0_to_20.csv", header=FALSE,
                         col.names=c("Age (y)", "Mass (kg)"),
                         check.names=FALSE)

      M_i_df2 = read.csv("human_female_mass_20_to_75.csv", header=FALSE,
                         col.names=c("Age (y)", "Mass (kg)"),
                         check.names=FALSE)
  } else if (sex=="Male"){
      M_i_df1 = read.csv("human_male_mass_0_to_20.csv", header=FALSE,
                         col.names=c("Age (y)", "Mass (kg)"),
                         check.names=FALSE)
      M_i_df2 = read.csv("human_male_mass_20_to_75.csv", header=FALSE,
                         col.names=c("Age (y)", "Mass (kg)"),
                         check.names=FALSE)
  } else {
      stop("Failed to assign filial BW.")
  }

  M_i_df1[,1] = M_i_df1[,1]*365 # Convert age from years to days
  M_i_df2[,1] = M_i_df2[,1]*365 # Convert age from years to days

  M_mf_in = cbind(times=M_m_df1[["Age"]], mass=M_m_df1[["Mass (kg)"]])

  # For the following arrays, since the value is constant (up to child-birth),
  #the exact time range doesn't matter, as the interpolation method will
  # extrapolate the same value to all other times.

  # Infant mass = 0, interpolation points up to moment of birth.
  M_i_in = cbind(times=c(0,1), mass=c(0,0))

  # Human breast milk ingestion rate = 0, interpolation points to birth.
  R_milk_in = cbind(times=c(0,1), rate=c(0, 0))

  # Infant dose = 0, interpolation points to birth.
  d_i_in = cbind(times=c(0,1), dose=c(0, 0))

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

  Vdaf_i_in = cbind(times = (Vdaf_input[,1]+25)*365,Vdaf_i_in = Vdaf_input[,2])

  DW_input = read.csv(waterin_m, header=FALSE,
                      col.names=c("Age (y)", "DW"),
                      check.names=FALSE)

  DW_m_in = cbind(times = DW_input[,1]*365,DW_m_in = DW_input[,2])

  DW_input = read.csv(waterin_i, header=FALSE,
                      col.names=c("Age (y)", "DW"),
                      check.names=FALSE)

  DW_i_in = cbind(times = (DW_input[,1]+25)*365,DW_i_in = DW_input[,2])

  df_events = rbind(df_CLr)

  # Run model simulation to birth
  out1 = run_model(mName, times, Y0=Y0, parms=parms, rtol=1e-6, atol=1e-6,
                   forcing = list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in, d_i_in,
                        Vdaf_m_in, Vdaf_i_in, DW_m_in, DW_i_in),
                   fcontrol=list(method="linear", rule=2, ties="ordered"),
                   event_list=list(data=df_events))

  # Now run post-pregnancy simulations, need to set up initial conditions and
  # time-dependent parameters first.

  t_start2 = 25*365; t_end2 = age_end*365; # Maternal age 25-age_end
                                        # infant age 0 - (age_end - 25)
  times2 = seq(from=t_start2, to=t_end2, by=0.25)
  # Maternal mass interpolation points after moment of birth.
  M_mf_in = cbind(times=M_m_df2[["Age"]], mass=M_m_df2[["Mass (kg)"]])

  # Infant mass interpolation points after moment of birth.
  # Uses same growth vs time as mother, but shift ages to reflect mother's age.
  M_i_in = cbind(times=c(M_i_df1[["Age (y)"]],M_i_df2[["Age (y)"]])+t_start2, mass=c(M_i_df1[["Mass (kg)"]],M_i_df2[["Mass (kg)"]]))

  # Ratio of concentration in mother to mother + fetus, assumed = 1
  r_m_mf_in = cbind(times=c(0,1), ratio=c(1,1))

  # Human breast milk ingestion rate and fraction from food interpolation points.
  r_milk = read.csv(milkin, header=FALSE,
                    col.names=c("Age (d)", "Rate (kg/d)", "Frac food"),
                    check.names=FALSE)

  R_milk_in = cbind(times=(r_milk[["Age (d)"]]+t_start2),
                    rate=r_milk[["Rate (kg/d)"]])

  if (p["P_milk"]>0) { # Lactation delivery occurring
    # Following sets infant direct dosing to ramp up begging at age 90 days to
    # approximately 58% of d_m at 365 days, then jumps to 100% of d_m.
    d_i_in = cbind(times=(r_milk[["Age (d)"]]+t_start2),
                 inf_dose=(r_milk[["Frac food"]]*p[["d_m"]]))
    # d_i_in[c(4,5),2]=c(0,0) # Un-comment to turn off ramp-up of infant direct
    # dose intake, so it just jumps from 0 to 100% at weaning.
  } else {
    d_i_in = cbind(times=c(0,1), inf_dose=(c(0,1)*p[["d_m"]]))
  }

    # Apportion the amount of substance in the mother at the moment of birth to
    # the mother and the infant(s) according to body mass. Use these values to
    # set the initial conditions for the post-birth simulation.

      # Concentration in mother at birth (mg/kg), converted from mg/mL with Vd
  C_m_birth = out1[[nrow(out1), "C_m"]] * parms[["Vd"]]
      # mg/kg in mother & infant after birth
  # C_m_birth = p[["d_m"]]*p[["Vd"]]/p[["CL"]]  # Un-comment to assume SS @ birth

  Y0["A_mf"] = C_m_birth * M_mf_in[1,"mass"]
  Y0["A_i"] =  out1[[nrow(out1), "C_i"]] * parms[["Vd"]] * M_i_in[1,"mass"] *
      Vdaf_input[Vdaf_input[,"Age (y)"]==0,"Vdaf"]
  # Y0["A_i"] =  C_m_birth * p[["r_f_m"]]* M_i_in[1,"mass"] # Un-comment to assume SS @ birth

      # Carry forward AUC and total ingested/excreted/lost calculations
  Y0["AUC_m"] = out1[[nrow(out1), "AUC_m"]]
  Y0["AUC_i"] = out1[[nrow(out1), "AUC_i"]]
  Y0["T_in"] = out1[[nrow(out1), "T_in"]]
  Y0["T_out"] = out1[[nrow(out1), "T_out"]] +
    C_m_birth*out1[[nrow(out1), "M_mf"]] - Y0["A_mf"] - Y0["A_i"]
      # Last part accounts for adds amount lost at birth: total in mother
      # and fetus just before birth minus amount in each just after birth.

  # Run simulation for times after birth.
  out2 = run_model(mName, times2, Y0=Y0, parms=parms, rtol=1e-6, atol=1e-6,
                   forcing=list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in, d_i_in,
                        Vdaf_m_in, Vdaf_i_in, DW_m_in, DW_i_in),
                   fcontrol=list(method="linear", rule=2, ties="ordered"),
                   event_list=list(data=df_events))
  out = rbind(out1[1:nrow(out1)-1, ], out2)

  return(out)   # Return output data frame.
}

## Legacy function
## human_sim <- function(times=NULL, p=NULL, M_m_df1=NULL, M_m_df2=NULL,
##                       M_m_df3=NULL, M_i_df1=NULL, M_i_df2=NULL, preg=TRUE,
##                       mName="pfas_tk") {
##   # If p is NULL, or if not all parameters are provided in p, use default
##   # default values for the parameters not provided.
##   p = init_human_p(p)

##   # If t is NULL, use default times for simulation.
##   if (is.null(times)) {
##     if (preg) {
##       t_start = -24.25 * 365
##       t_end = p[["t_gest"]] + p[["t_lact"]]
##       times = seq(from=t_start, to=t_end, by=0.25)
##     } else {
##       t_start = 0
##       t_end = 78 * 365    # Simulate 78 years.
##       times = seq(from=t_start, to=t_end, by=1)
##     }
##   }

##   # If not provided, obtain time course data for mass of human mother.
##   if (is.null(M_m_df1)){
##     M_m_df1 = read.csv("human_female_mass_0_to_20.csv", header=FALSE,
##                        col.names=c("Age (y)", "Mass (kg)"),
##                        check.names=FALSE)
##   }
##   if (is.null(M_m_df2)) {
##     if (preg) {
##       M_m_df2 = read.csv("human_female_mass_20_to_25_pregnancy.csv",
##                          header=FALSE,
##                          col.names=c("Age (y)", "Mass (kg)"),
##                          check.names=FALSE)
##     } else {
##       M_m_df2 = read.csv("human_female_mass_20_to_75.csv",
##                          header=FALSE,
##                          col.names=c("Age (y)", "Mass (kg)"),
##                          check.names=FALSE)
##     }
##   }
##   if (is.null(M_m_df3) && preg) {
##     M_m_df3 = read.csv("human_female_mass_25_to_26_postnatal.csv",
##                        header=FALSE, col.names=c("Age (y)", "Mass (kg)"),
##                        check.names=FALSE)
##   }

##   # If not provided, obtain time course data for mass of human infant
##   # (through adulthood).
##   if (preg) {
##     if (is.null(M_i_df1)) {
##       M_i_df1 = read.csv("human_female_mass_0_to_20.csv", header=FALSE,
##                          col.names=c("Age (y)", "Mass (kg)"),
##                          check.names=FALSE)
##     }
##     if (is.null(M_i_df2)) {
##       M_i_df2 = read.csv("human_female_mass_20_to_75.csv", header=FALSE,
##                          col.names=c("Age (y)", "Mass (kg)"),
##                          check.names=FALSE)
##     }
##   }

##   # Maternal mass interpolation points up to moment of birth. Assume the
##   # moment of birth never occurs when preg is FALSE.
##   M_m_df = rbind(M_m_df1, M_m_df2)
##   if (preg) {
##     t_m_vec = (M_m_df[["Age (y)"]] - 24.25) * 365
##   } else {
##     t_m_vec = M_m_df[["Age (y)"]] * 365
##   }
##   M_m_vec = M_m_df[["Mass (kg)"]]
##   M_m_in = cbind(times=t_m_vec, mass=M_m_vec)

##   tparms = c(min(t_m_vec), max(t_m_vec))
##   # Infant mass interpolation points up to moment of birth.
##   M_i_in = cbind(times=tparms, mass=c(0, 0))

##   # Human breast milk ingestion rate interpolation points.
##   R_milk_in = cbind(times=tparms, rate=c(0, 0))

##   # Ratio of in mother to mother + fetus
##   r_m_mf_in = cbind(times=tparms, ratio=c(1,1))

##   # Infant dose = 0, interpolation points to birth.
##   d_i_in = cbind(times=c(0,1), dose=c(0, 0))

##   # Clr is cearance for menstruating women.
##   df_CLr = data.frame(var=rep("CLr",3), # CLr has 3 entries
##                       time  = c(0,       12.4, 50)*365,
##                       value = c(0, p[["CLr"]],  0),
##                       method=rep("replace",3))
##   d_CLr = cbind(times=c(0,1), CLr=c(0, 0))

##   # Get parameters and initial values of state variables.
##   parms = initParms(newParms = c(p[c("CL","Vd","F_abs","P_milk","r_f_m","t_con",
##                                      "food_dose","pow_dose","pow_milk")],"n_i"=1))
##   Y0 = initStates(parms)

##   Y0["d_m"] = p[["d_m"]] # Dose (mg/kg/d)

##   Y0["A_mf"] = p[["A_m_init"]]

##   # Run simulation for times leading up to birth (if preg is TRUE) or all
##   # times (if preg is FALSE).
##   if (preg) {
##     times1 = times[times <= p[["t_gest"]]]
##   } else {
##     times1 = times
##   }

##   # Run model simulation to birth
##   out1 = run_model(mName, times, Y0, parms, rtol=1e-6, atol=1e-6,
##                    list(M_m_in, M_i_in, R_milk_in, r_m_mf_in, d_i_in),
##                    list(method="linear", rule=2, ties="ordered"),
##                    list(data=df_CLr))

##   # If this simulation does not incorporate any lactation/nursing period, stop
##   # here.
##   if (preg) {
##     if (max(times) <= p[["t_gest"]]) {
##       return(out1)
##     }
##   } else {
##     return(out1)
##   }

##   # Otherwise, continue with a simulation of the lactation/nursing period...

##   # Maternal mass interpolation points after moment of birth.
##   M_m_df = M_m_df3
##   t_m_vec = (M_m_df[["Age (y)"]] - 24.25) * 365
##   M_m_vec = M_m_df[["Mass (kg)"]]
##   M_m_in = cbind(times=t_m_vec, mass=M_m_vec)

##   # Infant mass interpolation points after moment of birth.
##   # Shift ages so that they reflect elapsed time since conception.
##   M_i_df = rbind(M_i_df1, M_i_df2)
##   t_i_vec = M_i_df[["Age (y)"]] * 365 + p[["t_gest"]]
##   m_i_vec = M_i_df[["Mass (kg)"]]
##   M_i_in = cbind(time=t_i_vec, mass=m_i_vec)

##   # Human breast milk ingestion rate interpolation points.
##   t_milk_0 = p[["t_gest"]]
##   t_milk_1e = t_milk_0 + 30
##   t_milk_2e = t_milk_0 + 90
##   t_milk_3e = t_milk_0 + 180
##   t_milk_4e = t_milk_0 + 365
##   t_milk_5e = t_milk_4e + 1e-6
##   r_milk_0 = 0.477
##   r_milk_1e = r_milk_0 +  2 * (0.510 - r_milk_0)
##   r_milk_2e = r_milk_1e + 2 * (0.690 - r_milk_1e)
##   r_milk_3e = r_milk_2e + 2 * (0.770 - r_milk_2e)
##   r_milk_4e = r_milk_3e + 2 * (0.620 - r_milk_3e)
##   r_milk_5e = 0.0
##   t_milk_vec = c(t_milk_0, t_milk_1e, t_milk_2e, t_milk_3e, t_milk_4e,
##                  t_milk_5e)
##   r_milk_vec = c(r_milk_0, r_milk_1e, r_milk_2e, r_milk_3e, r_milk_4e,
##                  r_milk_5e)
##   R_milk_in = cbind(times=t_milk_vec, rate=r_milk_vec)

##   # Apportion the amount of substance in the mother at the moment of birth to
##   # the mother and the infant(s) according to body mass. Use these values to
##   # set the initial conditions for the post-birth simulation.
##   M_m_birth = M_m_df3[[1, "Mass (kg)"]]
##   M_i_birth = m_i_vec[1]
##   C_m_birth = out1[[nrow(out1), "C_m"]]
##   Y0["A_m"] = C_m_birth * M_m_birth
##   Y0["A_i"] = C_m_birth * M_i_birth
##   Y0["AUC_m"] = out1[[nrow(out1), "AUC_m"]]
##   Y0["AUC_i"] = out1[[nrow(out1), "AUC_i"]]
##   Y0["d_m"] = out1[[nrow(out1), "d_m"]]
##   Y0["d_i"] = out1[[nrow(out1), "d_i"]]
##   Y0["T_in"] = out1[[nrow(out1), "T_in"]]
##   Y0["T_out"] = out1[[nrow(out1), "T_out"]] +
##     C_m_birth * (out1[[nrow(out1), "M_m"]] - M_m_birth - M_i_birth)

##   # Create an "event" so that direct dosing of the infant (filial) human
##   # begins at weaning and continues through the rest of the simulation.
##   df_dose = data.frame(var=c("d_i"),
##                        time=c(p[["t_gest"]] + p[["t_lact"]]),
##                        value=c(p[["d_m"]]),
##                        method=c("replace"))

##   # Run simulation for times after birth.
##   times2 = times[times >= p["t_gest"]]
##   out2 = run_model("lipophilic_tk", times2, Y0, parms,
##                    list(M_m_in, M_i_in, R_milk_in),
##                    list(method="linear", rule=2, ties="ordered"),
##                    list(data=df_dose))

##   # Assemble output data frame.
##   out = rbind(out1[1:nrow(out1)-1, ], out2)

##   # Return output data frame.
##   return(out)
## }

## # Legacy function
## init_rat_p <- function(new_p = NULL) {
##   # Default parameter values for rat simulations.
##   # Default parameter values for human simulations.
##   p = c(t_con = 90,  # Day of conception (from start of experiment).
##         t_gest = 22, # Duration (d) of pregnancy/gestation.
##         t_lact = 21,  # Duration (d) of lactation/nursing.
##         n_f = 10,     # Effective number of fetuses.
##         n_i = 10,     # Number of infants.
##         CL = 6.1,      # Clearance (mL/kg/d) of substance in this species.
##         Vd = 453,       # Volume of distribution (mL/kg) of substance in this species.
##         P_milk = 0.03,  # Milk/maternal blood concentration partition coefficient.
##         r_f_m = 1.0,    # Ratio of fetal and maternal concentrations.
##         F_abs = 1.0,    # Fraction of nominal dose that is absorbed.
##         food_dose = 0,  # A boolean flag. Dose is delivered in food if nonzero
##         d_m = 0,        # Dose rate (mg/kg/d or mg/kg) to mother.
##         pow_dose = 1.0, # When 1, d_m is mg/kg/d, when 0, d_m is mg/d
##         pow_milk = 0.0, # when 1, milk transfer rate is kg/kg/d, when 0, is kg/d
##         A_m_init = 0)   # Initial amount (mg) in mother.

##  # p <- c(t_gest = 22,        # Duration (d) of pregnancy/gestation.
##   #       t_lact = 21,        # Duration (d) of lactation/nursing.
##   #       M_m_1 = 0.25,       # Maternal mass 1 (kg).
##   #       M_m_2 = 0.273,      # Maternal mass 2 (kg).
##   #       t_m_1 = 1,          # Time (d) since conception for maternal mass 1.
##   #       t_m_2 = 44,         # Time (d) since conception for maternal mass 2.
##   #       M_i_1 = 0.0066,     # Infant mass 1 (kg).
##   #       M_i_2 = 0.014,      # Infant mass 2 (kg).
##   #       M_i_3 = 0.028,      # Infant mass 3 (kg).
##   #       M_i_4 = 0.25,       # Infant mass 4 (kg).
##   #       M_i_5 = NaN,        # Infant mass 5 (kg). Not used for rat.
##   #       M_i_6 = NaN,        # Infant mass 6 (kg). Not used for rat.
##   #       t_i_1 = 3,          # Time (d) since birth for infant mass 1.
##   #       t_i_2 = 10,         # Time (d) since birth for infant mass 2.
##   #       t_i_3 = 17,         # Time (d) since birth for infant mass 3.
##   #       t_i_4 = 75,         # Time (d) since birth for infant mass 4.
##   #       t_i_5 = NaN,        # Time (d) since birth for infant mass 5.
##   #       t_i_6 = NaN,        # Time (d) since birth for infant mass 6.
##   #       r_milk_0 = 0.001,   # Rate (kg/d) of milk ingestion at birth.
##   #       r_milk_1 = 0.003,   # Rate (kg/d) of milk ingestion during week 1.
##   #       r_milk_2 = 0.0054,  # Rate (kg/d) of milk ingestion during week 2.
##   #       r_milk_3 = 0.0059,  # Rate (kg/d) of milk ingestion during week 3.
##   #       r_milk_4 = NaN,     # Rate (kg/d) of milk ingestion during week 4.
##   #       t_m_start = 0,      # Time (d) since conception maternal dose starts.
##   #       t_m_end = 1,        # Time (d) since conception maternal dose ends.
##   #       A_m_init = 0,       # Initial amount (mg) in mother.
##   #       A_i_init = 0,       # Initial amount (mg) in infant.
##   #       d_i = 0,            # Dose rate (mg/kg/d or mg/kg) to infant.
##   #       t_i_start = 0,      # Time (d) since conception infant dose starts.
##   #       t_i_end = 1)        # Time (d) since conception infant dose ends.

##   # If "new_p" is not NULL, replace default parameter values with values from
##   # that vector.
##   if (!is.null(new_p)) {
##     if (!all(names(new_p) %in% c(names(p)))) {
##       stop("Illegal parameter name in 'new_p'.")
##     } else {
##       p[names(new_p)] <- new_p
##     }
##   }

##   return(p)
## }

## # Legacy function
## init_mouse_p <- function(new_p = NULL) {
##   # Default parameter values for mouse simulations.
##   p <- c(t_gest = 18,        # Duration (d) of pregnancy/gestation.
##          t_lact = 21,        # Duration (d) of lactation/nursing.
##          M_m_1 = 0.0255,     # Maternal mass 1 (kg).
##          M_m_2 = 0.0318,     # Maternal mass 2 (kg).
##          t_m_1 = 1,          # Time (d) since conception for maternal mass 1.
##          t_m_2 = 25,         # Time (d) since conception for maternal mass 2.
##          n_f = 6,            # Effective number of fetuses.
##          n_i = 6,            # Number of infants.
##          M_i_1 = 0.0014,     # Infant mass 1 (kg).
##          M_i_2 = 0.00608,    # Infant mass 2 (kg).
##          M_i_3 = 0.00885,    # Infant mass 3 (kg).
##          M_i_4 = 0.03,       # Infant mass 4 (kg).
##          M_i_5 = NaN,        # Infant mass 5 (kg). Not used for mouse.
##          M_i_6 = NaN,        # Infant mass 6 (kg). Not used for mouse.
##          t_i_1 = 1,          # Time (d) since birth for infant mass 1.
##          t_i_2 = 10,         # Time (d) since birth for infant mass 2.
##          t_i_3 = 18,         # Time (d) since birth for infant mass 3.
##          t_i_4 = 45,         # Time (d) since birth for infant mass 4.
##          t_i_5 = NaN,        # Time (d) since birth for infant mass 5.
##          t_i_6 = NaN,        # Time (d) since birth for infant mass 6.
##          pBW_milk = NaN,     # Per body mass rate of milk ingestion (kg/kg/d).
##          r_milk_0 = 0.0001,  # Rate (kg/d) of milk ingestion at birth.
##          r_milk_1 = 0.0003,  # Rate (kg/d) of milk ingestion during week 1.
##          r_milk_2 = 0.00054, # Rate (kg/d) of milk ingestion during week 2.
##          r_milk_3 = 0.00059, # Rate (kg/d) of milk ingestion during week 3.
##          r_milk_4 = NaN,     # Rate (kg/d) of milk ingestion during week 4.
##          CL = 6.1,           # Clearance (mL/kg/d) of substance in this species.
##          Vd = 453,           # Volume of distribution (mL/kg) of substance in this species.
##          P_milk = 1.0,       # Milk/maternal blood concentration partition coefficient.
##          r_f_m = 1.0,        # Ratio of fetal and maternal concentrations.
##          F_abs = 1.0,        # Fraction of nominal dose that is absorbed.
##          food_dose = 0,      # A boolean flag. Dose is delivered in food if
##          #   nonzero.
##          d_m = 0,            # Dose rate (mg/kg/d or mg/kg) to mother.
##          t_m_start = 0,      # Time (d) since conception maternal dose starts.
##          t_m_end = 1,        # Time (d) since conception maternal dose ends.
##          A_m_init = 0,       # Initial amount (mg) in mother.
##          A_i_init = 0,       # Initial amoung (mg) in infant.
##          d_i = 0,            # Dose rate (mg/kg/d or mg/kg) to infant.
##          t_i_start = 0,      # Time (d) since conception infant dose starts.
##          t_i_end = 1)        # Time (d) since conception infant dose ends.

##   # If "new_p" is not NULL, replace default parameter values with values from
##   # that vector.
##   if (!is.null(new_p)) {
##     if (!all(names(new_p) %in% c(names(p)))) {
##       stop("Illegal parameter name in 'new_p'.")
##     } else {
##       p[names(new_p)] <- new_p
##     }
##   }

##   return(p)

##   # If "new_p" is not NULL, replace default parameter values with values from
##   # that vector.
##   if (!is.null(new_p)) {
##     if (!all(names(new_p) %in% c(names(p)))) {
##       stop("Illegal parameter name in 'new_p'.")
##     } else {
##       p[names(new_p)] <- new_p
##     }
##   }

##   return(p)
## }

## # Legacy function
## init_monkey_p <- function(new_p = NULL) {
##   # Default parameter values for human simulations.
##   p = c(t_gest = 168,       # Duration (d) of pregnancy/gestation.
##         t_lact = 154,       # Duration (d) of lactation/nursing.
##         M_m_1 = 7.1,        # Maternal mass 1 (kg).
##         M_m_2 = 7.1,        # Maternal mass 2 (kg).
##         t_m_1 = 1,          # Time (d) since conception for maternal mass 1.
##         t_m_2 = 2,          # Time (d) since conception for maternal mass 2.
##         n_f = 1,            # Effective number of fetuses.
##         n_i = 1,            # Number of infants.
##         M_i_1 = 0.4775,     # Infant mass 1 (kg).
##         M_i_2 = 0.9375,     # Infant mass 2 (kg).
##         M_i_3 = 1.4275,     # Infant mass 3 (kg).
##         M_i_4 = 2.1925,     # Infant mass 4 (kg).
##         M_i_5 = 3.4275,     # Infant mass 5 (kg).
##         M_i_6 = 5.270,      # Infant mass 6 (kg).
##         t_i_1 = 1,          # Time (d) since birth for infant mass 1.
##         t_i_2 = 91,         # Time (d) since birth for infant mass 2.
##         t_i_3 = 182,        # Time (d) since birth for infant mass 3.
##         t_i_4 = 365,        # Time (d) since birth for infant mass 4.
##         t_i_5 = 730,        # Time (d) since birth for infant mass 5.
##         t_i_6 = 1095,       # Time (d) since birth for infant mass 6.
##         pBW_milk = 0.1,     # Per body mass rate of milk ingestion (kg/kg/d).
##         CL = 6.1,           # Clearance (mL/kg/d) of substance in this species.
##         Vd = 453,           # Volume of distribution (mL/kg) of substance in this species.
##         P_milk = 1.0,       # Milk/maternal blood concentration partition coefficient.
##         r_f_m = 1.0,        # Ratio of fetal and maternal concentrations.
##         F_abs = 1.0,        # Fraction of nominal dose that is absorbed.
##         food_dose = 0,      # A boolean flag. Dose is delivered in food if
##         #   nonzero.
##         d_m = 0,            # Dose rate (mg/kg/d or mg/kg) to mother.
##         t_m_start = 0,      # Time (d) since conception maternal dose starts.
##         t_m_end = 1,        # Time (d) since conception maternal dose ends.
##         A_m_init = 0,       # Initial amount (mg) in mother.
##         A_i_init = 0,       # Initial amoung (mg) in infant.
##         d_i = 0,            # Dose rate (mg/kg/d or mg/kg) to infant.
##         t_i_start = 0,      # Time (d) since conception infant dose starts.
##         t_i_end = 1)        # Time (d) since conception infant dose ends.

##   # If "new_p" is not NULL, replace default parameter values with values from
##   # that vector.
##   if (!is.null(new_p)) {
##     if (!all(names(new_p) %in% c(names(p)))) {
##       stop("Illegal parameter name in 'new_p'.")
##     } else {
##       p[names(new_p)] <- new_p
##     }
##   }

##   return(p)
## }

## #Legacy function
## rat_sim <- function(times=NULL, p=NULL) {
##   # If p is NULL, or if not all parameters are provided in p, use default
##   # values for the parameters not provided.
##   p = init_rat_p(p)

##   # If t is NULL, use default times for simulation.
##   if (is.null(times)) {
##     t_start = min(0, p[["t_m_start"]])
##     t_end = max(p[["t_gest"]] + p[["t_lact"]], p[["t_i_end"]])
##     times = seq(from=t_start, to=t_end, by=0.1)
##   }

##   # Print an error message if there is a problem with the input parameters
##   # relevant to dam mass and rat infant mass.
##   try(if (p[["t_m_1"]] >= p[["t_m_2"]])
##     stop("Time t_m_1 must be before time t_m_2."))
##   try(if (p[["M_m_1"]] <= 0) stop ("M_m_1 must be positive."))
##   try(if (p[["M_m_2"]] <= 0) stop ("M_m_2 must be positive."))
##   try(if (p[["t_i_1"]] > p[["t_i_2"]])
##     stop("Time t_i_1 must be before time t_i_2."))
##   try(if (p[["t_i_2"]] > p[["t_i_3"]])
##     stop("Time t_i_2 must be before time t_i_3."))
##   try(if (p[["t_i_3"]] > p[["t_i_4"]])
##     stop("Time t_i_3 must be before time t_i_4."))
##   try(if (p[["M_i_1"]] <= 0) stop ("M_i_1 must be positive."))
##   try(if (p[["M_i_2"]] <= 0) stop ("M_i_2 must be positive."))
##   try(if (p[["M_i_3"]] <= 0) stop ("M_i_3 must be positive."))
##   try(if (p[["M_i_4"]] <= 0) stop ("M_i_4 must be positive."))
##   if (is.finite(p[["M_i_5"]])) {
##     try(if (p[["t_i_4"]] > p[["t_i_5"]])
##       stop("Time t_i_4 must be before time t_i_5."))
##     try(if (p[["t_i_5"]] > p[["t_i_6"]])
##       stop("Time t_i_5 must be before time t_i_6."))
##     try(if (p[["M_i_5"]] <= 0) stop ("M_i_5 must be positive."))
##     try(if (p[["M_i_6"]] <= 0) stop ("M_i_6 must be positive."))
##   }

##   # Calculate total mass of pups at birth.
##   M_pups = p[["n_f"]] * p[["M_i_1"]]

##   # If the simulation starts before birth, run the pregnancy simulation.
##   if (times[1] <= p[["t_gest"]]) {
##     # Generate points to interpolate for calculations of dam mass up to moment of
##     # birth.
##     if (p[["t_m_2"]] <= p[["t_gest"]]) { # Case 1: t_m_2 <= t_gest.
##       if (p[["M_m_2"]] < (p[["M_m_1"]] + M_pups)) {
##         M_m_birth = p[["M_m_1"]] + M_pups  # Mass just before birth.
##       } else {
##         M_m_birth = p[["M_m_2"]]           # Mass just before birth.
##       }
##       t_m_vec = c(p[["t_m_1"]], p[["t_m_2"]], p[["t_gest"]])
##       M_m_vec = c(p[["M_m_1"]], p[["M_m_2"]], M_m_birth)
##     } else if (p[["t_m_1"]] <= p[["t_gest"]]) { # Case 2: t_m_1 <= t_gest < t_m_2.
##       M_m_birth = p[["M_m_2"]] + M_pups  # Mass just before birth.
##       t_m_vec = c(p[["t_m_1"]], p[["t_gest"]])
##       M_m_vec = c(p[["M_m_1"]], M_m_birth)
##     } else {
##       stop("Time t_m_1 must be before time t_gest.")
##     }

##     # Rat dam mass interpolation points up to moment of birth.
##     M_mf_in = cbind(times=t_m_vec, mass=M_m_vec)

##     # Rat infant mass interpolation points up to moment of birth.
##     M_i_in = cbind(times=c(0, 1), mass=c(0, 0))

##     # Rat breast milk ingestion rate interpolation points.
##     R_milk_in = cbind(times=c(0, 1), rate=c(0, 0))

##     # Dosing events. Turn dosing on and off at specified times.
##     df_dose = data.frame(var=c("d_m", "d_m", "d_i", "d_i"),
##                          time=c(p[["t_m_start"]], p[["t_m_end"]],
##                                 p[["t_i_start"]], p[["t_i_end"]]),
##                          value=c(p[["d_m"]], 0, p[["d_i"]], 0),
##                          method=rep("replace", 4))
##     df_dose = df_dose[order(df_dose["time"]), ]
##     rownames(df_dose) = 1:nrow(df_dose)

##     # Mass of mother and pups just after birth.
##     M_m_b = M_m_birth - M_pups

##     # Ratio of concentrations in mother and mother + fetus(es).
##     r_m_mf_0 = 1.0
##     r_m_mf_b = (M_m_b + M_pups) / (M_m_b + p[["r_f_m"]] * M_pups)
##     r_m_mf_in = cbind(times=c(0, p[["t_gest"]]), ratio=c(r_m_mf_0, r_m_mf_b))

##     # Get parameters and initial values of state variables.
##     parms = initParms()
##     parms["CL"] = p[["CL"]]
##     parms["Vd"] = p[["Vd"]]
##     parms["r_f_m"] = p[["r_f_m"]]
##     parms["F_abs"] = p[["F_abs"]]
##     parms["n_i"] = p[["n_i"]]
##     parms["food_dose"] = p[["food_dose"]]
##     parms = initParms(parms)   # Update calculated parameters.
##     Y0 = initStates(parms)
##     Y0["A_mf"] = p[["A_m_init"]]
##     Y0["A_i"] = 0
##     Y0["T_in"] = Y0[["A_mf"]]

##     # Run simulation for times leading up to birth.
##     times1 = times[times <= p[["t_gest"]]]
##     out1 = run_model("lipophilic_tk", times1, Y0, parms,
##                      list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in),
##                      list(method="linear", rule=2, ties="ordered"),
##                      list(data=df_dose))

##     # If this simulation does not incorporate any lactation/nursing period,
##     # stop here.
##     if (max(times) <= p[["t_gest"]]) {
##       return(out1)
##     }
##   }

##   # Continue with a simulation of the lactation/nursing period if necessary...

##   # If there are results from a pregnancy simulation (in "out1"), use those to
##   # set initial conditions.
##   if (times[1] <= p[["t_gest"]]) {
##     # Apportion the amount of substance in the mother at the moment of birth to
##     # the mother and the infant(s) according to body mass. Use these values to
##     # set the initial conditions for the post-birth simulation.
##     M_m_birth = out1[[nrow(out1), "M_mf"]] - M_pups
##     M_i_birth = p[["M_i_1"]]
##     C_m_birth = out1[[nrow(out1), "C_m"]]
##     C_i_birth = out1[[nrow(out1), "C_i"]]
##     Y0["A_mf"] = C_m_birth * M_m_birth
##     Y0["A_i"] = C_i_birth * M_i_birth
##     Y0["AUC_m"] = out1[[nrow(out1), "AUC_m"]]
##     Y0["AUC_i"] = out1[[nrow(out1), "AUC_i"]]
##     Y0["d_m"] = out1[[nrow(out1), "d_m"]]
##     Y0["d_i"] = out1[[nrow(out1), "d_i"]]
##     Y0["T_in"] = out1[[nrow(out1), "T_in"]]
##     Y0["T_out"] = out1[[nrow(out1), "T_out"]] +
##       C_i_birth * M_i_birth * (p[["n_f"]] - p[["n_i"]])
##   } else {
##     # Get parameters and initial values of state variables.
##     parms = initParms()
##     parms["CL"] = p[["CL"]]
##     parms["Vd"] = p[["Vd"]]
##     parms["r_f_m"] = p[["r_f_m"]]
##     parms["F_abs"] = p[["F_abs"]]
##     parms["n_i"] = p[["n_i"]]
##     parms["food_dose"] = p[["food_dose"]]
##     parms = initParms(parms)   # Update calculated parameters.
##     Y0 = initStates(parms)
##     Y0["A_mf"] = p[["A_m_init"]]
##     Y0["A_i"] = p[["A_i_init"]]
##     Y0["T_in"] = Y0[["A_mf"]] + Y0[["A_i"]]

##     # Set dosing data frame.
##     # Dosing events. Turn dosing on and off at specified times.
##     df_dose = data.frame(var=c("d_m", "d_m", "d_i", "d_i"),
##                          time=c(p[["t_m_start"]], p[["t_m_end"]],
##                                 p[["t_i_start"]], p[["t_i_end"]]),
##                          value=c(p[["d_m"]], 0, p[["d_i"]], 0),
##                          method=rep("replace", 4))
##     df_dose = df_dose[order(df_dose["time"]), ]
##     rownames(df_dose) = 1:nrow(df_dose)

##     # Set r_m_mf_in parameter (which is not used).
##     r_m_mf_in = cbind(times=c(0, 1), ratio=c(0, 0))
##   }


##   # Generate points to interpolate for calculations of dam mass after moment of
##   # birth.
##   if (p[["t_m_2"]] <= p[["t_gest"]]) { # Case 1: t_m_2 <= t_gest.
##     if (p[["M_m_2"]] < (p[["M_m_1"]] + M_pups)) {
##       M_m_birth = p[["M_m_1"]]           # Mass just after birth.
##     } else {
##       M_m_birth = p[["M_m_2"]] - M_pups  # Mass just after birth.
##     }
##     t_m_vec = c(p[["t_gest"]], p[["t_gest"]] + 1)
##     M_m_vec = c(M_m_birth, M_m_birth)
##   } else if (p[["t_m_1"]] <= p[["t_gest"]]) { # Case 2: t_m_1 <= t_gest < t_m_2.
##     M_m_birth = p[["M_m_2"]]             # Mass just after birth.
##     t_m_vec = c(p[["t_gest"]], p[["t_gest"]] + 1)
##     M_m_vec = c(M_m_birth, M_m_birth)
##   } else {
##     stop("Time t_m_1 must be before time t_gest.")
##   }

##   # Rat dam mass interpolation points after moment of birth.
##   M_mf_in = cbind(times=t_m_vec, mass=M_m_vec)

##   # Rat infant mass interpolation points.
##   if (is.nan(p[["M_i_5"]])) {
##     t_i_vec = p[["t_gest"]] + c(0, p[["t_i_1"]], p[["t_i_2"]], p[["t_i_3"]],
##                                 p[["t_i_4"]])
##     M_i_vec = c(p[["M_i_1"]], p[["M_i_1"]], p[["M_i_2"]], p[["M_i_3"]],
##                 p[["M_i_4"]])
##   } else {
##     t_i_vec = p[["t_gest"]] + c(0, p[["t_i_1"]], p[["t_i_2"]], p[["t_i_3"]],
##                                 p[["t_i_4"]], p[["t_i_5"]], p[["t_i_6"]])
##     M_i_vec = c(p[["M_i_1"]], p[["M_i_1"]], p[["M_i_2"]], p[["M_i_3"]],
##                 p[["M_i_4"]], p[["M_i_5"]], p[["M_i_6"]])
##   }
##   M_i_in = cbind(times=t_i_vec, mass=M_i_vec)

##   # Rat breast milk ingestion rate interpolation points.
##   if (is.nan(p[["pBW_milk"]])) {
##     t_milk_0 = p[["t_gest"]]
##     t_milk_1e = t_milk_0 + 7
##     t_milk_2e = t_milk_0 + 14
##     t_milk_3e = t_milk_0 + 21
##     r_milk_0 = p[["r_milk_0"]]
##     r_milk_1e = r_milk_0 + 7.0 * (p[["r_milk_1"]] - p[["r_milk_0"]]) / 3.5
##     r_milk_2e = r_milk_1e + 7.0 * (p[["r_milk_2"]] - r_milk_1e) / 3.5
##     r_milk_3e = r_milk_2e + 7.0 * (p[["r_milk_3"]] - r_milk_2e) / 3.5
##     if (is.nan(p[["r_milk_4"]])) {
##       t_milk_4e = t_milk_0 + p[["t_lact"]]
##       r_milk_4e = 0.0
##       t_milk_vec = c(t_milk_0, t_milk_1e, t_milk_2e, t_milk_3e, t_milk_4e)
##       r_milk_vec = c(r_milk_0, r_milk_1e, r_milk_2e, r_milk_3e, r_milk_4e)
##     } else {
##       t_milk_4e = t_milk_0 + 28
##       r_milk_4e = r_milk_3e + 7.0 * (p[["r_milk_4"]] - r_milk_3e) / 3.5
##       t_milk_5e = t_milk_0 + p[["t_lact"]]
##       r_milk_5e = 0.0
##       t_milk_vec = c(t_milk_0, t_milk_1e, t_milk_2e, t_milk_3e, t_milk_4e,
##                      t_milk_5e)
##       r_milk_vec = c(r_milk_0, r_milk_1e, r_milk_2e, r_milk_3e, r_milk_4e,
##                      r_milk_5e)
##     }
##   } else { # Milk consumption will be based on infant mass.
##     t_wean = p[["t_gest"]] + p[["t_lact"]]
##     M_wean = approx(t_i_vec, M_i_vec, xout=t_wean)$y
##     t_milk_vec = c(t_i_vec, t_wean)
##     r_milk_vec = p[["pBW_milk"]] * c(M_i_vec, M_wean)
##     r_milk_vec = r_milk_vec[t_milk_vec <= t_wean]
##     t_milk_vec = t_milk_vec[t_milk_vec <= t_wean]
##     t_milk_vec = c(t_milk_vec, t_wean)
##     r_milk_vec = c(r_milk_vec, 0.0)
##   }
##   R_milk_in = cbind(times=t_milk_vec, rate=r_milk_vec)

##   # Run simulation for times after birth.
##   times2 = times[times >= p["t_gest"]]
##   out2 = run_model("lipophilic_tk", times2, Y0, parms,
##                    list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in),
##                    list(method="linear", rule=2, ties="ordered"),
##                    list(data=df_dose))

##   # Assemble output data frame. Combine pre-birth and post-birth simulations if
##   # necessary.
##   if (times[1] <= p[["t_gest"]]) {
##     out = rbind(out1[1:nrow(out1)-1, ], out2)
##   } else {
##     out = out2
##   }

##   # Return output data frame.
##   return(out)
## }

## # Legacy function
## mouse_sim <- function(times=NULL, p=NULL) {
##   # If p is NULL, or if not all parameters are provided in p, use default
##   # values for the parameters not provided.
##   p = init_mouse_p(p)

##   # Run the simulation as for rat and return the result.
##   out = rat_sim(times=times, p=p)
## }

## # Legacy function
## mink_sim <- function(times=NULL, p=NULL) {
##   # If p is NULL, or if not all parameters are provided in p, use default
##   # values for the parameters not provided.
##   p = init_mink_p(p)

##   # Run the simulation as for rat and return the result.
##   out = rat_sim(times=times, p=p)
## }

## # Legacy function
## monkey_sim <- function(times=NULL, p=NULL) {
##   # If p is NULL, or if not all parameters are provided in p, use default
##   # values for the parameters not provided.
##   p = init_monkey_p(p)

##   # Run the simulation as for rat and return the result.
##   out = rat_sim(times=times, p=p)
## }

## # Legacy function
## rat_demo1 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(A_m_init = 100)     # Initial amount (mg) in mother.

##   # Run simulation.
##   out <- rat_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## rat_demo2 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(d_m = 6,            # Dose rate (mg/kg/d) to mother.
##         t_m_start = 6,      # Time (d) since conception maternal dose starts.
##         t_m_end = 28)       # Time (d) since conception maternal dose ends.

##   # Run simulation
##   out <- rat_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## rat_demo3 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Adult (non-developmental study) parameters.
##   r_dose = 6                # Dose rate (mg/kg/d) to adult animal.
##   t_dose = 22               # Duration (d) of dose.
##   t_exp = 28                # Duration (d) of experiment
##   F_abs = 1.0               # Fraction of dose absorbed.
##   r_f_m = 1.0               # Ratio of concentrations in fetus and mother.
##   # Setting this value to 1.0 will cause r_m_mf to be
##   # calculated as 1.0 in the function "rat_sim".

##   # Define simulation parameters.
##   p = c(d_m = r_dose,               # Dose rate (mg/kg/d) to "mother".
##         t_m_start = -t_exp,         # Time (d) since conception "maternal" dose
##         #  starts.
##         t_m_end = -t_exp + t_dose,  # Time (d) since conception "maternal" dose
##         #  ends.
##         F_abs = F_abs,
##         r_f_m = r_f_m)

##   # Define simulation times.
##   times = seq(from=-t_exp, to=0, by=0.1)

##   # Run simulation.
##   out <- rat_sim(times=times, p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1] + t_exp, out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }

##   # # Calculate average adult body burden over entire simulation.
##   # p = init_rat_p(p)
##   # DM5 = comp_dose_metric(out[ , 1] + t_exp, out[ , "C_m"],
##   # out[ , "AUC_m"], p[["t_gest"]], p[["t_lact"]], 5)

##   # # Print dose metric 5.
##   # cat("Average adult body burden =", DM5, "mg/kg")
## }

## # Legacy function
## rat_demo4 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(d_i = 6,            # Dose rate (mg/kg/d) to infant.
##         t_i_start = 25,     # Time (d) since conception infant dose starts.
##         t_i_end = 30)       # Time (d) since conception infant dose ends.

##   # Run simulation.
##   out <- rat_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }

##   # # Calculate dose metric 2 (average concentration over entire period of interest).
##   # p = init_rat_p(p)
##   # DM2 = comp_dose_metric(out[ , "time"], out[ , "C_i"], out[ , "AUC_i"],
##   # p["t_gest"], p["t_lact"], 2)

##   # # Print dose metric 2.
##   # cat("DM2 =", DM2 , "mg/kg")
## }

## # Legacy function
## rat_demo5 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(d_m = 6,            # Dose rate (mg/kg/d) to mother.
##         t_m_start = -10,    # Time (d) since conception maternal dose starts.
##         t_m_end = 28)       # Time (d) since conception maternal dose ends.

##   # Run simulation
##   out <- rat_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## mouse_demo1 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(A_m_init = 100)     # Initial amount (mg) in mother.

##   # Run simulation.
##   out <- mouse_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## mouse_demo2 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(d_m = 6,            # Dose rate (mg/kg/d) to mother.
##         t_m_start = 6,      # Time (d) since conception maternal dose starts.
##         t_m_end = 28)       # Time (d) since conception maternal dose ends.

##   # Run simulation.
##   out <- mouse_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## mink_demo1 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(A_m_init = 100)     # Initial amount (mg) in mother.

##   # Run simulation.
##   out <- mink_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## mink_demo2 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(d_m = 6,            # Dose rate (mg/kg/d) to mother.
##         t_m_start = 6,      # Time (d) since conception maternal dose starts.
##         t_m_end = 70)       # Time (d) since conception maternal dose ends.

##   # Run simulation.
##   out <- mink_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## monkey_demo1 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(A_m_init = 100)     # Initial amount (mg) in mother.

##   # Run simulation.
##   out <- monkey_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## monkey_demo2 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(d_m = 6,            # Dose rate (mg/kg/d) to mother.
##         t_m_start = 6,      # Time (d) since conception maternal dose starts.
##         t_m_end = 200)      # Time (d) since conception maternal dose ends.

##   # Run simulation.
##   out <- monkey_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## human_sim_old <- function(times=NULL, p=NULL, M_m_df1=NULL, M_m_df2=NULL,
##                           M_m_df3=NULL) {
##   # If p is NULL, or if not all parameters are provided in p, use default
##   # default values for the parameters not provided.
##   p = init_human_p(p)

##   # If t is NULL, use default times for simulation.
##   if (is.null(times)) {
##     t_start = -24.25 * 365
##     t_end = p[["t_gest"]] + p[["t_lact"]]
##     times = seq(from=t_start, to=t_end, by=0.25)
##   }

##   # If not provided, obtain time course data for mass of human mother.
##   if (is.null(M_m_df1)) {
##     M_m_df1 = read.csv("human_female_mass_0_to_20.csv", header=FALSE,
##                        col.names=c("Age (y)", "Mass (kg)"), check.names=FALSE)
##   }
##   if (is.null(M_m_df2)) {
##     M_m_df2 = read.csv("human_female_mass_20_to_25_pregnancy.csv",
##                        header=FALSE, col.names=c("Age (y)", "Mass (kg)"),
##                        check.names=FALSE)
##   }
##   if (is.null(M_m_df3)) {
##     M_m_df3 = read.csv("human_female_mass_25_to_26_postnatal.csv",
##                        header=FALSE, col.names=c("Age (y)", "Mass (kg)"),
##                        check.names=FALSE)
##   }

##   # Maternal mass interpolation points up to moment of birth.
##   M_m_df = rbind(M_m_df1, M_m_df2)
##   t_m_vec = (M_m_df[["Age (y)"]] - 24.25) * 365
##   M_m_vec = M_m_df[["Mass (kg)"]]
##   M_m_in = cbind(times=t_m_vec, mass=M_m_vec)

##   # Infant mass interpolation points up to moment of birth.
##   M_i_in = cbind(times=c(0, 1), mass=c(0, 0))

##   # Human breast milk ingestion rate interpolation points.
##   R_milk_in = cbind(times=c(0, 1), rate=c(0, 0))

##   # Get parameters and initial values of state variables.
##   parms = initParms()
##   parms["CL"] = p[["CL"]]
##   parms["Vd"] = p[["Vd"]]
##   parms["F_abs"] = p[["F_abs"]]
##   parms["n_i"] = 1
##   parms["food_dose"] = p[["food_dose"]]
##   parms = initParms(parms)   # Update calculated parameters.
##   Y0 = initStates(parms)
##   Y0["A_m"] = p[["A_m_init"]]
##   Y0["T_in"] = Y0[["A_m"]]
##   Y0["d_m"] = p[["d_m"]]

##   # Run simulation for times leading up to birth.
##   times1 = times[times <= p[["t_gest"]]]
##   out1 = run_model("lipophilic_tk", times1, Y0, parms,
##                    list(M_m_in, M_i_in, R_milk_in),
##                    list(method="linear", rule=2, ties="ordered"))

##   # If this simulation does not incorporate any lactation/nursing period, stop
##   # here.
##   if (max(times) <= p[["t_gest"]]) {
##     return(out1)
##   }

##   # Otherwise, continue with a simulation of the lactation/nursing period...

##   # Maternal mass interpolation points after moment of birth.
##   M_m_df = M_m_df3
##   t_m_vec = (M_m_df[["Age (y)"]] - 24.25) * 365
##   M_m_vec = M_m_df[["Mass (kg)"]]
##   M_m_in = cbind(times=t_m_vec, mass=M_m_vec)

##   # Infant mass interpolation points after moment of birth.
##   m_i_1000 = 7.4 + (1000 - 272.5) * (9.2 - 7.4) / (272.5 - 135)
##   t_i_vec = c(0, 15, 60, 135, 272.5, 1000) + p[["t_gest"]]
##   m_i_vec = c(3.36, 4.8, 5.9, 7.4, 9.2, m_i_1000)

##   # Human breast milk ingestion rate interpolation points.
##   t_milk_0 = p[["t_gest"]]
##   t_milk_1e = t_milk_0 + 30
##   t_milk_2e = t_milk_0 + 90
##   t_milk_3e = t_milk_0 + 180
##   t_milk_4e = t_milk_0 + 365
##   t_milk_5e = t_milk_4e + 1e-6
##   r_milk_0 = 0.477
##   r_milk_1e = r_milk_0 +  2 * (0.510 - r_milk_0)
##   r_milk_2e = r_milk_1e + 2 * (0.690 - r_milk_1e)
##   r_milk_3e = r_milk_2e + 2 * (0.770 - r_milk_2e)
##   r_milk_4e = r_milk_3e + 2 * (0.620 - r_milk_3e)
##   r_milk_5e = 0.0
##   t_milk_vec = c(t_milk_0, t_milk_1e, t_milk_2e, t_milk_3e, t_milk_4e,
##                  t_milk_5e)
##   r_milk_vec = c(r_milk_0, r_milk_1e, r_milk_2e, r_milk_3e, r_milk_4e,
##                  r_milk_5e)
##   R_milk_in = cbind(times=t_milk_vec, rate=r_milk_vec)

##   # Apportion the amount of substance in the mother at the moment of birth to
##   # the mother and the infant(s) according to body mass. Use these values to
##   # set the initial conditions for the post-birth simulation.
##   M_m_birth = M_m_df3[[1, "Mass (kg)"]]
##   M_i_birth = m_i_vec[1]
##   C_m_birth = out1[[nrow(out1), "C_m"]]
##   Y0["A_m"] = C_m_birth * M_m_birth
##   Y0["A_i"] = C_m_birth * M_i_birth
##   Y0["AUC_m"] = out1[[nrow(out1), "AUC_m"]]
##   Y0["AUC_i"] = out1[[nrow(out1), "AUC_i"]]
##   Y0["d_m"] = out1[[nrow(out1), "d_m"]]
##   Y0["d_i"] = out1[[nrow(out1), "d_i"]]
##   Y0["T_in"] = out1[[nrow(out1), "T_in"]]
##   Y0["T_out"] = out1[[nrow(out1), "T_out"]] +
##     C_m_birth * (out1[[nrow(out1), "M_m"]] - M_m_birth - M_i_birth)


##   # Run simulation for times after birth.
##   times2 = times[times >= p["t_gest"]]
##   out2 = run_model("lipophilic_tk", times2, Y0, parms,
##                    list(M_m_in, M_i_in, R_milk_in),
##                    list(method="linear", rule=2, ties="ordered"))

##   # Assemble output data frame.
##   out = rbind(out1[1:nrow(out1)-1, ], out2)

##   # Return output data frame.
##   return(out)
## }

## # Legacy function
## human_demo1 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(A_m_init = 100)       # Initial amount (mg) in mother.

##   # Run simulation.
##   t_start = -24.25 * 365
##   t_end = 78.75 * 365          # Simulate to 78.75 years post-conception.
##   times = seq(from=t_start, to=t_end, by=0.25)
##   out <- human_sim(times=times, p=p)
##   # out <- human_sim(p=p)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1] / 365 + 24.25, out[ , idx], type='l',
##          xlab="Age of Mother (y)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## human_demo2 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(d_m = 0.0063)            # Dose rate (mg/kg/d) to mother.

##   # Run simulation.
##   t_start = -24.25 * 365
##   t_end = 30.0 * 365            # Simulate to 30 years post-conception.
##   times = seq(from=t_start, to=t_end, by=0.25)
##   out <- human_sim(times=times, p=p)
##   # out <- human_sim(p=p)

##   # Plot all output variables vs. time.
##   conception_idx = match(0, out[ , 1])
##   for (idx in 2:ncol(out)) {
##     plot(out[conception_idx:nrow(out), 1] / 365 + 24.25,
##          out[conception_idx:nrow(out), idx],
##          type='l', xlab="Age of Mother (y)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## human_demo3 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(A_m_init = 100)       # Initial amount (mg) in mother.

##   # Run simulation.
##   out <- human_sim(p=p, preg=FALSE)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1] / 365, out[ , idx], type='l',
##          xlab="Age of Mother (y)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## human_demo4 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Define simulation parameters.
##   p = c(d_m = 0.0063)            # Dose rate (mg/kg/d) to mother.

##   # Run simulation.
##   out <- human_sim(p=p, preg=FALSE)

##   # Plot all output variables vs. time.
##   for (idx in 2:ncol(out)) {
##     plot(out[ , 1] / 365, out[ , idx], type='l',
##          xlab="Age of Mother (y)",
##          ylab=colnames(out)[idx])
##   }
## }

## # Legacy function
## comp_dose_metric <- function(times, conc, AUC, t_gest, t_lact, DM_code) {
##   conception_idx = match(-t_gest, times)
##   birth_idx = match(0, times)
##   wean_idx = match(t_lact, times)

##   if (DM_code == 1) {         # Peak concentration during gestation/nursing.
##     out = (max(conc[conception_idx:wean_idx]))
##   } else if (DM_code == 2) {  # Average concentration during gestation/nursing.
##     out = (AUC[wean_idx] - AUC[conception_idx]) / (times[wean_idx]
##                                                    - times[conception_idx])
##   } else if (DM_code == 3) {  # Average concentration during gestation.
##     out = (AUC[birth_idx] - AUC[conception_idx]) / (times[birth_idx]
##                                                     - times[conception_idx])
##   } else if (DM_code == 4) {  # Average concentration during nursing.
##     out = (AUC[wean_idx] - AUC[birth_idx]) / (times[wean_idx]
##                                               - times[birth_idx])
##   } else if (DM_code == 5) {  # Average concentration during simulation.
##     out = (AUC[length(AUC)] - AUC[conception_idx]) / (times[length(times)]
##                                                       - times[conception_idx])
##   }
##   else {
##     stop("Dose metric code is not recognized.")
##   }
##   return(out)
## }

## # Legacy function
## cost_fun_old <- function(d_m, DM_code, DM_val, M_m_df1=NULL, M_m_df2=NULL,
##                          M_m_df3=NULL) {
##   # Define simulation parameters.
##   p = c(d_m=d_m)            # Dose rate (mg/kg/d) to human mother.
##   p = init_human_p(p)

##   # Run simulation.
##   out <- human_sim(p=p, M_m_df1=M_m_df1, M_m_df2=M_m_df2, M_m_df3=M_m_df3)

##   # Compute human dose metric.
##   DM_human = comp_dose_metric(out[ , "time"], out[ , "C_i"], out[ , "AUC_i"],
##                               p["t_gest"], p["t_lact"], DM_code)

##   # Cost is square of the difference between the human dose metric and the
##   # nominal dose metric value.
##   cost = (DM_human - DM_val) ** 2
##   return(cost)
## }

## # Legacy function
## cost_fun <- function(d_m, DM_code, DM_val, M_m_df1=NULL, M_m_df2=NULL,
##                      M_m_df3=NULL, M_i_df1=NULL, M_i_df2=NULL,
##                      preg=TRUE, human_age=NULL) {
##   # Define simulation parameters.
##   p = c(d_m=d_m)            # Dose rate (mg/kg/d) to human mother.
##   p = init_human_p(p)

##   # Define final age if preg is FALSE and human_age is non-Null.
##   if (preg) {
##     if (DM_code < 5) {
##       times=NULL
##     } else {
##       if (is.null(human_age)) {
##         t_end = 78 * 365          # Lifetime of filial human (d).
##       }
##       else {
##         t_end = human_age         # Lifetime of filial human (d).
##       }
##       t_start = -24.25 * 365        # Birth of mother (d).
##       times = seq(from=t_start, to=t_end, by=0.25)
##     }
##   } else {
##     if (is.null(human_age)) {
##       times = NULL
##     } else {
##       times = seq(from=0, to=human_age, by=1)
##     }
##   }
##   # Run simulation.
##   out <- human_sim(times=times, p=p, M_m_df1=M_m_df1, M_m_df2=M_m_df2,
##                    M_m_df3=M_m_df3, M_i_df1=M_i_df1, M_i_df2=M_i_df2,
##                    preg=preg)

##   # Compute human dose metric.
##   DM_human = comp_dose_metric(out[ , "time"], out[ , "C_i"], out[ , "AUC_i"],
##                               p[["t_gest"]], p[["t_lact"]], DM_code)

##   # Cost is square of the difference between the human dose metric and the
##   # nominal dose metric value.
##   cost = (DM_human - DM_val) ** 2
##   return(cost)
## }

## # Legacy function
## comp_hed <- function(DM_code, DM_val, M_m_df1=NULL, M_m_df2=NULL,
##                      M_m_df3=NULL, M_i_df1=NULL, M_i_df2=NULL,
##                      preg=TRUE, human_age=NULL) {
##   opt_result = optimize(cost_fun, c(0, 1e6), DM_code, DM_val,
##                         M_m_df1=M_m_df1, M_m_df2=M_m_df2, M_m_df3=M_m_df3,
##                         M_i_df1=M_i_df1, M_i_df2=M_i_df2, preg=preg,
##                         human_age=human_age)
##   return(opt_result$minimum)
## }


## # hed_demo <- function() {
## #   # Load model.
## #   load_model("lipophilic_tk")
## #
## #   # Define parameters for rat experiment.
## #   p_rat = c(d_m = 6,          # Dose rate (mg/kg/d) to mother.
## #             t_m_start = 6,    # Time (d) since conception maternal dose starts.
## #             t_m_end = 28)     # Time (d) since conception maternal dose ends.
## #   p_rat = init_rat_p(p_rat)
## #
## #   # Run simulation of rat experiment.
## #   out_rat <- rat_sim(p=p_rat)
## #
## #   # Compute dose metric #3: average concentration in infant during gestation
## #   # period.
## #   DM_code = 3
## #   DM_rat = comp_dose_metric(out_rat[ , "time"], out_rat[ , "C_i"],
## #                             out_rat[ , "AUC_i"], p_rat["t_gest"],
## #                             p_rat["t_lact"], DM_code)
## #
## #   # Compute human equivalent dose.
## #   M_m_df1 = read.csv("human_female_mass_0_to_20.csv", header=FALSE,
## #                      col.names=c("Age (y)", "Mass (kg)"), check.names=FALSE)
## #   M_m_df2 = read.csv("human_female_mass_20_to_25_pregnancy.csv",
## #                      header=FALSE, col.names=c("Age (y)", "Mass (kg)"),
## #                      check.names=FALSE)
## #   M_m_df3 = read.csv("human_female_mass_25_to_26_postnatal.csv",
## #                      header=FALSE, col.names=c("Age (y)", "Mass (kg)"),
## #                      check.names=FALSE)
## #   HED = comp_hed(DM_code, DM_rat, M_m_df1, M_m_df2, M_m_df3)
## #
## #   # Run simulation of human exposure at HED.
## #   p_human = c(d_m=HED)
## #   p_human = init_human_p(p_human)
## #   out_human <- human_sim(p=p_human, M_m_df1=M_m_df1, M_m_df2=M_m_df2,
## #                          M_m_df3=M_m_df3)
## #
## #   # Compute dose metric for human.
## #   DM_human = comp_dose_metric(out_human[ , "time"], out_human[ , "C_i"],
## #                               out_human[ , "AUC_i"], p_human["t_gest"],
## #                               p_human["t_lact"], DM_code)
## #
## #   cat("Rat dose", p_rat[["d_m"]], "mg/kg/d produces average concentration",
## #       DM_rat, "mg/kg in infant \nduring gestation.\n")
## #   cat("Human dose", p_human[["d_m"]], "mg/kg/d produces average concentration",
## #       DM_human, "mg/kg in infant \nduring gestation.\n")
## # }

## # Legacy function
## hed_demo1 <- function() {
##   # Load model.
##   load_model("lipophilic_tk")

##   # Get human maternal mass dataframes.
##   M_m_df1 = read.csv("human_female_mass_0_to_20.csv", header=FALSE,
##                      col.names=c("Age (y)", "Mass (kg)"), check.names=FALSE)
##   M_m_df2 = read.csv("human_female_mass_20_to_25_pregnancy.csv",
##                      header=FALSE, col.names=c("Age (y)", "Mass (kg)"),
##                      check.names=FALSE)
##   M_m_df3 = read.csv("human_female_mass_25_to_26_postnatal.csv",
##                      header=FALSE, col.names=c("Age (y)", "Mass (kg)"),
##                      check.names=FALSE)

##   # Get human infant mass dataframes.
##   M_i_df1 = read.csv("human_female_mass_0_to_20.csv", header=FALSE,
##                      col.names=c("Age (y)", "Mass (kg)"), check.names=FALSE)
##   M_i_df2 = read.csv("human_female_mass_20_to_75.csv", header=FALSE,
##                      col.names=c("Age (y)", "Mass (kg)"), check.names=FALSE)

##   # Define parameters for rat experiment.
##   p_rat = c(d_m = 6,          # Dose rate (mg/kg/d) to mother.
##             t_m_start = 6,    # Time (d) since conception maternal dose starts.
##             t_m_end = 28,     # Time (d) since conception maternal dose ends.
##             t_i_end = 50)     # Time (d) since conception infant dose ends.
##   p_rat = init_rat_p(p_rat)

##   # Run simulation of rat experiment.
##   out_rat <- rat_sim(p=p_rat)

##   # # Plot all output variables vs. time.
##   # out = out_rat
##   # for (idx in 2:ncol(out)) {
##   #     plot(out[ , 1], out[ , idx], type='l', xlab="Time (d)",
##   #          ylab=colnames(out)[idx])
##   # }

##   # Dose metric codes and corresponding descriptions
##   DM_code_vec = c(1, 2, 3, 4, 5)
##   desc_vec1 = c("peak concentration",
##                 "average concentration",
##                 "average concentration",
##                 "average concentration",
##                 "average concentration")
##   desc_vec2 = c("fetus/infant",
##                 "fetus/infant",
##                 "fetus",
##                 "infant",
##                 "fetus/infant")
##   desc_vec3 = c("gestation and nursing",
##                 "gestation and nursing",
##                 "gestation",
##                 "nursing",
##                 "the entire simulation")
##   for (idx in seq(length(DM_code_vec))) {
##     DM_code = DM_code_vec[idx]
##     DM_rat = comp_dose_metric(out_rat[ , "time"], out_rat[ , "C_i"],
##                               out_rat[ , "AUC_i"], p_rat["t_gest"],
##                               p_rat["t_lact"], DM_code)

##     # Compute human equivalent dose.
##     HED = comp_hed(DM_code, DM_rat, M_m_df1=M_m_df1, M_m_df2=M_m_df2,
##                    M_m_df3=M_m_df3, M_i_df1=M_i_df1, M_i_df2=M_i_df2)

##     # Run simulation of human exposure at HED.
##     p_human = c(d_m=HED)
##     p_human = init_human_p(p_human)
##     out_human <- human_sim(p=p_human, M_m_df1=M_m_df1, M_m_df2=M_m_df2,
##                            M_m_df3=M_m_df3, M_i_df1=M_i_df1, M_i_df2=M_i_df2)

##     # Compute dose metric for human.
##     DM_human = comp_dose_metric(out_human[ , "time"], out_human[ , "C_i"],
##                                 out_human[ , "AUC_i"], p_human["t_gest"],
##                                 p_human["t_lact"], DM_code)

##     cat("Dose Metric #", DM_code, ":\n", sep="")
##     cat("   -Rat dose", p_rat[["d_m"]], "mg/kg/d produces", desc_vec1[idx],
##         DM_rat, "mg/kg in", desc_vec2[idx], "\n    during", desc_vec3[idx],
##         "\n")
##     cat("   -Human dose", p_human[["d_m"]], "mg/kg/d produces", desc_vec1[idx],
##         DM_human, "mg/kg in", desc_vec2[idx], "\n    during", desc_vec3[idx],
##         "\n")
##   }
## }
