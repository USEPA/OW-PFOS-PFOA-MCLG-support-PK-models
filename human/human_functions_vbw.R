human_child_v <- function(child_age=3,sex=NULL,chem=NULL,dose=0.001){
    if(is.null(sex)){
        print("Sex of child not specified.")
        return(NULL)
    }
    if (chem=="PFOA"){
        p_PFOA=c("CL"=0.085/1000,"Vd"=0.170, "CLr"=0.085/1000, "P_milk"=0.0577, "r_f_m"=0.783, "pow_milk"=1) # PFOA parameters
        p=c(p_PFOA,"d_m"=dose)# 1 ug/kg/d
        r_m_mf="r_m_mf_pfoa_v.csv"
    } else if (chem=="PFOS"){
        p_PFOS=c("CL"=0.081/1000,"Vd"=0.230, "CLr"=0.081/1000, "P_milk"=0.0138, "r_f_m"=0.454, "pow_milk"=1) # PFOS parameters
        p=c(p_PFOS,"d_m"=dose)# 1 ug/kg/d
        r_m_mf="r_m_mf_pfos_v.csv"
    } else {
        print("Specify PFOA or PFOS.")
        return(NULL)
    }
    out <- human_mother_child_v(age_end=child_age+25,p=p,sex=sex,mName="pfas_tk",
                                milkin="r_milk_verner1.csv",r_m_mf=r_m_mf)
    out <- cbind(out,(out[,"time"]/365)-25)
    colnames(out)[ncol(out)] <- "years"
    return(out)
}

human_mother_child_v <- function(age_start=0, a_start=0, age_end=30, p=NULL, sex="Female",
                               mName="pfas_tk", milkin="r_milk.csv", r_m_mf=NULL) {
  t_start = age_start*365
  t_end = 25*365 # 25 years of age, end of pregnancy, converted to days
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


  M_m_df1 = read.csv("human_female_mass_0_to_25_v.csv", header=FALSE,
                    col.names=c("Age", "Mass (kg)"),
                    check.names=FALSE)

  M_m_df1[,1] = M_m_df1[,1]*365 # Convert age from years to days

  M_m_df2 =read.csv("human_female_mass_25_to_35_v.csv",
                  header=FALSE, col.names=c("Age", "Mass (kg)"),
                  check.names=FALSE)
  M_m_df2[,1] = M_m_df2[,1]*365 # Convert age from years to days

  if (sex=="Female"){
      M_i_df1 = read.csv("human_female_mass_0_to_20_v.csv", header=FALSE,
                         col.names=c("Age (y)", "Mass (kg)"),
                         check.names=FALSE)

      M_i_df2 = read.csv("human_female_mass_20_to_75_v.csv", header=FALSE,
                         col.names=c("Age (y)", "Mass (kg)"),
                         check.names=FALSE)
  } else if (sex=="Male"){
      M_i_df1 = read.csv("human_male_mass_0_to_20_v.csv", header=FALSE,
                         col.names=c("Age (y)", "Mass (kg)"),
                         check.names=FALSE)
      M_i_df2 = read.csv("human_male_mass_20_to_75_v.csv", header=FALSE,
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
  if(is.null(r_m_mf)){
      r_m_mf_in = cbind(times=c(0,24.25*365,t_end), ratio=c(1,1,1))
  } else {
      r_m_mf_data = read.csv(r_m_mf, header=FALSE,
                         col.names=c("Age (y)", "Ratio"),
                         check.names=FALSE)
      r_m_mf_in = cbind(times=((r_m_mf_data[["Age (y)"]]*365)+t_start),
                        rate=r_m_mf_data[["Ratio"]])
  }


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

  waterin_m="water_95_lact.csv"
  waterin_i="water_95_bf.csv"

  DW_input = read.csv(waterin_m, header=FALSE,
                      col.names=c("Age (y)", "DW"),
                      check.names=FALSE)

  DW_m_in = cbind(times = DW_input[,1]*365,DW_m_in = DW_input[,2])

  DW_input = read.csv(waterin_i, header=FALSE,
                      col.names=c("Age (y)", "DW"),
                      check.names=FALSE)

  DW_i_in = cbind(times = (DW_input[,1]+25)*365,DW_i_in = DW_input[,2])

  df_events = rbind(df_CLr)

  ## Run model simulation to birth
  out1 = run_model(mName, times, Y0, parms, rtol=1e-6, atol=1e-6,
                   list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in, d_i_in,
                        Vdaf_m_in, Vdaf_i_in, DW_m_in, DW_i_in),
                   list(method="linear", rule=2, ties="ordered"),
                   list(data=df_CLr))

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
  Y0["A_i"] =  out1[[nrow(out1), "C_i"]] * parms[["Vd"]] * M_i_in[1,"mass"]
  ## Y0["A_mf"] = out1[[nrow(out1),"A_mf"]]*M_mf_in[1,"mass"]/(M_mf_in[1,"mass"]+M_i_in[1,"mass"])
  ## Y0["A_i"] = out1[[nrow(out1),"A_mf"]]*M_i_in[1,"mass"]/(M_mf_in[1,"mass"]+M_i_in[1,"mass"])
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
  out2 = run_model(mName, times2, Y0, parms, rtol=1e-6, atol=1e-6,
                   list(M_mf_in, M_i_in, R_milk_in, r_m_mf_in, d_i_in,
                        Vdaf_m_in, Vdaf_i_in, DW_m_in, DW_i_in),
                   list(method="linear", rule=2, ties="ordered"),
                   list(data=df_CLr))
  out = rbind(out1[1:nrow(out1)-1, ], out2)

  return(out)   # Return output data frame.
}
