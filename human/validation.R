##------------------------------------------------------------------------------
## validation.R
## File to compare output of human developmental model to selected datasets.
## For validation, PFOA/S exposure levels were set to match maternal serum levels at delivery.
## Maternal levels in the Mogensen 2015 study were back-calculated from the levels they reported
## in children at delivery, which they estimated from maternal levels. Model predictions are
## compared to measured values at 6 and 19 months (Fromme 2010), 36 mo. (Granum 2013),
## and at 11, 19, and 60 mo. (Mogensen 2015).
##
##------------------------------------------------------------------------------

## Load required functions and libraries.
source("compile_pfas_tk.R")
source("sim_functions.R")
source("human_functions.R")

## Define external data
## PFOA

## Fromme 2010 Study (Table 1 of Verner 2015)
Cm_Fromme_PFOA <- 2.4
Ci_6 <- 8.7
Ci_6sd <- 5.6
Ci_19 <- 5.7
Ci_19sd <- 2.5
fromme_data <- data.frame(age=c(6,19),"Mean PFOA"=c(Ci_6,Ci_19),"SD PFOA"=c(Ci_6sd,Ci_19sd))

## Nowegian Study (MOBA) (Table 1 of Verner 2015) (Granum 2013 and Magnus 2006)
Cm_MOBA_PFOA <- 1.1
Ci_36 <- 2.8
Ci_36sd <- 1.0
moba_data <- data.frame(age=c(36),"Mean PFOA"=c(Ci_36),"SD PFOA"=c(Ci_36sd))

##PFOS data
## Fromme Study (Table 1 of Verner 2015)
Cm_Fromme_PFOS <- 3.5
Ci_6 <- 3.6
Ci_6sd <- 2.1
Ci_19 <- 2.4
Ci_19sd <- 1.5

fromme_data <- cbind(fromme_data,"Mean PFOS"=c(Ci_6,Ci_19),"SD PFOS"=c(Ci_6sd,Ci_19sd))

## Nowegian Study (MOBA) (Table 1 of Verner 2015)
Cm_MOBA_PFOS <- 5.4
Ci_36 <- 5.0
Ci_36sd <- 2.0

moba_data <- cbind(moba_data,"Mean PFOS"=c(Ci_36),"SD PFOS"=c(Ci_36sd))

## Mogensen 2015 (means/sd extracted from Figure 1)
mogensen_data <- read.csv("Mogensen Fig1.csv")[1:4,6:10]
## For mogensen, child levels at age 0 are estimate by multiplying maternal levels by 0.34 (PFOA) or 0.72 (PFOS).# These are backwards in Mogensen!!!!! See Needham, 2011 (original source)
Ci0_mogensen_PFOA <- mogensen_data[mogensen_data[,"Age"]==0,3]
Ci0_mogensen_PFOS <- mogensen_data[mogensen_data[,"Age"]==0,2]
Cm_mogensen_PFOA <- Ci0_mogensen_PFOA/0.72
Cm_mogensen_PFOS <- Ci0_mogensen_PFOS/0.34

## Trims out 1st point based on maternal serum levels
mogensen_data <- mogensen_data[2:4,]

DW=FALSE
if(DW){
    dose_scale <- 0.2E5# For DW exposure
} else {
    dose_scale <- 1E3## For /kg exposure
}

## Use reverse dosimetry on maternal levels at delivery to set Fromme exposure
Fromme_PFOA_dose <- child_rev(Cm_Fromme_PFOA,age=0,sex="Male",chem="PFOA",maternal=T,DW=DW,
                              waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
## Run model with Fromme_PFOA_dose.
Fromme_PFOA_data <- human_child(child_age=19/12,sex="Male",chem="PFOA",dose=Fromme_PFOA_dose,DW=DW,
                                waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
## Add a column with C_i normalized by dose (and some scalar)
Fromme_PFOA_data <- data.frame(Fromme_PFOA_data,C_i_n=Fromme_PFOA_data[,"C_i"]/(Fromme_PFOA_dose*dose_scale))

## Optional plots of Fromme data
## par(mfrow=c(1,1),mar=c(4.1,4.1,2.1,2.1),lwd=2)
## plot(Fromme_PFOA_data[,"years"],Fromme_PFOA_data[,"C_i"],type='l',xlim=c(0,19/12),ylim=c(0,max(Fromme_PFOA_data[,"C_i"],fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"])))
## lines(Fromme_PFOA_data2[,"years"],Fromme_PFOA_data2[,"C_i"],col=2)
## points(fromme_data[,"age"]/12,fromme_data[,"Mean.PFOA"])
## arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,y0=fromme_data[,"Mean.PFOA"]-fromme_data[,"SD.PFOA"],y1=fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"],code=3,angle=90,length=0.2)

## plot(Fromme_PFOA_data2[,"years"],Fromme_PFOA_data2[,"C_i"],type='l',xlim=c(0,19/12),ylim=c(0,max(Fromme_PFOA_data2[,"C_i_n"],fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"])))
## lines(Fromme_PFOA_data2[,"years"],Fromme_PFOA_data2[,"C_i"],type='l',xlim=c(0,19/12),ylim=c(0,max(Fromme_PFOA_data2[,"C_i"],fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"])),col=2)
## points(fromme_data[,"age"]/12,fromme_data[,"Mean.PFOA"])
## arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,y0=fromme_data[,"Mean.PFOA"]-fromme_data[,"SD.PFOA"],y1=fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"],code=3,angle=90,length=0.2)

## plot(Fromme_PFOA_data[,"years"],Fromme_PFOA_data[,"C_m"],type='l',xlim=c(-1,1))
## points(0,Cm_Fromme_PFOA)

## Calculation for MOBA_PFOA as was performed previously for Fromme_PFOA
MOBA_PFOA_dose <- child_rev(Cm_MOBA_PFOA,age=0,sex="Male",chem="PFOA",maternal=T,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOA_data <- human_child(child_age=36/12,sex="Male",chem="PFOA",dose=MOBA_PFOA_dose,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOA_data <- data.frame(MOBA_PFOA_data,C_i_n=MOBA_PFOA_data[,"C_i"]/(MOBA_PFOA_dose*dose_scale))

## Optional plots
## plot(MOBA_PFOA_data[,"years"],MOBA_PFOA_data[,"C_i"],type='l',xlim=c(0,36/12),ylim=c(0,max(MOBA_PFOA_data[,"C_i"],fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"])))
## points(moba_data[,"age"]/12,moba_data[,"Mean.PFOA"])
## arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,y0=moba_data[,"Mean.PFOA"]-moba_data[,"SD.PFOA"],y1=moba_data[,"Mean.PFOA"]+moba_data[,"SD.PFOA"],code=3,angle=90,length=0.2)

## Calculation for Mogensen_PFOA as was performed previously for Fromme_PFOA. Reverse dosimitry was based on measured maternal data. Their estimate of child level (which they present in the paper) was based on multiplying the measured maternal level by a ratio of cord to maternal serum.
Mogensen_PFOA_dose <- child_rev(Cm_mogensen_PFOA,age=0,sex="Male",chem="PFOA",maternal=T,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOA_data <- human_child(child_age=60/12,sex="Male",chem="PFOA",dose=Mogensen_PFOA_dose,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOA_data <- data.frame(Mogensen_PFOA_data,C_i_n=Mogensen_PFOA_data[,"C_i"]/(Mogensen_PFOA_dose*dose_scale))

## Optional plot
## plot(Mogensen_PFOA_data[,"years"],Mogensen_PFOA_data[,"C_i"],type='l',xlim=c(0,60/12),ylim=c(0,max(Mogensen_PFOA_data[,"C_i"],fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"])))
## points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOA"])
## arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,y0=mogensen_data[,"Mean.PFOA"]-mogensen_data[,"SD.PFOA"],y1=mogensen_data[,"Mean.PFOA"]+mogensen_data[,"SD.PFOA"],code=3,angle=90,length=0.2)

## Analogous PFOS calculations:
Fromme_PFOS_dose <- child_rev(Cm_Fromme_PFOS,age=0,sex="Male",chem="PFOS",maternal=T,DW=DW,waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
Fromme_PFOS_data <- human_child(child_age=19/12,sex="Male",chem="PFOS",dose=Fromme_PFOS_dose,DW=DW,waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
Fromme_PFOS_data <- data.frame(Fromme_PFOS_data,C_i_n=Fromme_PFOS_data[,"C_i"]/(Fromme_PFOS_dose*dose_scale))

## plot(Fromme_PFOS_data[,"years"],Fromme_PFOS_data[,"C_i"],type='l',xlim=c(0,19/12),ylim=c(0,max(Fromme_PFOS_data[,"C_i"],fromme_data[,"Mean PFOS"]+fromme_data[,"SD PFOS"])))
## points(fromme_data[,"age"]/12,fromme_data[,"Mean PFOS"])
## arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,y0=fromme_data[,"Mean PFOS"]-fromme_data[,"SD PFOS"],y1=fromme_data[,"Mean PFOS"]+fromme_data[,"SD PFOS"],code=3,angle=90,length=0.2)

## plot(Fromme_PFOS_data[,"years"],Fromme_PFOS_data[,"C_m"],type='l',xlim=c(-1,1))
## points(0,Cm_Fromme_PFOS)

MOBA_PFOS_dose <- child_rev(Cm_MOBA_PFOS,age=0,sex="Male",chem="PFOS",maternal=T,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOS_data <- human_child(child_age=36/12,sex="Male",chem="PFOS",dose=MOBA_PFOS_dose,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOS_data <- data.frame(MOBA_PFOS_data,C_i_n=MOBA_PFOS_data[,"C_i"]/(MOBA_PFOS_dose*dose_scale))

## plot(MOBA_PFOS_data[,"years"],MOBA_PFOS_data[,"C_i"],type='l',xlim=c(0,36/12),ylim=c(0,max(MOBA_PFOS_data[,"C_i"],moba_data[,"Mean PFOS"]+moba_data[,"SD PFOS"])))
## points(moba_data[,"age"]/12,moba_data[,"Mean PFOS"])
## arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,y0=moba_data[,"Mean PFOS"]-moba_data[,"SD PFOS"],y1=moba_data[,"Mean PFOS"]+moba_data[,"SD PFOS"],code=3,angle=90,length=0.2)

Mogensen_PFOS_dose <- child_rev(Cm_mogensen_PFOS,age=0,sex="Male",chem="PFOS",maternal=T,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOS_data <- human_child(child_age=60/12,sex="Male",chem="PFOS",dose=Mogensen_PFOS_dose,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOS_data <- data.frame(Mogensen_PFOS_data,C_i_n=Mogensen_PFOS_data[,"C_i"]/(Mogensen_PFOS_dose*dose_scale))

## plot(Mogensen_PFOS_data[,"years"],Mogensen_PFOS_data[,"C_i"],type='l',xlim=c(0,60/12),ylim=c(0,max(Mogensen_PFOS_data[,"C_i"],mogensen_data[,"Mean.PFOS"]+mogensen_data[,"SD.PFOS"])))
## points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOS"])
## arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,y0=mogensen_data[,"Mean.PFOS"]-mogensen_data[,"SD.PFOS"],y1=mogensen_data[,"Mean.PFOS"]+mogensen_data[,"SD.PFOS"],code=3,angle=90,length=0.2)


## Plot for PFOA and PFOS that is normalized by dose
## ## svg(filename="validation.svg",height=7,width=14)
## par(mfrow=c(1,2),mar=c(4.1,4.1,2.1,2.1),lwd=2)
## ## Dose normalized plot - PFOA
## plot(Fromme_PFOA_data[,"years"],Fromme_PFOA_data[,"C_i_n"],type='l',xlim=c(0,60/12),ylim=c(0,40),xlab="Age (yr)",ylab="PFOA (ng/ml)",col=2)
## lines(MOBA_PFOA_data[,"years"],MOBA_PFOA_data[,"C_i_n"],col=4)
## lines(Mogensen_PFOA_data[,"years"],Mogensen_PFOA_data[,"C_i_n"],col=3,lty=3)
## points(fromme_data[,"age"]/12,fromme_data[,"Mean.PFOA"]/(Fromme_PFOA_dose*dose_scale),col=2,pch=16)
## arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,
##        y0=(fromme_data[,"Mean.PFOA"]-fromme_data[,"SD.PFOA"])/(Fromme_PFOA_dose*dose_scale),
##        y1=(fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"])/(Fromme_PFOA_dose*dose_scale),code=3,angle=90,length=0.2,col=2)
## points(moba_data[,"age"]/12,moba_data[,"Mean.PFOA"]/(MOBA_PFOA_dose*dose_scale),col=4,pch=17)
## arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,
##        y0=(moba_data[,"Mean.PFOA"]-moba_data[,"SD.PFOA"])/(MOBA_PFOA_dose*dose_scale),
##        y1=(moba_data[,"Mean.PFOA"]+moba_data[,"SD.PFOA"])/(MOBA_PFOA_dose*dose_scale),code=3,angle=90,length=0.2,col=4)
## points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOA"]/(Mogensen_PFOA_dose*dose_scale),col=3,pch=18)
## arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,
##        y0=(mogensen_data[,"Mean.PFOA"]-mogensen_data[,"SD.PFOA"])/(Mogensen_PFOA_dose*dose_scale),
##        y1=(mogensen_data[,"Mean.PFOA"]+mogensen_data[,"SD.PFOA"])/(Mogensen_PFOA_dose*dose_scale),code=3,angle=90,length=0.2,col=3)
## legend('topright',legend=c("Fromme, 2010","MOBA Cohort","Mogensen 2015"),pch=c(16,17,18),col=c(2,4,3))

## ## Dose normalized plot - PFOS
## plot(Fromme_PFOS_data[,"years"],Fromme_PFOS_data[,"C_i_n"],type='l',xlim=c(0,60/12),ylim=c(0,30),xlab="Age (yr)",ylab="PFOS (ng/ml)",col=2)
## lines(MOBA_PFOS_data[,"years"],MOBA_PFOS_data[,"C_i_n"],col=4)
## lines(Mogensen_PFOS_data[,"years"],Mogensen_PFOS_data[,"C_i_n"],col=3,lty=3)
## points(fromme_data[,"age"]/12,fromme_data[,"Mean PFOS"]/(Fromme_PFOS_dose*dose_scale),col=2,pch=16)
## arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,
##        y0=(fromme_data[,"Mean PFOS"]-fromme_data[,"SD PFOS"])/(Fromme_PFOS_dose*dose_scale),
##        y1=(fromme_data[,"Mean PFOS"]+fromme_data[,"SD PFOS"])/(Fromme_PFOS_dose*dose_scale),code=3,angle=90,length=0.2,col=2)
## points(moba_data[,"age"]/12,moba_data[,"Mean PFOS"]/(MOBA_PFOS_dose*dose_scale),col=4,pch=17)
## arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,
##        y0=(moba_data[,"Mean PFOS"]-moba_data[,"SD PFOS"])/(MOBA_PFOS_dose*dose_scale),
##        y1=(moba_data[,"Mean PFOS"]+moba_data[,"SD PFOS"])/(MOBA_PFOS_dose*dose_scale),code=3,angle=90,length=0.2,col=4)
## points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOS"]/(Mogensen_PFOS_dose*dose_scale),col=3,pch=18)
## arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,
##        y0=(mogensen_data[,"Mean.PFOS"]-mogensen_data[,"SD.PFOS"])/(Mogensen_PFOS_dose*dose_scale),
##        y1=(mogensen_data[,"Mean.PFOS"]+mogensen_data[,"SD.PFOS"])/(Mogensen_PFOS_dose*dose_scale),code=3,angle=90,length=0.2,col=3)
## legend('topright',legend=c("Fromme, 2010","MOBA Cohort","Mogensen 2015"),pch=c(16,17,18),col=c(2,4,3))
## ## dev.off()

## Plot for PFOA and PFOS that is not normalized by dose. This is used in the Appendix.
## svg(filename="validation.svg",height=7,width=14)
png(filename="validation.png",height=7,width=14,units="in",res=600)
par(mfrow=c(1,2),mar=c(4.1,4.1,2.1,2.1),lwd=2)
## Dose unnormalized plot - PFOA
plot(Fromme_PFOA_data[,"years"],Fromme_PFOA_data[,"C_i"],type='l',xlim=c(0,60/12),ylim=c(0,15),xlab="Age (yr)",ylab="PFOA (ng/ml)",col=2)
lines(MOBA_PFOA_data[,"years"],MOBA_PFOA_data[,"C_i"],col=4)
lines(Mogensen_PFOA_data[,"years"],Mogensen_PFOA_data[,"C_i"],col=3,lty=3)
points(fromme_data[,"age"]/12,fromme_data[,"Mean.PFOA"],col=2,pch=16)
arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,
       y0=(fromme_data[,"Mean.PFOA"]-fromme_data[,"SD.PFOA"]),
       y1=(fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"]),code=3,angle=90,length=0.2,col=2)
points(moba_data[,"age"]/12,moba_data[,"Mean.PFOA"],col=4,pch=17)
arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,
       y0=(moba_data[,"Mean.PFOA"]-moba_data[,"SD.PFOA"]),
       y1=(moba_data[,"Mean.PFOA"]+moba_data[,"SD.PFOA"]),code=3,angle=90,length=0.2,col=4)
points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOA"],col=3,pch=18)
arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,
       y0=(mogensen_data[,"Mean.PFOA"]-mogensen_data[,"SD.PFOA"]),
       y1=(mogensen_data[,"Mean.PFOA"]+mogensen_data[,"SD.PFOA"]),code=3,angle=90,length=0.2,col=3)
legend('topright',legend=c("Fromme, 2010","MOBA Cohort","Mogensen, 2015"),pch=c(16,17,18),col=c(2,4,3))

## Dose unnormalized plot - PFOS
plot(Fromme_PFOS_data[,"years"],Fromme_PFOS_data[,"C_i"],type='l',xlim=c(0,60/12),ylim=c(0,40),xlab="Age (yr)",ylab="PFOS (ng/ml)",col=2)
lines(MOBA_PFOS_data[,"years"],MOBA_PFOS_data[,"C_i"],col=4)
lines(Mogensen_PFOS_data[,"years"],Mogensen_PFOS_data[,"C_i"],col=3,lty=3)
points(fromme_data[,"age"]/12,fromme_data[,"Mean PFOS"],col=2,pch=16)
arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,
       y0=(fromme_data[,"Mean PFOS"]-fromme_data[,"SD PFOS"]),
       y1=(fromme_data[,"Mean PFOS"]+fromme_data[,"SD PFOS"]),code=3,angle=90,length=0.2,col=2)
points(moba_data[,"age"]/12,moba_data[,"Mean PFOS"],col=4,pch=17)
arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,
       y0=(moba_data[,"Mean PFOS"]-moba_data[,"SD PFOS"]),
       y1=(moba_data[,"Mean PFOS"]+moba_data[,"SD PFOS"]),code=3,angle=90,length=0.2,col=4)
points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOS"],col=3,pch=18)
arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,
       y0=(mogensen_data[,"Mean.PFOS"]-mogensen_data[,"SD.PFOS"]),
       y1=(mogensen_data[,"Mean.PFOS"]+mogensen_data[,"SD.PFOS"]),code=3,angle=90,length=0.2,col=3)
legend('topright',legend=c("Fromme, 2010","MOBA Cohort","Mogensen, 2015"),pch=c(16,17,18),col=c(2,4,3))
dev.off()

## Wrangling of data to plot all 3 studies together:
fromme_data$pre.PFOA <- c(human_child_conc(child_age=fromme_data[1,1]/12,sex="Male",chem="PFOA",dose=Fromme_PFOA_dose,DW=DW,
                                           waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv"),
                          human_child_conc(child_age=fromme_data[2,1]/12,sex="Male",chem="PFOA",dose=Fromme_PFOA_dose,DW=DW,
                                           waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv"))
moba_data$pre.PFOA <- c(human_child_conc(child_age=moba_data[1,1]/12,sex="Male",chem="PFOA",dose=MOBA_PFOA_dose,DW=DW,
                                           waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv"))
mogensen_data$pre.PFOA <- c(human_child_conc(child_age=mogensen_data[1,1]/12,sex="Male",chem="PFOA",dose=Mogensen_PFOA_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"),
                            human_child_conc(child_age=mogensen_data[2,1]/12,sex="Male",chem="PFOA",dose=Mogensen_PFOA_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"),
                            human_child_conc(child_age=mogensen_data[3,1]/12,sex="Male",chem="PFOA",dose=Mogensen_PFOA_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"))

fromme_data$pre.PFOS <- c(human_child_conc(child_age=fromme_data[1,1]/12,sex="Male",chem="PFOS",dose=Fromme_PFOS_dose,DW=DW,
                                           waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv"),
                          human_child_conc(child_age=fromme_data[2,1]/12,sex="Male",chem="PFOS",dose=Fromme_PFOS_dose,DW=DW,
                                           waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv"))
moba_data$pre.PFOS <- c(human_child_conc(child_age=moba_data[1,1]/12,sex="Male",chem="PFOS",dose=MOBA_PFOS_dose,DW=DW,
                                           waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv"))
mogensen_data$pre.PFOS <- c(human_child_conc(child_age=mogensen_data[1,1]/12,sex="Male",chem="PFOS",dose=Mogensen_PFOS_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"),
                            human_child_conc(child_age=mogensen_data[2,1]/12,sex="Male",chem="PFOS",dose=Mogensen_PFOS_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"),
                            human_child_conc(child_age=mogensen_data[3,1]/12,sex="Male",chem="PFOS",dose=Mogensen_PFOS_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"))

pl_data_pfoa <- data.frame(observed=c(fromme_data[,2],moba_data[,2],mogensen_data[,3]),
                           observed_uci=c(fromme_data[,2]+fromme_data[,3],moba_data[,2]+moba_data[,3],mogensen_data[,3]+mogensen_data[,5]),
                           observed_lci=c(fromme_data[,2]-fromme_data[,3],moba_data[,2]-moba_data[,3],mogensen_data[,3]-mogensen_data[,5]),
                           predicted=c(fromme_data[,6],moba_data[,6],mogensen_data[,6]))
pl_data_pfos <- data.frame(observed=c(fromme_data[,4],moba_data[,4],mogensen_data[,2]),
                           observed_uci=c(fromme_data[,4]+fromme_data[,5],moba_data[,4]+moba_data[,5],mogensen_data[,2]+mogensen_data[,4]),
                           observed_lci=c(fromme_data[,4]-fromme_data[,5],moba_data[,4]-moba_data[,5],mogensen_data[,2]-mogensen_data[,4]),
                           predicted=c(fromme_data[,7],moba_data[,7],mogensen_data[,7]))

## Plot observed vs predicted values # In Appendix
## svg(filename="validation_CvC.svg",height=7,width=14)
png(filename="validation_CvC.png",height=7,width=14,units="in",res=600)
par(mfrow=c(1,2),mar=c(4.1,4.1,2.1,2.1),lwd=2)
plot(pl_data_pfoa$observed,pl_data_pfoa$predicted,xlim=c(0,12),ylim=c(0,12),xlab="Observed PFOA (ng/ml)",ylab="Predicted PFOA (ng/ml)",xaxs='i',yaxs='i')
abline(a=0,b=1)
plot(pl_data_pfos$observed,pl_data_pfos$predicted,xlim=c(0,30),ylim=c(0,30),xlab="Observed PFOS (ng/ml)",ylab="Predicted PFOS (ng/ml)",xaxs='i',yaxs='i')
abline(a=0,b=1)
dev.off()

## Save relative changes for pfoa/pfos with constant Vd
error_pfoa_vdcon <- (pl_data_pfoa$predicted-pl_data_pfoa$observed)/pl_data_pfoa$observed
error_pfos_vdcon <- (pl_data_pfos$predicted-pl_data_pfos$observed)/pl_data_pfos$observed

##------------------------------------------------------------------------
DW=TRUE #With exposure in terms of a constant drinking water level.

## Use reverse dosimetry on maternal levels at delivery to set Fromme exposure
Fromme_PFOA_dose <- child_rev(Cm_Fromme_PFOA,age=0,sex="Male",chem="PFOA",maternal=T,DW=DW,
                              waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
## Run model with Fromme_PFOA_dose.
Fromme_PFOA_data <- human_child(child_age=19/12,sex="Male",chem="PFOA",dose=Fromme_PFOA_dose,DW=DW,
                                waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
## Add a column with C_i normalized by dose (and some scalar)
Fromme_PFOA_data <- data.frame(Fromme_PFOA_data,C_i_n=Fromme_PFOA_data[,"C_i"]/(Fromme_PFOA_dose*dose_scale))

## Calculation for MOBA_PFOA as was performed previously for Fromme_PFOA
MOBA_PFOA_dose <- child_rev(Cm_MOBA_PFOA,age=0,sex="Male",chem="PFOA",maternal=T,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOA_data <- human_child(child_age=36/12,sex="Male",chem="PFOA",dose=MOBA_PFOA_dose,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOA_data <- data.frame(MOBA_PFOA_data,C_i_n=MOBA_PFOA_data[,"C_i"]/(MOBA_PFOA_dose*dose_scale))

## Calculation for Mogensen_PFOA as was performed previously for Fromme_PFOA. Reverse dosimitry was based on measured maternal data. Their estimate of child level (which they present in the paper) was based on multiplying the measured maternal level by a ratio of cord to maternal serum.
Mogensen_PFOA_dose <- child_rev(Cm_mogensen_PFOA,age=0,sex="Male",chem="PFOA",maternal=T,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOA_data <- human_child(child_age=60/12,sex="Male",chem="PFOA",dose=Mogensen_PFOA_dose,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOA_data <- data.frame(Mogensen_PFOA_data,C_i_n=Mogensen_PFOA_data[,"C_i"]/(Mogensen_PFOA_dose*dose_scale))

## Analogous PFOS calculations (with DW=TRUE):
Fromme_PFOS_dose <- child_rev(Cm_Fromme_PFOS,age=0,sex="Male",chem="PFOS",maternal=T,DW=DW,waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
Fromme_PFOS_data <- human_child(child_age=19/12,sex="Male",chem="PFOS",dose=Fromme_PFOS_dose,DW=DW,waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
Fromme_PFOS_data <- data.frame(Fromme_PFOS_data,C_i_n=Fromme_PFOS_data[,"C_i"]/(Fromme_PFOS_dose*dose_scale))

MOBA_PFOS_dose <- child_rev(Cm_MOBA_PFOS,age=0,sex="Male",chem="PFOS",maternal=T,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOS_data <- human_child(child_age=36/12,sex="Male",chem="PFOS",dose=MOBA_PFOS_dose,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOS_data <- data.frame(MOBA_PFOS_data,C_i_n=MOBA_PFOS_data[,"C_i"]/(MOBA_PFOS_dose*dose_scale))

Mogensen_PFOS_dose <- child_rev(Cm_mogensen_PFOS,age=0,sex="Male",chem="PFOS",maternal=T,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOS_data <- human_child(child_age=60/12,sex="Male",chem="PFOS",dose=Mogensen_PFOS_dose,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOS_data <- data.frame(Mogensen_PFOS_data,C_i_n=Mogensen_PFOS_data[,"C_i"]/(Mogensen_PFOS_dose*dose_scale))

## Plot for PFOA and PFOS that is not normalized by dose. This is used in the Appendix for the DW exposure example.
## svg(filename="validation.svg",height=7,width=14)
png(filename="validationDW.png",height=7,width=14,units="in",res=600)
par(mfrow=c(1,2),mar=c(4.1,4.1,2.1,2.1),lwd=2)
## Dose unnormalized plot - PFOA
plot(Fromme_PFOA_data[,"years"],Fromme_PFOA_data[,"C_i"],type='l',xlim=c(0,60/12),ylim=c(0,15),xlab="Age (yr)",ylab="PFOA (ng/ml)",col=2)
lines(MOBA_PFOA_data[,"years"],MOBA_PFOA_data[,"C_i"],col=4)
lines(Mogensen_PFOA_data[,"years"],Mogensen_PFOA_data[,"C_i"],col=3,lty=3)
points(fromme_data[,"age"]/12,fromme_data[,"Mean.PFOA"],col=2,pch=16)
arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,
       y0=(fromme_data[,"Mean.PFOA"]-fromme_data[,"SD.PFOA"]),
       y1=(fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"]),code=3,angle=90,length=0.2,col=2)
points(moba_data[,"age"]/12,moba_data[,"Mean.PFOA"],col=4,pch=17)
arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,
       y0=(moba_data[,"Mean.PFOA"]-moba_data[,"SD.PFOA"]),
       y1=(moba_data[,"Mean.PFOA"]+moba_data[,"SD.PFOA"]),code=3,angle=90,length=0.2,col=4)
points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOA"],col=3,pch=18)
arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,
       y0=(mogensen_data[,"Mean.PFOA"]-mogensen_data[,"SD.PFOA"]),
       y1=(mogensen_data[,"Mean.PFOA"]+mogensen_data[,"SD.PFOA"]),code=3,angle=90,length=0.2,col=3)
legend('topright',legend=c("Fromme, 2010","MOBA Cohort","Mogensen, 2015"),pch=c(16,17,18),col=c(2,4,3))

## Dose unnormalized plot - PFOS
plot(Fromme_PFOS_data[,"years"],Fromme_PFOS_data[,"C_i"],type='l',xlim=c(0,60/12),ylim=c(0,40),xlab="Age (yr)",ylab="PFOS (ng/ml)",col=2)
lines(MOBA_PFOS_data[,"years"],MOBA_PFOS_data[,"C_i"],col=4)
lines(Mogensen_PFOS_data[,"years"],Mogensen_PFOS_data[,"C_i"],col=3,lty=3)
points(fromme_data[,"age"]/12,fromme_data[,"Mean PFOS"],col=2,pch=16)
arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,
       y0=(fromme_data[,"Mean PFOS"]-fromme_data[,"SD PFOS"]),
       y1=(fromme_data[,"Mean PFOS"]+fromme_data[,"SD PFOS"]),code=3,angle=90,length=0.2,col=2)
points(moba_data[,"age"]/12,moba_data[,"Mean PFOS"],col=4,pch=17)
arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,
       y0=(moba_data[,"Mean PFOS"]-moba_data[,"SD PFOS"]),
       y1=(moba_data[,"Mean PFOS"]+moba_data[,"SD PFOS"]),code=3,angle=90,length=0.2,col=4)
points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOS"],col=3,pch=18)
arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,
       y0=(mogensen_data[,"Mean.PFOS"]-mogensen_data[,"SD.PFOS"]),
       y1=(mogensen_data[,"Mean.PFOS"]+mogensen_data[,"SD.PFOS"]),code=3,angle=90,length=0.2,col=3)
legend('topright',legend=c("Fromme, 2010","MOBA Cohort","Mogensen, 2015"),pch=c(16,17,18),col=c(2,4,3))
dev.off()

stop()
##--------------------------------------------------------------------------------------------------------------
## The default model has constant vd. For variable vd, Vdaf_input has to be set to load vd_child.csv in human_adult (human_functions.R) and human_mother_child (sim_functions.R)
source("sim_functions.R")# Change Vdaf_input in line 261
source("human_functions.R")# Change Vdaf_input in line 156

DW=FALSE #With exposure in terms of a constant drinking water level.

## Use reverse dosimetry on maternal levels at delivery to set Fromme exposure
Fromme_PFOA_dose <- child_rev(Cm_Fromme_PFOA,age=0,sex="Male",chem="PFOA",maternal=T,DW=DW,
                              waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
## Run model with Fromme_PFOA_dose.
Fromme_PFOA_data <- human_child(child_age=19/12,sex="Male",chem="PFOA",dose=Fromme_PFOA_dose,DW=DW,
                                waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
## Add a column with C_i normalized by dose (and some scalar)
Fromme_PFOA_data <- data.frame(Fromme_PFOA_data,C_i_n=Fromme_PFOA_data[,"C_i"]/(Fromme_PFOA_dose*dose_scale))

## Calculation for MOBA_PFOA as was performed previously for Fromme_PFOA
MOBA_PFOA_dose <- child_rev(Cm_MOBA_PFOA,age=0,sex="Male",chem="PFOA",maternal=T,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOA_data <- human_child(child_age=36/12,sex="Male",chem="PFOA",dose=MOBA_PFOA_dose,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOA_data <- data.frame(MOBA_PFOA_data,C_i_n=MOBA_PFOA_data[,"C_i"]/(MOBA_PFOA_dose*dose_scale))

## Calculation for Mogensen_PFOA as was performed previously for Fromme_PFOA. Reverse dosimitry was based on measured maternal data. Their estimate of child level (which they present in the paper) was based on multiplying the measured maternal level by a ratio of cord to maternal serum.
Mogensen_PFOA_dose <- child_rev(Cm_mogensen_PFOA,age=0,sex="Male",chem="PFOA",maternal=T,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOA_data <- human_child(child_age=60/12,sex="Male",chem="PFOA",dose=Mogensen_PFOA_dose,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOA_data <- data.frame(Mogensen_PFOA_data,C_i_n=Mogensen_PFOA_data[,"C_i"]/(Mogensen_PFOA_dose*dose_scale))

## Analogous PFOS calculations (with DW=TRUE):
Fromme_PFOS_dose <- child_rev(Cm_Fromme_PFOS,age=0,sex="Male",chem="PFOS",maternal=T,DW=DW,waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
Fromme_PFOS_data <- human_child(child_age=19/12,sex="Male",chem="PFOS",dose=Fromme_PFOS_dose,DW=DW,waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv")
Fromme_PFOS_data <- data.frame(Fromme_PFOS_data,C_i_n=Fromme_PFOS_data[,"C_i"]/(Fromme_PFOS_dose*dose_scale))

MOBA_PFOS_dose <- child_rev(Cm_MOBA_PFOS,age=0,sex="Male",chem="PFOS",maternal=T,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOS_data <- human_child(child_age=36/12,sex="Male",chem="PFOS",dose=MOBA_PFOS_dose,DW=DW,waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv")
MOBA_PFOS_data <- data.frame(MOBA_PFOS_data,C_i_n=MOBA_PFOS_data[,"C_i"]/(MOBA_PFOS_dose*dose_scale))

Mogensen_PFOS_dose <- child_rev(Cm_mogensen_PFOS,age=0,sex="Male",chem="PFOS",maternal=T,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOS_data <- human_child(child_age=60/12,sex="Male",chem="PFOS",dose=Mogensen_PFOS_dose,DW=DW,waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv")
Mogensen_PFOS_data <- data.frame(Mogensen_PFOS_data,C_i_n=Mogensen_PFOS_data[,"C_i"]/(Mogensen_PFOS_dose*dose_scale))

## Plot for PFOA and PFOS that is not normalized by dose. This is used in the Appendix for the variable Vd example
## svg(filename="validation.svg",height=7,width=14)
png(filename="validationVdvar.png",height=7,width=14,units="in",res=600)
par(mfrow=c(1,2),mar=c(4.1,4.1,2.1,2.1),lwd=2)
## Dose unnormalized plot - PFOA
plot(Fromme_PFOA_data[,"years"],Fromme_PFOA_data[,"C_i"],type='l',xlim=c(0,60/12),ylim=c(0,15),xlab="Age (yr)",ylab="PFOA (ng/ml)",col=2)
lines(MOBA_PFOA_data[,"years"],MOBA_PFOA_data[,"C_i"],col=4)
lines(Mogensen_PFOA_data[,"years"],Mogensen_PFOA_data[,"C_i"],col=3,lty=3)
points(fromme_data[,"age"]/12,fromme_data[,"Mean.PFOA"],col=2,pch=16)
arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,
       y0=(fromme_data[,"Mean.PFOA"]-fromme_data[,"SD.PFOA"]),
       y1=(fromme_data[,"Mean.PFOA"]+fromme_data[,"SD.PFOA"]),code=3,angle=90,length=0.2,col=2)
points(moba_data[,"age"]/12,moba_data[,"Mean.PFOA"],col=4,pch=17)
arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,
       y0=(moba_data[,"Mean.PFOA"]-moba_data[,"SD.PFOA"]),
       y1=(moba_data[,"Mean.PFOA"]+moba_data[,"SD.PFOA"]),code=3,angle=90,length=0.2,col=4)
points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOA"],col=3,pch=18)
arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,
       y0=(mogensen_data[,"Mean.PFOA"]-mogensen_data[,"SD.PFOA"]),
       y1=(mogensen_data[,"Mean.PFOA"]+mogensen_data[,"SD.PFOA"]),code=3,angle=90,length=0.2,col=3)
legend('topright',legend=c("Fromme, 2010","MOBA Cohort","Mogensen, 2015"),pch=c(16,17,18),col=c(2,4,3))

## Dose unnormalized plot - PFOS
plot(Fromme_PFOS_data[,"years"],Fromme_PFOS_data[,"C_i"],type='l',xlim=c(0,60/12),ylim=c(0,40),xlab="Age (yr)",ylab="PFOS (ng/ml)",col=2)
lines(MOBA_PFOS_data[,"years"],MOBA_PFOS_data[,"C_i"],col=4)
lines(Mogensen_PFOS_data[,"years"],Mogensen_PFOS_data[,"C_i"],col=3,lty=3)
points(fromme_data[,"age"]/12,fromme_data[,"Mean PFOS"],col=2,pch=16)
arrows(x0=fromme_data[,"age"]/12,x1=fromme_data[,"age"]/12,
       y0=(fromme_data[,"Mean PFOS"]-fromme_data[,"SD PFOS"]),
       y1=(fromme_data[,"Mean PFOS"]+fromme_data[,"SD PFOS"]),code=3,angle=90,length=0.2,col=2)
points(moba_data[,"age"]/12,moba_data[,"Mean PFOS"],col=4,pch=17)
arrows(x0=moba_data[,"age"]/12,x1=moba_data[,"age"]/12,
       y0=(moba_data[,"Mean PFOS"]-moba_data[,"SD PFOS"]),
       y1=(moba_data[,"Mean PFOS"]+moba_data[,"SD PFOS"]),code=3,angle=90,length=0.2,col=4)
points(mogensen_data[,"Age"]/12,mogensen_data[,"Mean.PFOS"],col=3,pch=18)
arrows(x0=mogensen_data[,"Age"]/12,x1=mogensen_data[,"Age"]/12,
       y0=(mogensen_data[,"Mean.PFOS"]-mogensen_data[,"SD.PFOS"]),
       y1=(mogensen_data[,"Mean.PFOS"]+mogensen_data[,"SD.PFOS"]),code=3,angle=90,length=0.2,col=3)
legend('topright',legend=c("Fromme, 2010","MOBA Cohort","Mogensen, 2015"),pch=c(16,17,18),col=c(2,4,3))
dev.off()

## Wrangling of data to plot all 3 studies together:
fromme_data$pre.PFOA <- c(human_child_conc(child_age=fromme_data[1,1]/12,sex="Male",chem="PFOA",dose=Fromme_PFOA_dose,DW=DW,
                                           waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv"),
                          human_child_conc(child_age=fromme_data[2,1]/12,sex="Male",chem="PFOA",dose=Fromme_PFOA_dose,DW=DW,
                                           waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv"))
moba_data$pre.PFOA <- c(human_child_conc(child_age=moba_data[1,1]/12,sex="Male",chem="PFOA",dose=MOBA_PFOA_dose,DW=DW,
                                           waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv"))
mogensen_data$pre.PFOA <- c(human_child_conc(child_age=mogensen_data[1,1]/12,sex="Male",chem="PFOA",dose=Mogensen_PFOA_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"),
                            human_child_conc(child_age=mogensen_data[2,1]/12,sex="Male",chem="PFOA",dose=Mogensen_PFOA_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"),
                            human_child_conc(child_age=mogensen_data[3,1]/12,sex="Male",chem="PFOA",dose=Mogensen_PFOA_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"))

fromme_data$pre.PFOS <- c(human_child_conc(child_age=fromme_data[1,1]/12,sex="Male",chem="PFOS",dose=Fromme_PFOS_dose,DW=DW,
                                           waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv"),
                          human_child_conc(child_age=fromme_data[2,1]/12,sex="Male",chem="PFOS",dose=Fromme_PFOS_dose,DW=DW,
                                           waterin_m="water_lact_fromme.csv",waterin_i="water_bf_fromme.csv",milkin="r_milk_bw_fromme.csv"))
moba_data$pre.PFOS <- c(human_child_conc(child_age=moba_data[1,1]/12,sex="Male",chem="PFOS",dose=MOBA_PFOS_dose,DW=DW,
                                           waterin_m="water_lact.csv",waterin_i="water_bf.csv",milkin="r_milk_bw.csv"))
mogensen_data$pre.PFOS <- c(human_child_conc(child_age=mogensen_data[1,1]/12,sex="Male",chem="PFOS",dose=Mogensen_PFOS_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"),
                            human_child_conc(child_age=mogensen_data[2,1]/12,sex="Male",chem="PFOS",dose=Mogensen_PFOS_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"),
                            human_child_conc(child_age=mogensen_data[3,1]/12,sex="Male",chem="PFOS",dose=Mogensen_PFOS_dose,DW=DW,
                                             waterin_m="water_lact_mogensen.csv",waterin_i="water_bf_mogensen.csv",milkin="r_milk_bw_mogensen.csv"))

pl_data_pfoa <- data.frame(observed=c(fromme_data[,2],moba_data[,2],mogensen_data[,3]),
                           observed_uci=c(fromme_data[,2]+fromme_data[,3],moba_data[,2]+moba_data[,3],mogensen_data[,3]+mogensen_data[,5]),
                           observed_lci=c(fromme_data[,2]-fromme_data[,3],moba_data[,2]-moba_data[,3],mogensen_data[,3]-mogensen_data[,5]),
                           predicted=c(fromme_data[,6],moba_data[,6],mogensen_data[,6]))
pl_data_pfos <- data.frame(observed=c(fromme_data[,4],moba_data[,4],mogensen_data[,2]),
                           observed_uci=c(fromme_data[,4]+fromme_data[,5],moba_data[,4]+moba_data[,5],mogensen_data[,2]+mogensen_data[,4]),
                           observed_lci=c(fromme_data[,4]-fromme_data[,5],moba_data[,4]-moba_data[,5],mogensen_data[,2]-mogensen_data[,4]),
                           predicted=c(fromme_data[,7],moba_data[,7],mogensen_data[,7]))

## Store a vector of percent differences for pfoa and pfos with variable and constant vd
error_pfoa_vdvar <- (pl_data_pfoa$predicted-pl_data_pfoa$observed)/pl_data_pfoa$observed
error_pfos_vdvar <- (pl_data_pfos$predicted-pl_data_pfos$observed)/pl_data_pfos$observed

## mean percent different and mean absolute percent difference for observed vs. predicted data.
perc_error <- rbind(c(mean(error_pfoa_vdvar),mean(error_pfos_vdvar),sum(c(mean(error_pfoa_vdvar),mean(error_pfos_vdvar)))),
                    c(mean(error_pfoa_vdcon),mean(error_pfos_vdcon),sum(c(mean(error_pfoa_vdcon),mean(error_pfos_vdcon)))))

abs_error <- rbind(c(mean(abs(error_pfoa_vdvar)),mean(abs(error_pfos_vdvar)),sum(c(mean(abs(error_pfoa_vdvar)),mean(abs(error_pfos_vdvar))))),
                   c(mean(abs(error_pfoa_vdcon)),mean(abs(error_pfos_vdcon)),sum(c(mean(abs(error_pfoa_vdcon)),mean(abs(error_pfos_vdcon))))))

write.csv(perc_error,"perc_error.csv")
write.csv(abs_error,"abs_error.csv")

