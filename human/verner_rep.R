##------------------------------------------------------------------------------
## verner_rep.R
## Code to replicate results of Verner model, which were obtained from running
## the original acslX model and saving the output in csv form.
##
##------------------------------------------------------------------------------

## Load required functions and libraries.
source("compile_pfas_tk.R")
source("sim_functions.R")
source("human_functions.R")

## Run the updated verner model
female_pfoa_child <- human_child(child_age=5,sex="Female",chem="PFOA")
female_pfos_child <- human_child(child_age=5,sex="Female",chem="PFOS")

## Load verner-specific functions & parameters and run the model with original parameter
source("human_functions_vbw.R")
female_pfoa_child_v <- human_child_v(child_age=10,sex="Female",chem="PFOA")
female_pfos_child_v <- human_child_v(child_age=10,sex="Female",chem="PFOS")

## Load saved output from running original acslX model and convert to mg/L
acsl_pfoa <- read.csv(file="./Verner_2015_PFOA.csv",header=F)
colnames(acsl_pfoa) <- c("years","C_m","C_i")

acsl_pfoa[,2] <- acsl_pfoa[,2]/1000 ## Converts ug/L in plasma to mg/L
acsl_pfoa[,3] <- acsl_pfoa[,3]/1000 ## Converts ug/L in plasma to mg/L

acsl_pfos <- read.csv(file="./Verner_2015_PFOS.csv",header=F)
colnames(acsl_pfos) <- c("years","C_m","C_i")

acsl_pfos[,2] <- acsl_pfos[,2]/1000 ## Converts ug/L in plasma to mg/L
acsl_pfos[,3] <- acsl_pfos[,3]/1000 ## Converts ug/L in plasma to mg/L

## Generate figure comparing outputs of updated Verner model, R Verner model, and acslX Verner model.
## svg(filename="pfas_tk_verner_rep.svg",height=7,width=14)
png(filename="pfas_tk_verner_rep.png",height=7,width=14,units="in",res=600)
par(mfrow=c(1,2),lwd=2,mar=c(4.1,4.1,2.1,2.1))
plot(female_pfoa_child[,"years"],female_pfoa_child[,"C_m"],type='l',
     xlab="Age (yr)",ylab="PFOA Plasma Concentration (mg/L)",xaxs='i',
     ylim=c(0,max(acsl_pfoa[,"C_i"])),xlim=c(-1,1.5))
lines(female_pfoa_child[,"years"],female_pfoa_child[,"C_i"],lty=2)
lines(female_pfoa_child_v[,"years"],female_pfoa_child_v[,"C_m"],col=2)
lines(female_pfoa_child_v[,"years"],female_pfoa_child_v[,"C_i"],col=2,lty=2)
lines(acsl_pfoa[,"years"]-25,acsl_pfoa[,"C_m"],col=4)
lines(acsl_pfoa[,"years"]-25,acsl_pfoa[,"C_i"],lty=2,col=4)

plot(female_pfos_child[,"years"],female_pfos_child[,"C_m"],type='l',
     xlab="Age (yr)",ylab="PFOS Plasma Concentration (mg/L)",xaxs='i',
     ylim=c(0,max(acsl_pfos[,"C_i"])),xlim=c(-1,1.5))
lines(female_pfos_child[,"years"],female_pfos_child[,"C_i"],lty=2)
lines(female_pfos_child_v[,"years"],female_pfos_child_v[,"C_m"],col=2)
lines(female_pfos_child_v[,"years"],female_pfos_child_v[,"C_i"],col=2,lty=2)
lines(acsl_pfos[,"years"]-25,acsl_pfos[,"C_m"],col=4)
lines(acsl_pfos[,"years"]-25,acsl_pfos[,"C_i"],lty=2,col=4)
legend(x="bottomright",legend=c("Maternal PFAS_TK","Child PFAS_TK","Maternal Verner","Child Verner",
                                "Maternal Rep.","Child Rep."),
       col=c(1,1,2,2,4,4),lty=c(1,2,1,2,1,2))
dev.off()
