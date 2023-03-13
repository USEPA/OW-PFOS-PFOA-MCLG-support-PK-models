##------------------------------------------------------------------------------
## LSA.R - File to perform local sensitivity analysis on selected parameters
##
##------------------------------------------------------------------------------

## Load required functions and libraries.
source("compile_pfas_tk.R")
source("sim_functions.R")
source("human_functions.R")

## A version of human_child that runs the model with 2 parameter values (for the parameter in var_name), the original value and a value 1% greater. Then the sensitivity is calculated as the relative change in serum concentration divided by the 1% change in the parameter.
child_LSA <- function(var_name,child_age=3,sex=NULL,chem=NULL,dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv",maternal=F){
    if(is.null(sex)){
        print("Sex of child not specified.")
        return(NULL)
    }
    if (chem=="PFOA"){
        Vd = 0.170
        t12 = 2.7
        Cl = (Vd*1000)*log(2)/(t12*365)
        ## p_PFOA=c("CL"=0.085/1000,"Vd"=0.170, "CLr"=0.085/1000, "P_milk"=0.0577, "r_f_m"=0.783, "pow_milk"=1) # PFOA parameters
        p_PFOA=c("CL"=0.120/1000,"Vd"=0.170, "CLr"=0.120/1000, "P_milk"=0.049, "r_f_m"=0.83, "pow_milk"=1) # PFOA parameters
        p=c(p_PFOA,"d_m"=dose)# 1 ug/kg/d
    } else if (chem=="PFOS"){
        Vd = 0.230
        t12 = 3.4
        Cl = (Vd*1000)*log(2)/(t12*365)
        ## p_PFOS=c("CL"=0.081/1000,"Vd"=0.230, "CLr"=0.081/1000, "P_milk"=0.0138, "r_f_m"=0.454, "pow_milk"=1) # PFOS parameters
        p_PFOS=c("CL"=0.128/1000,"Vd"=0.230, "CLr"=0.128/1000, "P_milk"=0.0160, "r_f_m"=0.40, "pow_milk"=1) # PFOS parameters
        p=c(p_PFOS,"d_m"=dose)# 1 ug/kg/d
    } else {
        print("Specify PFOA or PFOS.")
        return(NULL)
    }
    if(var_name%in%names(p)){
        p_var <- p
        p_var[var_name] <- p[var_name]*1.01
    } else if (var_name=="t12"){
        p_var <- p
        t12_var <- t12*1.01
        Cl_var = (Vd*1000)*log(2)/(t12_var*365)
        p_var["CL"] = Cl_var/1000
        p_var["CLr"] = Cl_var/1000
    } else if(var_name=="Vd1"){# Vd and effect of Vd on clearance
        p_var <- p
        Vd_var <- p["Vd"]*1.01
        Cl_var <- Vd_var*log(2)/(t12*365)
        p_var["Vd"] = Vd_var
        p_var["CL"] = Cl_var
        p_var["CLr"] = Cl_var
    } else if(var_name=="Vd2"){# Effect of Vd with constant excretion
        p_var <- p
        Vd_var <- p["Vd"]*1.01
        ## Cl_var <- Vd_var*log(2)/(t12*365)
        p_var["Vd"] = Vd_var
        ## p_var["CL"] = Cl_var
        ## p_var["CLr"] = Cl_var
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
    if(maternal){
        out_C <- tail(out[,"C_m"],1)
    } else {
        out_C <- tail(out[,"C_i"],1)
    }
    if(var_name%in%names(p)|var_name=="t12"|var_name=="Vd1"|var_name=="Vd2"){
        out_var <- human_mother_child(age_end=run_child_age+25,p=p_var,sex=sex,mName="pfas_tk",milkin=milkin,waterin_m=waterin_m,waterin_i=waterin_i,DW=DW)

    } else if(var_name=="milkin"){
        milkin= paste(substring(milkin,1,nchar(milkin)-4),"_lsa.csv",sep='')
        out_var <- human_mother_child(age_end=run_child_age+25,p=p,sex=sex,mName="pfas_tk",milkin=milkin,waterin_m=waterin_m,waterin_i=waterin_i,DW=DW)
    }
    out_var <- cbind(out_var,(out_var[,"time"]/365)-25)
    colnames(out_var)[ncol(out_var)] <- "years"
    if(child_age<=0){
        out_var <- out_var[out_var[,"years"]<=child_age,]
    }
    if(maternal){
        out_C_var <- tail(out_var[,"C_m"],1)
    } else {
        out_C_var <- tail(out_var[,"C_i"],1)
    }
    sens <- ((out_C_var-out_C)/out_C)/0.01
    return(sens)
}

### Code below calculated sensitivity for the parameters in par_list for a variety of sex, age, and chemical categories.
par_list <- c("t12","Vd1","P_milk","r_f_m","milkin")
pfoa_male_lsa1 <- sapply(par_list,child_LSA,child_age=1,sex="Male",chem="PFOA",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")
pfos_male_lsa1 <- sapply(par_list,child_LSA,child_age=1,sex="Male",chem="PFOS",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")
pfoa_female_lsa1 <- sapply(par_list,child_LSA,child_age=1,sex="Female",chem="PFOA",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")
pfos_female_lsa1 <- sapply(par_list,child_LSA,child_age=1,sex="Female",chem="PFOS",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")

pfoa_male_lsa5 <- sapply(par_list,child_LSA,child_age=5,sex="Male",chem="PFOA",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")
pfos_male_lsa5 <- sapply(par_list,child_LSA,child_age=5,sex="Male",chem="PFOS",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")
pfoa_female_lsa5 <- sapply(par_list,child_LSA,child_age=5,sex="Female",chem="PFOA",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")
pfos_female_lsa5 <- sapply(par_list,child_LSA,child_age=5,sex="Female",chem="PFOS",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")

pfoa_male_lsa0 <- sapply(par_list,child_LSA,child_age=0,sex="Male",chem="PFOA",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")
pfos_male_lsa0 <- sapply(par_list,child_LSA,child_age=0,sex="Male",chem="PFOS",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")
pfoa_female_lsa0 <- sapply(par_list,child_LSA,child_age=0,sex="Female",chem="PFOA",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")
pfos_female_lsa0 <- sapply(par_list,child_LSA,child_age=0,sex="Female",chem="PFOS",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv")

pfoa_male_lsaM <- sapply(par_list,child_LSA,child_age=0,sex="Male",chem="PFOA",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv",maternal=T)
pfos_male_lsaM <- sapply(par_list,child_LSA,child_age=0,sex="Male",chem="PFOS",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv",maternal=T)
pfoa_female_lsaM <- sapply(par_list,child_LSA,child_age=0,sex="Female",chem="PFOA",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv",maternal=T)
pfos_female_lsaM <- sapply(par_list,child_LSA,child_age=0,sex="Female",chem="PFOS",dose=0.001,DW=FALSE,waterin_m="water_95_lact.csv",waterin_i="water_95_bf.csv",milkin="r_milk_bw_95.csv",maternal=T)

## Fancier names for par_list for plotting
long_names <- c("Half-life", "Volume of Distribution", "Milk:Maternal Serum Ratio","Cord Blood:Maternal Serum Ratio", "Breastmilk Intake per kg BW")

##  Save svg files of plots of sensitivity coefficients.
svg("LSA1.svg",width=10,height=6.5)
par(las=2,cex=1.4,cex.main=2,lwd=3)
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_male_lsa1,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2,2), main="Male PFOA")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_male_lsa1,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2,2), main="Male PFOS")
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
axis(1,cex.axis=1.5,las=1)
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_female_lsa1,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2,2), main="Female PFOA")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_female_lsa1,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2,2), main="Female PFOS")
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
dev.off()

svg("LSA5.svg",width=10,height=6.5)
par(las=2,cex=1.4,cex.main=2,lwd=3)
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_male_lsa5,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2,2), main="Male PFOA")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_male_lsa5,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2,2), main="Male PFOS")
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
axis(1,cex.axis=1.5,las=1)
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_female_lsa5,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2,2), main="Female PFOA")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_female_lsa5,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2,2), main="Female PFOS")
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
dev.off()

svg("LSA0.svg",width=10,height=6.5)
par(las=2,cex=1.4,cex.main=2,lwd=3)
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_male_lsa0,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2,2), main="Male PFOA")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_male_lsa0,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2,2), main="Male PFOS")
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
axis(1,cex.axis=1.5,las=1)
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_female_lsa0,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2,2), main="Female PFOA")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_female_lsa0,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2,2), main="Female PFOS")
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
dev.off()

svg("LSAM.svg",width=10,height=6.5)
par(las=2,cex=1.4,cex.main=2,lwd=3)
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_male_lsaM,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2,2), main="Male PFOA")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_male_lsaM,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2,2), main="Male PFOS")
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
axis(1,cex.axis=1.5,las=1)
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_female_lsaM,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2,2), main="Female PFOA")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_female_lsaM,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2,2), main="Female PFOS")
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
dev.off()

svg("pfoa_male_LSA.svg",width=10,height=6.5)
par(las=2,cex=1.4,cex.main=2,lwd=3)
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_male_lsaM,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2.5,2.5), main="Maternal at Delivery")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfoa_male_lsa0,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2.5,2.5), main="Cord Blood")
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
axis(1,cex.axis=1.5,las=1)
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_male_lsa1,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2.5,2.5), main="Child at 1 yr")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfoa_male_lsa5,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2.5,2.5), main="Child at 5 yr")
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
dev.off()

svg("pfoa_female_LSA.svg",width=10,height=6.5)
par(las=2,cex=1.4,cex.main=2,lwd=3)
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_female_lsaM,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2.5,2.5), main="Maternal at Delivery")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfoa_female_lsa0,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2.5,2.5), main="Cord Blood")
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
axis(1,cex.axis=1.5,las=1)
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfoa_female_lsa1,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2.5,2.5), main="Child at 1 yr")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfoa_female_lsa5,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2.5,2.5), main="Child at 5 yr")
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
dev.off()

svg("pfos_male_LSA.svg",width=10,height=6.5)
par(las=2,cex=1.4,cex.main=2,lwd=3)
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfos_male_lsaM,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2.5,2.5), main="Maternal at Delivery")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_male_lsa0,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2.5,2.5), main="Cord Blood")
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
axis(1,cex.axis=1.5,las=1)
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfos_male_lsa1,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2.5,2.5), main="Child at 1 yr")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_male_lsa5,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2.5,2.5), main="Child at 5 yr")
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
dev.off()

svg("pfos_female_LSA.svg",width=10,height=6.5)
par(las=2,cex=1.4,cex.main=2,lwd=3)
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfos_female_lsaM,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2.5,2.5), main="Maternal at Delivery")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_female_lsa0,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2.5,2.5), main="Cord Blood")
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
axis(1,cex.axis=1.5,las=1)
par(mar=c(2.1,0.4,2.1,14.1))
plot(0,ylim=c(1,14),type='n',axes=F)
par(mar=c(3.1,0.4,2.1,0.4))
a <- barplot(pfos_female_lsa1,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
             xlim=c(-2.5,2.5), main="Child at 1 yr")
axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
barplot(pfos_female_lsa5,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
        xlim=c(-2.5,2.5), main="Child at 5 yr")
axis(1,cex.axis=1.5,las=1)
abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
dev.off()


## par(las=2,cex=1.4,cex.main=2,lwd=3)
## layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
## par(mar=c(2.1,0.4,2.1,14.1))
## plot(0,ylim=c(1,14),type='n',axes=F)
## par(mar=c(3.1,0.4,2.1,0.4))
## a <- barplot(pfoa_male_lsa,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
##              xlim=c(-2,2), main="Male PFOA")
## axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
## axis(1,cex.axis=1.5,las=1)
## abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
## barplot(pfos_male_lsa,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
##         xlim=c(-2,2), main="Male PFOS")
## abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
## axis(1,cex.axis=1.5,las=1)
## par(mar=c(2.1,0.4,2.1,14.1))
## plot(0,ylim=c(1,14),type='n',axes=F)
## par(mar=c(3.1,0.4,2.1,0.4))
## a <- barplot(pfoa_female_lsa,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
##              xlim=c(-2,2), main="Female PFOA")
## axis(2,a,labels=long_names,lty=0,line=-1.5,cex.axis=1.5)
## axis(1,cex.axis=1.5,las=1)
## abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
## barplot(pfos_female_lsa,col=0,horiz=T,yaxt='n',xaxt='n',space=0.5,
##         xlim=c(-2,2), main="Female PFOS")
## axis(1,cex.axis=1.5,las=1)
## abline(v=c(0,1,-1),lwd=1,lty=c(1,2,2))
