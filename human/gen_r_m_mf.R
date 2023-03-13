##------------------------------------------------------------------------------
## gen_r_m_mf.R
## Generates an interpolated value of r_m_mf to reproduce the acslX model
##
##------------------------------------------------------------------------------

interp <- function(time,times,var){
    if(all(time<times)){
        return(var[1])
    } else if(all(time>=times)){
        return(var[length(var)])
    }
    index <- tail(which(time>=times),1)
    val <- var[index]*(times[index+1]-time)/(times[index+1]-times[index]) + var[index+1]*(1-(times[index+1]-time)/(times[index+1]-times[index]))
    return(val)
}

M_df = read.csv("./Verner_2015_PFOA_mass.csv", header=FALSE,# Maternal and child BW is not chem dependent.
                   col.names=c("years", "M_m","M_i"),
                   check.names=FALSE)

r_f_m <- 0.783
r_m_mf <- ifelse(M_df[,1]<=25,(M_df[,2]+M_df[,3])/(M_df[,2]+r_f_m*M_df[,3]),1)

write.csv(cbind(M_df[,1],r_m_mf),file="./r_m_mf_pfoa_v.csv")

r_f_m <- 0.454
r_m_mf <- ifelse(M_df[,1]<=25,(M_df[,2]+M_df[,3])/(M_df[,2]+r_f_m*M_df[,3]),1)

write.csv(cbind(M_df[,1],r_m_mf),file="./r_m_mf_pfos_v.csv")

## M_m_df1 = read.csv("human_female_mass_0_to_25_preg.csv", header=FALSE,
##                    col.names=c("Age", "Mass (kg)"),
##                    check.names=FALSE)
## M_m_df1[,1] = M_m_df1[,1]*365 # Convert age from years to days

## M_m_df2 =read.csv("human_female_mass_25_to_35_postnatal.csv",
##                   header=FALSE, col.names=c("Age", "Mass (kg)"),
##                   check.names=FALSE)
## M_m_df2[,1] = M_m_df2[,1]*365 # Convert age from years to days

