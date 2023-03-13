## ## Code to compile RMCSim package locally
## Sys.setenv(PATH = paste("C:/Rtools40/usr/bin",Sys.getenv("PATH"), sep=";"))
## library(devtools)
## devtools::install_local("C:/Users/MDZIERLE/Documents/RMCSim/",force=T)

## Script to compile the pfas_k.model using local RMCSim package
library(RMCSim)
Sys.setenv(PATH = paste("C:/Rtools40/usr/bin",Sys.getenv("PATH"), sep=";"))
RMCSim::compile_model("pfas_tk")
RMCSim::load_model("pfas_tk")

## Load deSolve.
library(deSolve)
