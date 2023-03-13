#------------------------------------------------------------------------------
# RMCSim.R
#
# This file contains R functions for compiling, loading, and running an ODE
# model encoded in the GNU MCSim model specification language. The "mod.exe"
# utility must be available in the user's PATH.
#
# Author: Dustin Kapraun, U.S. EPA, November 2018.
#------------------------------------------------------------------------------

# Load deSolve.
library(deSolve)

compile_model <- function(mName) {
    # This function translates a model that has been defined in an MCSim model
    # (".model") file into the C language (i.e., a ".c" file). It then compiles the
    # model to create an object code (".o") file and a dynamic linked library
    # (".dll") file, as well as an R script ("_inits.R") containing several R
    # functions that can be used for initializing model states and parameters.
    #
    # Inputs:
    #   mName: String containing the name of the MCSim model. Exclude the file name
    #   suffix ".model". If the function is called from a working directory other
    #   than the one containing the ".model" file, the full path of the ".model"
    #   file should be provided.
    
    # Construct names of required files and objects from mName.
    model_file = paste(mName, ".model", sep="")
    c_file = paste(mName, "_model.c", sep="")
    dll_name = paste(mName, "_model", sep="")
    dll_file = paste(dll_name, .Platform$dynlib.ext, sep="")
    
    # Unload DLL if it has been loaded.
    if (is.loaded("derivs", PACKAGE=dll_name)) {
        dyn.unload(dll_file)
    }
    
    # Create a C model file (ending with ".c") and an R parameter
    # initialization file (ending with "_inits.R") from the GNU MCSim model
    # definition file (ending with ".model"). Using the "-R" option generates
    # code compatible with functions in the R deSolve package.
    system(paste("mod -R ", model_file, " ", c_file, sep = ""))
    
    # Compile the C model to obtain "mName_model.o" and "mName_model.dll".
    system(paste("R CMD SHLIB ", c_file, sep = ""))
}


load_model <- function(mName) {
    # This function loads a model that is defined using DLL and "_inits.R" files
    # based on translation of an MCSim model (".model") file. The R functions
    # "initParms" and "initStates" and the R vector "Outputs" are defined and
    # assigned to the "global" environment.
    #
    # Inputs:
    #   mName: String containing the name of the MCSim model. Exclude the file name
    #   suffix ".model". If the function is called from a working directory other
    #   than the one containing the ".model" file, the full path of the ".model"
    #   file should be provided.
    
    # Construct names of required files and objects from mName.
    dll_name = paste(mName, "_model", sep="")
    dll_file = paste(dll_name, .Platform$dynlib.ext, sep="")
    inits_file = paste(dll_name, "_inits.R", sep="")
    
    # Load the compiled model (DLL).
    dyn.load(dll_file)
    
    # Run script that defines initialization functions.
    source(inits_file)
    
    # Assign initialization functions and list of output variable names to the
    # "global" environment.
    assign("initParms", initParms, envir=.GlobalEnv)
    assign("initStates", initStates, envir=.GlobalEnv)
    assign("Outputs", Outputs, envir=.GlobalEnv)
}


run_model <- function(mName, times, Y0=NULL, parms=NULL, rtol=1e-6, atol=1e-6,
                      forcing=NULL, fcontrol=NULL, event_list=NULL, method="lsoda") {
    # This function runs a model that is defined using DLL and "_inits.R" files
    # based on translation of an MCSim model (".model") file. The model must first
    # be compiled (using the function "compile_model") and loaded into the R
    # workspace (using the function "load_model").
    #
    # Inputs:
    #   mName: String containing the name of the MCSim model. Exclude the file name
    #   suffix ".model". If the function is called from a working directory other
    #   than the one containing the ".model" file, the full path of the ".model"
    #   file should be provided.
    #
    #   times: Vector containing the times for which the model simulation results
    #   are desired. This vector should contain a minimum of two values (an initial
    #   time and an end time.
    #
    #   Y0: Vector containing the initial values of the state variables in the
    #   model. If this argument is NULL, the default initial values of the state
    #   variables (as returned by the function "initStates") will be used.
    #
    #   parms: Vector containing the values of the model parameters. If this
    #   argument is NULL, the default values of the parameters (as returned by the
    #   function "initParms") will be used.
    #
    #   r_tol: Relative error tolerance. See the documentation for the R package
    #   "deSolve".
    #
    #   a_tol: Absolute error tolerance. See the documentation for the R package
    #   "deSolve".
    #
    #   forcing: List containing values of input variables to be passed to the ODE
    #   solver. See the documentation for the R package "deSolve".
    #
    #   fcontrol: List containing arguments that describe how interpolations of the
    #   forcing functions should be performed. See the documentation for the R
    #   package "deSolve".
    #
    #   event_list: List containing events to be passed to the ODE solver. See
    #   the documentation for the R package deSolve.
    #
    # Returns:
    #   out: Dataframe containing values of all state variables (named in the
    #   argument "Y0") and all output variables (named in the global environment
    #   variable "Outputs") at all times (defined in the argument "times").
    
    # Construct DLL name from mName.
    dll_name = paste(mName, "_model", sep="")
    
    # If parameter values are not provided, use default values.
    if (is.null(parms)) {
        parms = initParms()
    }
    
    # If initial values for state variables are not provided, use default
    # values.
    if (is.null(Y0)) {
        Y0 = initStates(parms)
    }
    
    # Solve the ODE system using the "ode" function from the package "deSolve".
    out = ode(Y0, times, func="derivs", parms=parms, rtol=rtol, atol=atol,
              dllname=dll_name, initforc="initforc", forcing=forcing,
              fcontrol=fcontrol, initfunc="initmod", nout=length(Outputs),
              outnames=Outputs, events=event_list, method=method)

    # Return the simulation output.
    return(out)
}