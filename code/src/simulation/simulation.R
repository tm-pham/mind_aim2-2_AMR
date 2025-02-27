# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Simulation study for multinomial logistic regression
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Organism: S. aureus
# ---------------------------------------------------------------------------- #
# Simulate data for multinomial logistic regression 
# Facility-level antibiotic use data is based on observed distributions
# A short description of the simulation study is in the Readme file for 
# each simulated dataset. 
# ============================================================================ #
# Load packages
library(dplyr)

# Load functions
source("mind_aim2-2_function_sim_multinom_regr.R")

# Load antibiotic use data 
orderly2::orderly_shared_resource("mind_aim2-2_abx_params_lognormal_FQL_ASBL_MACR_2007_2008.RData")
orderly2::orderly_shared_resource("plotting_template.R")
orderly2::orderly_shared_resource("mind_aim2-2_function_mblogit_sim_coef_plot.R")

################################################################################
# Set global variables
orderly_parameters(p = 1, 
                   org_short = "SA", 
                   data_id = NULL, 
                   n_facilities = 70, 
                   n_obs_per_facility = 150, 
                   gamma = NULL, 
                   folder_path = NULL)

################################################################################
# Parameters for simulation
params <- list(sim_id = data_id, 
               org_short = org_short, 
               time_period = 2007:2008, 
               n_facilities = n_facilities,                         # Number of facilities 
               n_obs_per_facility = n_obs_per_facility,                  # Number of observations by facility
               beta0 = c(0.05, -2, -3, -2, -3, -2, -1.5), # Intercepts corresponding to each antibiotic resistance pattern
               beta = c(0.0001, 0.0000001, 0.0000001),    # Coefficients for time (year, day, month)
               gamma = rep(gamma, length(abx_params)),            # Coefficients for facility-level previous abx use
               sigma = rep(0, 6),                         # Coefficients for facility-level covariates (census region, facility rurality)
               rho = rep(0, 6),                           # Coefficients for patient-level covariates 
               random_effect_sd = 0.0001,                 # Random effects parameter 
               antibiograms = c("R-R-R", "R-R-S", "R-S-R", "S-R-R", "S-S-R"), # Antibiotic resistance patterns, order should coincide with order of abx_distr
               reference = "S-S-S"
)

################################################################################
# Generate data
cat("This script runs a multinomial logistic regression simulation based on the following data set: ", data_id, "\n")
sim_data <- generate.data(params, 
                          seed = p,
                          abx_params)
