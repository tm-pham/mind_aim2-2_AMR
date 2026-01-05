# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Simulation study for multinomial logistic regression
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# ---------------------------------------------------------------------------- #
# Run file for simulation study
# 1. Generate synthetic data for multinomial logistic regression
# 2. Save simulated data as .RData file
# ============================================================================ #
# Load function to simulate data for multinomial logistic regression (generate.data)
source(here::here("simulation","00_config_simulation.R"))
source(here::here("simulation","01_function_sim_multinom_regr.R"))
# Load parameters for facility-level antibiotic use distributions
load(here::here("data","mind_aim2-2_abx_params.RData"))

# ------------------------------------------------------------------------------
# Parameters for simulation
# ------------------------------------------------------------------------------
params <- list(sim_id = "sim_multinom_regr",              # Simulation ID
               org_short = "SA", 
               time_period = 2007:2008,                   # Time period for simulation
               n_facilities = 70,               # Number of facilities 
               n_obs_per_facility = 150,   # Number of observations by facility
               beta0 = rep(-0.5, 4),          # Intercepts corresponding to each antibiotic resistance pattern
               beta = c(0.0001, 0.0000001, 0.0000001),    # Coefficients for time (year, day, month)
               gamma = rep(0.001, length(abx_params)),    # Coefficients for facility-level previous abx use
               sigma = rep(0, 6),                         # Coefficients for facility-level covariates (census region, facility rurality)
               rho = rep(0, 6),                           # Coefficients for patient-level covariates 
               random_effect_sd = 0.0001,                 # Random effects parameter 
               antibiograms = c("R-R-R", "R-R-S", "R-S-R", "S-R-R", "S-S-R"), # Antibiotic resistance patterns, order should coincide with order of abx_distr
               reference = "S-S-S"
)

# ------------------------------------------------------------------------------
# Generate data
# ------------------------------------------------------------------------------
cat("This script runs a multinomial logistic regression simulation with the following parameters:", params, "\n")
sim_data <- generate.data(params, 
                          seed = p,
                          abx_params)
