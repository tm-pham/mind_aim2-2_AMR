# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Simulation study for multinomial logistic regression
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# ---------------------------------------------------------------------------- #
# Run file for simulation study
# 1. Generate synthetic data for multinomial logistic regression
# 2. Save simulated data as .RData file
# In this simulation, we use abx use distributions where we keep the mean fixed 
# and vary the variance. The respective distributions are generated from 
# observed fluoroquinolone antibiotic prescribing data. 
# ============================================================================ #
# Load function to simulate data for multinomial logistic regression (generate.data)
source(here::here("code","simulation","00_config_simulation.R"))
source(here::here("code", "simulation","01_function_sim_multinom_regr.R"))

# ------------------------------------------------------------------------------
# Generate data
# ------------------------------------------------------------------------------
cat("This script runs a multinomial logistic regression simulation.\n")
cat("Generating", n_sim, "synthetic antimicrobial use and resistance data.\n")

# Parent folder for this simulation
data_id = paste0("sim", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
sim_ids <- rep(NA, n_sim)

for(i in 1:length(gamma_distr)){
  sim_ids[i] <- i
  params <- list(data_id = data_id, 
                 sim_id = sim_ids[i],              # Simulation ID
                 org_short = "SA", 
                 time_period = 2007:2008,                   # Time period for simulation
                 n_facilities = 70,               # Number of facilities 
                 n_obs_per_facility = 150,   # Number of observations by facility
                 beta0 = rep(-0.5, 5),          # Intercepts corresponding to each antibiotic resistance pattern
                 beta = c(0.0001, 0.0000001, 0.0000001),    # Coefficients for time (year, day, month)
                 gamma = gamma_distr[i],    # Coefficients for facility-level previous abx use
                 sigma = rep(0, 6),                         # Coefficients for facility-level covariates (census region, facility rurality)
                 rho = rep(0, 6),                           # Coefficients for patient-level covariates 
                 random_effect_sd = 0.0001,                 # Random effects parameter 
                 antibiograms = c("R-R-R", "R-R-S", "R-S-R", "S-R-R", "S-S-R"), # Antibiotic resistance patterns, order should coincide with order of abx_distr
                 reference = "S-S-S"
  )
  
  sim_data <- generate.data(params, 
                            seed = i,
                            abx_params, 
                            prefix = PREFIX)
}

# Save the simulation IDs
save(data_id, sim_ids, file = here::here("results", "simulation", data_id, "sim_ids.RData"))




