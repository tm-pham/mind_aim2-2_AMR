# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Simulation study for multinomial logistic regression for scenario 1
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ============================================================================ #
library(dplyr)

################################################################################
# Generate distributions for regression coefficients for abx use (gamma)
orderly_parameters(mean = 0.0025, sd = 0.1, 
                   n_sim = 5, 
                   org_short = "SA", data_id = "SA_20072021_4", 
                   n_facilities = 70, n_obs_per_facility = 150)
gamma_distr <- rlnorm(n_sim, mean = log(mean), sd = sd)

# Make a new orderly task
sim_ids <- vector("character", length(gamma_distr))
for(i in 1:length(gamma_distr)){
  sim_ids[i] <- orderly2::orderly_run("simulation", 
                                      location = "./simulation/",
                                      parameters = list(p = i, 
                                                        org_short = org_short, 
                                                        data_id = data_id,
                                                        n_facilities = n_facilities, 
                                                        n_obs_per_facility = n_obs_per_facility, 
                                                        gamma = gamma_distr[i], 
                                                        folder_path = "/Users/tm-pham/academia/hsph/mind/publications/aim2-2/data_output/"))
}

# Save the simulation IDs
save(data_id, sim_ids, file = "sim_ids.RData")