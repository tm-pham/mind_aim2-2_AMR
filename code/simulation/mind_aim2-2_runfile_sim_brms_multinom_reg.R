# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Simulation study for multinomial logistic regression (brms package)
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Organism: S. aureus
# ---------------------------------------------------------------------------- #
# 
# ============================================================================ #
setwd("/Users/tm-pham/academia/hsph/mind/")
library(ggplot2)
library(dplyr)
library(brms)

# Load functions
source("code/plotting_template.R")
source("code/aim2-2/simulation/mind_aim2-2_function_sim_multinom_regr.R")
source("code/aim2-2/simulation/mind_aim2-2_function_mblogit_sim_coef_plot.R")
source("data/SA/aim2-2/SA_20072021_2/mind_aim2-2_SA_20072021_2_params.R") # df_input_comb


# Load dataset
load("data_output/simReg/SA/SA_20072021_2_1/mind_aim2-2_simReg_SA_20072021_2_1.RData")
# load("data_output/simReg/SA_20072021_2_2/mind_aim2-2_simReg_SA_20072021_2_2.RData")
sim_data <- simReg$data

# Set priors 
priors <- c(
  set_prior("normal(0, 9/4)", class = "b", coef = "beta_lactams", dpar = paste0("mu", c("RRRR", "RRSS", "RSSR", "SSRR", "SRRR", "RRSR", "SSSR"))),
  set_prior("normal(0, 9/4)", class = "b", coef = "fluoroquinolones", dpar = paste0("mu", c("RRRR", "RRSS", "RSSR", "SSRR", "SRRR", "RRSR", "SSSR"))),
  set_prior("normal(0, 9/4)", class = "b", coef = "lincosamides", dpar = paste0("mu", c("RRRR", "RRSS", "RSSR", "SSRR", "SRRR", "RRSR", "SSSR"))),
  set_prior("normal(0, 9/4)", class = "b", coef = "macrolides", dpar = paste0("mu", c("RRRR", "RRSS", "RSSR", "SSRR", "SRRR", "RRSR", "SSSR"))),
  set_prior("normal(0, 9/4)", class = "Intercept", dpar = paste0("mu", c("RRRR", "RRSS", "RSSR", "SSRR", "SRRR", "RRSR", "SSSR"))), 
  set_prior("inv_gamma(0.01, 0.01)", class = "sd", dpar = paste0("mu", c("RRRR", "RRSS", "RSSR", "SSRR", "SRRR", "RRSR", "SSSR")), 
            group = "facility")
)


formula <- outcome ~ 1 + (1 | facility) + year + day + month + fluoroquinolones + lincosamides + macrolides + beta_lactams
start_time <- Sys.time()
brms_fit <- brms.fit.model(data = sim_data, 
                           formula = formula, 
                           family = "categorical", 
                           priors = priors, 
                           n_iter = 5000, 
                           n_cores = 4, 
                           n_chains = 4, 
                           seed = 123, 
                           folder_path_prefix = paste0("results/", org_short, "/aim2-2/simReg/"))
end_time <- Sys.time()
(run_time <- end_time - start_time)
save(brms_fit, run_time, file = paste0("results/SA/aim2-2/simReg/mind_aim2-2_simReg_", unique(sim_data$sim_id), "_brms_fit.RData"))


# Next simulation
load("data_output/simReg/SA_20072021_2_3/mind_aim2-2_simReg_SA_20072021_2_3.RData")
sim_data <- simReg$data

formula <- outcome ~ 1 + (1 | facility) + year + day + month + fluoroquinolones + lincosamides + macrolides + beta_lactams
start_time <- Sys.time()
brms_fit <- brms.fit.model(data = sim_data, 
                           formula = formula, 
                           family = "categorical", 
                           priors = priors, 
                           n_iter = 5000, 
                           n_cores = 4, 
                           n_chains = 4, 
                           seed = 123, 
                           folder_path_prefix = paste0("results/", org_short, "/aim2-2/simReg/"))
end_time <- Sys.time()
(run_time <- end_time - start_time)
save(brms_fit, run_time, file = paste0("results/SA/aim2-2/simReg/mind_aim2-2_simReg_", unique(sim_data$sim_id), "_brms_fit.RData"))

# Next simulation
# load("data_output/simReg/SA_20072021_2_4/mind_aim2-2_simReg_SA_20072021_2_4.RData")
# sim_data <- simReg$data
# 
# formula <- outcome ~ 1 + (1 | facility) + year + day + month + fluoroquinolones + lincosamides + macrolides + beta_lactams
# start_time <- Sys.time()
# brms_fit <- brms.fit.model(data = sim_data, 
#                            formula = formula, 
#                            family = "categorical", 
#                            priors = priors, 
#                            n_iter = 5000, 
#                            n_cores = 4, 
#                            n_chains = 4, 
#                            seed = 123, 
#                            folder_path_prefix = paste0("results/", org_short, "/aim2-2/simReg/"))
# end_time <- Sys.time()
# (run_time <- end_time - start_time)
# save(brms_fit, run_time, file = paste0("results/SA/aim2-2/simReg/mind_aim2-2_simReg_", unique(sim_data$sim_id), "_brms_fit.RData"))
# 
