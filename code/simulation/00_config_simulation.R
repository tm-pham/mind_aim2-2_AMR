# ============================================================================ #
# Project: MInD Aim 2.2
# Title: Configuration file for simulation study
# Author: Thi Mui Pham, mui.k.pham@mgail.com
# ============================================================================ #
# Load packages
library(dplyr)
library(here)
library(ggplot2)

# PATHS ------------------------------------------------------------------------
PROJECT_FOLDER <- "/Users/tm-pham/academia/hsph/mind/publications/aim2-2/"
setwd(PROJECT_FOLDER)

DATA <- "/Users/tm-pham/academia/hsph/mind/publications/aim2-2/data"
RESULTS <- "/Users/tm-pham/academia/hsph/mind/publications/aim2-2/results/simulation"
FIGURES <- "/Users/tm-pham/academia/hsph/mind/publications/aim2-2/figures/simulation"

# PARAMETERS FOR SIMULATION STUDY ----------------------------------------------
# Parameters for facility-level antibiotic use distributions
abx_params <- readRDS(here::here("data","simulation_abx_params_lognormal_FQL_variance.RDS"))

# Log-normal distribution parameters for regression coefficient gamma (facility-level previous abx use)
mean = 0.0025
sd = 0.1

# Number of simulations
n_sim = 5

# Generate gamma distribution for facility-level previous abx use coefficient
gamma_distr <- rlnorm(n_sim, mean = log(mean), sd = sd)

# File name prefix for simulation data
PREFIX = "mind_aim2-2_simData_"

# Antibiotic labels string
abx_labels = c("Drug 1", "Drug 2", "Drug 3")



