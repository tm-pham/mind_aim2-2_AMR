# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Simulation analysis
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# ---------------------------------------------------------------------------- #
# Description: Run files for simulation study in order
# ============================================================================ #
# Set up environment -----------------------------------------------------------
PROJECT_FOLDER <- "/Users/tm-pham/academia/hsph/mind/publications/aim2-2/"
setwd(PROJECT_FOLDER)
here::i_am("code/simulation/runfile_simulation_analysis.R")

# Load configuration settings --------------------------------------------------
cat("Loading configuration settings...\n")
source(here::here("code","simulation","00_config_simulation.R"))
cat("✓ Configuration settings loaded\n")

# Load all custom functions ----------------------------------------------------
cat("Loading custom functions...\n")
function_files <- list.files(
  here("code", "simulation"),
  pattern = "01_function*",
  full.names = TRUE
)
invisible(lapply(function_files, source))
cat("✓", length(function_files), "function files loaded\n")

# Generate simulation data -----------------------------------------------------
source(here::here("code", "simulation", "02_runfile_simulation_data.R"))
cat("✓ Simulation data generated\n")
# Fit models -------------------------------------------------------------------
source(here::here("code", "simulation", "03_runfile_mblogit_fit.R"))
cat("✓ Models fitted to simulation data\n")
# Plot results -----------------------------------------------------------------
source(here::here("code", "simulation", "04_create_output.R")) 
cat("✓ Summary statistics created and results plotted.\n")

