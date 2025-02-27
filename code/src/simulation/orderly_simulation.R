# File for initilializing Orderly repository
library(orderly2)
setwd("/Users/tm-pham/academia/hsph/mind/publications/aim2-2/code/")

orderly_run("simulation_loop_scenario2", 
            parameters = list(mean = 0.0025, sd = 0.1, 
                              n_sim = 1000, 
                              org_short = "SA", data_id = "SA_20072021_4", 
                              n_facilities = 70, n_obs_per_facility = 1000))

orderly_run("simulation_mblogit_fit")

orderly_run("simulation_plots")




