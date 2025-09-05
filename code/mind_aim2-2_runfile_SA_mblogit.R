# ============================================================================ #
# Project: MInD Aim 2.2
# Run multinomial logistic regression
# Organism: S. aureus
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ---------------------------------------------------------------------------- #
# Drugs: FQL, Antistaphylococcal Beta-Lactams, Macrolides
# ============================================================================ #
# Working directory
setwd("path/to/folder")

# Load functions
source("code/plotting_template.R")
source("code/mind_aim2-2_function_mblogit_coef_manuscript_plot.R")

# Load packages
library(mclogit)
library(dplyr)

################################################################################
# Load data 
# Staphylococcus aureus
data_id <- "SA_20072021_4"
org <- "Staphylococcus aureus"
org_short <- "SA"

df_input <- read.csv(paste0("path/to/folder/data/mind_aim2-2_", data_id, "_df_input_multinom_regr.csv", header = T))

################################################################################
# Set reference categories 
df_input <- df_input %>% 
  mutate(year_numeric = as.numeric(date_year),
         facility_rurality = factor(facility_rurality, levels = c("U","R")), 
         complexitylevel = factor(complexitylevel, levels = c("1a", "1b", "1c", "2", "3", "Excl")), 
         census_region = factor(census_region, levels = c("south", "west", "midwest", "northeast", "pr and outlying areas")), 
         sampletype_spec = factor(sampletype_spec, levels = c("SPUTUM", "OTHER", "BLOOD", "WOUND", "FOOT", "SWAB", "SKIN")), 
         gender = factor(gender, levels = c("M", "F")), 
         race = factor(race, levels = c("white", 
                                        "black or african american", 
                                        "missing", 
                                        "mixed", 
                                        "native hawaiian or other pacific islander", 
                                        "american indian or alaska native", 
                                        "asian"))) %>% 
  filter(!complexitylevel%in%c("Excl", "3"), antibiogram_name!="Other") %>%
  mutate(antibiogram_name = factor(antibiogram_name, levels = c("AB_2", paste0("AB_", c(1,3,4,5,6)))))

################################################################################
# Fit model
fit <- mclogit::mblogit(formula = antibiogram_name ~ 1 + dayofweek + monthofyear + year_numeric  + 
                          rate_prev_num_antistaphbetalactams + 
                          rate_prev_num_class_macr + 
                          rate_prev_num_class_flq +
                          AB_1_importation_rate + 
                          AB_2_importation_rate +
                          AB_3_importation_rate +
                          AB_4_importation_rate +
                          AB_5_importation_rate +
                          AB_6_importation_rate +
                          census_region +
                          facility_rurality + 
                          complexitylevel + 
                          median_cci + 
                          median_los + 
                          perc_black + 
                          perc_female
                        , 
                        random = ~1|sta6a,
                        data = df_input)

suffix = "main_analysis"
save(fit, file = paste0("results/", org_short, "/aim2-2/", data_id, "/mind_aim2-2_", data_id, "_mblogit_", suffix, ".RData"))

################################################################################
# Plot
mblogit.coef.plot(model = fit, 
                  org = org_short, 
                  abx_order = c("flq", "ASBL",  "macr"),
                  abx_plot_names = c("FQL", 
                                     "ASBL", 
                                     "MACR"),
                  plot_width = 6, 
                  plot_height = 14,
                  suffix = paste0("figure3_per100d_"), 
                  exponentiate = T,
                  outcome_option = 2, 
                  round_digits = 3, 
                  figure_path =  paste0("path/to/folder/figures/", org_short, "/aim2-2/", data_id, "/"),
                  data_path =  paste0("path/to/folder/results/", org_short, "/aim2-2/", data_id, "/"))
