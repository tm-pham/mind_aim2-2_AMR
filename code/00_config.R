# ============================================================================ #
# Project: MInD Aim 2.2
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# Title: Configuration file for figures (Epidemics conference)
# ============================================================================ #

################################################################################
# Required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggh4x) # facet_wrap_nested

################################################################################
# Paths
FIGURES <- "/Users/tm-pham/academia/hsph/mind/publications/aim2-2/figures"
RESULTS <-  "/Users/tm-pham/academia/hsph/mind/publications/aim2-2/results"


source("plotting_template.R")


################################################################################
# Order of pathogens
bugs_ordered <- c(
  "Staphylococcus aureus", 
  "Escherichia coli", 
  "Klebsiella pneumoniae",
  "Pseudomonas aeruginosa"
)

################################################################################
# Phenotype order
FQL_phenotype_order <- c(
  "R-R-R",
  "R-R-S",
  "R-S-R",
  "R-S-S",
  "S-R-R",
  "S-R-S",
  "S-S-R",
  "S-S-S"
)

ASBL_phenotype_order <- c(
  "R-R-R",
  "R-R-S",
  "S-R-R",
  "R-S-R",
  "S-R-S",
  "R-S-S",
  "S-S-R",
  "S-S-S"
)

GC_phenotype_order <- c(
  "R-R-R",
  "R-R-S",
  "S-R-R",
  "S-R-S",
  "R-S-R",
  "R-S-S",
  "S-S-R",
  "S-S-S"
)

CPM_phenotype_order <- c(
  "R-R-R",
  "R-R-S",
  "S-R-R",
  "S-R-S",
  "R-S-R",
  "S-S-R",
  "R-S-S",
  "S-S-S"
)

################################################################################
# Combination Organism-phenotype order
org_ab_order <- c(
  # Staphylococcus aureus
  paste0("Staphylococcus aureus__", FQL_phenotype_order),
  # Escherichia coli
  paste0("Escherichia coli__", FQL_phenotype_order),
  # Klebsiella pneumoniae
  paste0("Klebsiella pneumoniae__", FQL_phenotype_order),
  # Pseudomonas aeruginosa
  paste0("Pseudomonas aeruginosa__", FQL_phenotype_order)
)


facet_labs <- c("<b><i>Staphylococcus aureus</i></b><br>FQL - ASBL - MACR",
                "<b><i>Escherichia coli</i></b><br>FQL - 3/4GC - BL/BLI",
                "<b><i>Klebsiella pneumoniae</i></b><br>FQL - 3/4GC - BL/BLI",
                "<b><i>Pseudomonas aeruginosa</i></b><br>FQL - BL/BLI - CPM")
facet_table <- as.data.frame(cbind(organismofinterest = bugs_ordered, 
                                   facet_labs = facet_labs))


################################################################################
# Define colors for each species-phenotype combination -------------------------
pathogen_colors <- c(
  "Staphylococcus aureus" = "#E69F00",
  "Escherichia coli" = "#00496FFF",
  "Klebsiella pneumoniae" = "#D46780FF",
  "Pseudomonas aeruginosa" = "#0E7C7BFF"
)

phenotype_colors <- c(
  ## Staphylococcus aureus
  "Staphylococcus aureus__R-R-R" = "#664700",
  "Staphylococcus aureus__R-R-S" = "#B27B00",
  "Staphylococcus aureus__R-S-R" = "#B27B00",
  "Staphylococcus aureus__S-R-R" = "#B27B00",
  "Staphylococcus aureus__R-S-S" = "#FFB000",
  "Staphylococcus aureus__S-R-S" = "#FFB000",
  "Staphylococcus aureus__S-S-R" = "#FFB000",
  "Staphylococcus aureus__S-S-S" = "#FFC84D",
  
  ## Escherichia coli
  "Escherichia coli__R-R-R" = "#004366",
  "Escherichia coli__R-R-S" = "#0075B2",
  "Escherichia coli__R-S-R" = "#0075B2",
  "Escherichia coli__S-R-R" = "#0075B2",
  "Escherichia coli__R-S-S" = "#00A8FF",
  "Escherichia coli__S-R-S" = "#00A8FF",
  "Escherichia coli__S-S-R" = "#00A8FF",
  "Escherichia coli__S-S-S" = "#4DC2FF",
  
  ## Klebsiella pneumoniae
  "Klebsiella pneumoniae__R-R-R" = "#6A0F38",
  "Klebsiella pneumoniae__R-R-S" = "#98154D",
  "Klebsiella pneumoniae__R-S-R" = "#98154D",
  "Klebsiella pneumoniae__S-R-R" = "#98154D",
  "Klebsiella pneumoniae__R-S-S" = "#C31E6E",
  "Klebsiella pneumoniae__S-R-S" = "#C31E6E",
  "Klebsiella pneumoniae__S-S-R" = "#C31E6E",
  "Klebsiella pneumoniae__S-S-S" = "#E661A5",
  
  ## Pseudomonas aeruginosa
  "Pseudomonas aeruginosa__R-R-R" = "#0A4F4E",
  "Pseudomonas aeruginosa__R-R-S" = "#0E6A69",
  "Pseudomonas aeruginosa__R-S-R" = "#0E6A69",
  "Pseudomonas aeruginosa__S-R-R" = "#0E6A69",
  "Pseudomonas aeruginosa__R-S-S" = "#158886",
  "Pseudomonas aeruginosa__S-R-S" = "#158886",
  "Pseudomonas aeruginosa__S-S-R" = "#158886",
  "Pseudomonas aeruginosa__S-S-S" = "#4AA3A2"
)

# Trend colors
trend_colors <- c(
  "Increasing" = "antiquewhite4",    
  "Decreasing" = "#FFFFFF",    
  "Stable" = "#FFFFFF"         
)
