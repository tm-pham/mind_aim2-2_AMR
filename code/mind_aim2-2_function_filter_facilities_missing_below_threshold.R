# ============================================================================ #
# Project: MInD Aim 2.2
# FUNCTION
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# Title: Filter facilities with a maximum proportion of missing values in x 
# consecutive years
# ---------------------------------------------------------------------------- #
# This file contains two functions: 
# fac.missing.years
# fac.missing.years.multiple.ab
# 
# This file is used in multiple other files, including:
# mind_aim2-2_function_patchwork_facility_filter.R
# mind_aim2-2_runfile_patchwork_facility_filter.R (indirectly)
# ============================================================================ #
setwd("path/to/folder/")

# Load packages
library(dplyr)

# ---------------------------------------------------------------------------- #
# Function: Facilities with missingness below a threshold for one antibiotic
# ----------------------------------------------------------------------------
# Input: 
# df_micro = Microbiology data set (similar to microbiology data from McKenna)
# org = Organism
# ab = antibiotic or antibiotic class
# threshold = maximal missingness proportion (between 0 and 1)
# time_unit = unit of time (e.g., month, quarter, year)
# n_t_consec = Number of consecutive time units
# 
# Output:
# fac_time_subset_list = List of facilities with missingness not exceeding 
# threshold per time period (indicated by years)
# 
# Example:
# ab = FQL_class
# time_unit = date_qrt
# n_t_consec = 16 
# threshold = 0.2
# 
# Output will be sta6a ids of facilities where the proportion of missing values 
# in susceptibility test results for FQL class is at most 20% in each quarter 
# for 16 consecutive quarters (i.e., 4 years)
# ---------------------------------------------------------------------------- #
facility.missingness.filter <- function(df_micro, 
                                        ab, 
                                        threshold = 0.01, 
                                        time_unit = "date_qrt",
                                        n_t_consec = 16){
  # Total time period in the required time units
  time_period <- sort(unique(df_micro[, time_unit]))
  
  # This variable contains the names of the time period subsets
  # Columns are the respective time period subsets
  time_subsets_ind <- sapply(1:(length(time_period)-n_t_consec), function(x) seq(x, x+n_t_consec-1))
  
  # Only choose subsets that start at the beginning of the year
  # Example: time_period[time_subsets_ind] = subset of the total time period
  if(time_unit=="date_qrt") time_subsets_ind <- time_subsets_ind[, seq(1, ncol(time_subsets_ind), by=4)]
  if(time_unit=="date_month") time_subsets_ind <- time_subsets_ind[, seq(1, ncol(time_subsets_ind), by=12)]
  
  # Data frame with proportion of missingness per facility
  # Remove isolates with missings in all relevant abx
  perc_miss_time <- df_micro %>% 
    select(sta6a, all_of(ab), eval(time_unit)) %>%
    group_by(eval(parse(text = time_unit)), sta6a) %>% 
    summarize(miss_perc = sum(!(eval(parse(text=ab))%in%c("S","R","I")))/n()) %>% 
    ungroup()
  colnames(perc_miss_time)[1] <- "t_unit"
  
  # Only relevant for time_unit=date_qrt
  remove_Q <- paste0(" Q", seq(1,4))
  
  fac_time_subset_list <- list()
  names_list <- vector(mode="character")
  
  if(is.null(ncol(time_subsets_ind))) time_subsets_ind <- as.data.frame(time_subsets_ind, ncol=1)
    
  for(i in 1:ncol(time_subsets_ind)){
    time_subset <- time_period[time_subsets_ind[,i]]
    # Name the subset after the years
    names_list[i] <- paste(unique(stringr::str_extract(time_subset, "\\d+")), collapse=" ")
    
    # Filter facilities where missingness is below threshold
    df_elig_fac <- perc_miss_time %>% 
      filter(t_unit%in%time_subset) %>% 
      ungroup() %>%
      group_by(sta6a) %>% 
      filter(miss_perc <= threshold) %>% 
      mutate(elig = n() >= n_t_consec) %>% 
      filter(elig == T) 
    
    fac_time_subset_list[[i]] <- unique(df_elig_fac$sta6a)
  }
  names(fac_time_subset_list) <- names_list
  
  return(fac_time_subset_list)
}

# ---------------------------------------------------------------------------- #
# Function: Facilities with missingness below a threshold for multiple abx and 
# organisms
# ----------------------------------------------------------------------------
# Output:
# List of eligible facilities per organism and per antibiotic (class) 
# ---------------------------------------------------------------------------- #
fac.missing.years.multiple.ab <- function(df_micro, 
                                          org_vec, 
                                          abx_list, 
                                          threshold = 0.01, # between 0 and 1
                                          time_unit = "date_qrt",
                                          n_t_consec = 16, 
                                          reduce_if_one = F){
  fac_org_list <- fac_list <-list()
  for(o in 1:length(org_vec)){
    org <- org_vec[o]
    cat(org, "\n")
    micro_temp <- df_micro %>% 
      filter(organismofinterest==org) %>%
      mutate(across(abx_list[[org]], ~ifelse(.=="N", NA, .))) %>%
      filter(!if_all(abx_list[[org]], ~ is.na(.)))
    print(paste0("Number of antibiotics: ", length(abx_list[[org]])))
    for(i in 1:length(abx_list[[org]])){
      ab <- abx_list[[org]][i]
      fac_org_list[[org]][[ab]] <- facility.missingness.filter(df_micro = micro_temp, 
                                                               ab = abx_list[[org]][i], 
                                                               threshold = threshold, 
                                                               time_unit = time_unit, 
                                                               n_t_consec = n_t_consec)
      
    }
  }
  
  if(reduce_if_one){
    if(length(org_vec)==1){
      return(fac_org_list[[1]])
    }else{
      return(fac_org_list)
    }
  }else{
    return(fac_org_list)
  }

}
