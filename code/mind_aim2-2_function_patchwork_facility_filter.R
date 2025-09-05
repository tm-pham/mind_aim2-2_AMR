# ============================================================================ #
# FUNCTION
# Title: Patchwork facility filtering approach
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ---------------------------------------------------------------------------- #
# This function filters facilities based on minimal missingness in their
# antibiotic susceptibility test results. 
# In this version, the resulting dataset will be a composite of different 
# facilities in different time segments.
# This function is applicable to a vector of organisms and a vector of 
# antibiotics (usually from the antibiogram panel). 
# ============================================================================ #
# Working directory
setwd("P:/ORD_Samore_202109019D/Mui/")

# Load packages
library(dplyr)

# Load functions
source("code/mind_aim2-2_function_filter_facilities_missing_below_threshold.R")


# Function to merge datasets
patchwork.filter <- function(df_micro,                  # Microbiology dataset
                             org_vec,                   # Vector with organism names
                             org_short_vec,             # Vector with abbreviations for organisms (for saving files)
                             abx_list,                  # Vector of antibiotics (usually panel of antibiogram)
                             threshold = 0.1,           # Maximum proportion of allowed missing values, between 0 and 1
                             time_unit = "date_qrt",    # Time unit for calculating missingness
                             n_t_consec = 4,            # Number of time segments over which missingness has to be less than threshold
                             time_period = c(2007,2021),# Overall time period as vector: first entry is start year, last entry is end year 
                             data_path = "P:/ORD_Samore_202109019D/Mui/data/", 
                             suffix = "", 
                             return_micro = T,          # Flag whether dataset should be returned
                             descriptive_stats = T){    # FLag whether tables with descriptives should be created     
  # Determine the time segements in years
  if(time_unit=="date_qrt") n_consec_years <- n_t_consec/4
  
  # List to save resulting microbiology subset
  patchwork_dataset <- list()
  
  # Loop through n_consec_years periods and apply filter
  for(start_year in seq(time_period[1], time_period[2], by = n_consec_years)){
    end_year <- start_year + n_consec_years - 1
    # This returns a list with eligible facilities for each antibiotic class
    # If there n antibiotic classes, then the list will have length n
    fac_list <- fac.missing.years.multiple.ab(df_micro %>% 
                                                filter(date_year >= start_year, 
                                                       date_year <= end_year), 
                                              org_vec = org_vec, 
                                              abx_list = abx_list, 
                                              threshold = threshold, 
                                              time_unit = time_unit, 
                                              n_t_consec = n_t_consec)
    
    fac_list_all <- list()
    time_names_list <- names(fac_list[[1]][[1]])
    for(o in 1:length(org_vec)){
      org <- org_vec[o]
      fac_list_all[[org]][[time_names_list]] <- Reduce(intersect, lapply(fac_list[[org]], function(x) x[[time_names_list]]))
    }
    
    for(i in 1:length(org_vec)){
      if(length(unique(start_year, end_year))> 1){
        time_segment_name <- paste(start_year, end_year, sep="-")
      }else{
        time_segment_name <- start_year
      }
      
      patchwork_dataset[[org]][[as.character(time_segment_name)]] <- df_micro %>% 
        filter(date_year >= start_year, date_year <= end_year, 
               sta6a%in%fac_list_all[[org]][[time_names_list]])
    }
    
  }
  
  patchwork_combined <- lapply(org_vec, function(x) do.call("rbind", patchwork_dataset[[x]]))
  names(patchwork_combined) <- org_vec
  # Save into RData
  for(i in 1:length(org_vec)){
    org_short <- org_short_vec[i]
    rownames(patchwork_combined[[i]]) <- NULL
    
    data_id <- paste0(org_short, "_", suffix)
    
    if(descriptive_stats){
      df_n_year <- patchwork_combined[[i]] %>% 
        group_by(date_year) %>% 
        summarize(n = n())
      
      df_n_year_total <- df_micro %>% 
        group_by(date_year) %>% 
        summarize(n_total = n())
      
      df_n_year <- left_join(df_n_year, df_n_year_total, by = "date_year") %>% 
        mutate(Percentage = n/n_total)
      
      dir.create(file.path(paste0(data_path, "/", org_short, "/aim2-2/", data_id, "/")), showWarnings = FALSE)
      
      write.table(df_n_year, file =  paste0(data_path, "/", org_short, "/aim2-2/", org_short, "_", suffix, "/mind_aim2-2_", org_short, "_", suffix, "_", n_consec_years, "y", "_patchwork_n_isolates.csv"),
                row.names = F, sep = ",")
      
      df_n_fac <- patchwork_combined[[i]] %>% 
        group_by(date_year) %>% 
        summarize(n_sta6a = length(unique(sta6a)))
      
      write.table(df_n_fac, file =  paste0(data_path, "/", org_short, "/aim2-2/", org_short, "_", suffix, "/mind_aim2-2_", org_short, "_", suffix, "_", n_consec_years, "y", "_patchwork_n_facilities.csv"),
                row.names = F, sep = ",")
      
    }
    
    # Meta data of patchwork data 
    meta_data <- list(threshold = threshold, time_unit = time_unit, n_t_consec = n_t_consec, abx_list[[org_vec[i]]])
    
    dir.create(file.path(paste0(data_path, "/", org_short, "/aim2-2/", data_id, "/")), showWarnings = FALSE)
    save(patchwork_combined, meta_data, file = paste0(data_path, "/", org_short, "/aim2-2/", org_short, "_", suffix, "/mind_aim2-2_", org_short, "_", suffix, "_", n_consec_years, "y", "_patchwork_data.RData"))
  }
  
  
  if(return_micro) return(patchwork_combined)
}



