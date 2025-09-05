# ============================================================================ #
# Project: MInD Aim 2.2
# FUNCTION
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# Title: Data preparation for multinomial logistic regression 
# ---------------------------------------------------------------------------- #
# Function for data preparation for multinomial logistic regression
# 1. Adds antibiotic use data to microbiology dataset
# Either restricts the dataset to the most common antibiograms (and removes all
# others) or includes an additional level called "Other" 
# 2. Adds importation term
# Note: This function is typically called after patchwork dataset is created. 
# ============================================================================ #
# Working directory
setwd("path/to/folder/")

# Load packages
source("code/packages.R")
source("code/mind_aim2-2_function_importation_terms.R")


data.prep.multinom.log.regr <- function(df, 
                                        df_micro = micro_data_inc30d, 
                                        org, 
                                        org_short, 
                                        bug_drugs, 
                                        data_id, 
                                        df_abx_prev, 
                                        lookback = "14d", 
                                        n_antibiograms = 6, 
                                        include_other = F, 
                                        start_date = "2007-03-31",
                                        onset_filter = "Hospital-onset", 
                                        create_importation = T, 
                                        imp_prev_days = 90, 
                                        suffix = "patchwork_2007_2021", 
                                        save_result = T){
  # -------------------------------------------------------------------------- #
  # 0. Add antibiograms (phenotypes) to microbiology dataset
  # Create antibiogram/phenotype table
  ab_table <- create.antibiograms(df_micro = df %>% 
                                    filter(organismofinterest%in%org, onset%in%onset_filter), 
                                  bugs = org, 
                                  bugs_short = org_short, 
                                  drugs = bug_drugs, 
                                  data_path = "data")
  print(head(ab_table, 10))
  
  df_micro_ab <- left_join(df, ab_table %>% 
                             dplyr::select(antibiogram_name, ab_name, bug_drugs[[org]]), 
                           by = bug_drugs[[org]])
  
  # -------------------------------------------------------------------------- #
  # 1. Add antibiotic use
  # Columns of antibiotic day level dataset that are relevant
  abx_colnames <- c("Sta6a", "datevalue", "dayofweek", "IsWeekend", "weekofyear", 
                    "monthofyear", "quarter", "FederalHolidayFlag", "num_patient_days", 
                    colnames(df_abx_prev)[str_detect(colnames(df_abx_prev), "prev_")])
  abx_prev_colnames <- colnames(df_abx_prev)[str_detect(colnames(df_abx_prev), "prev_")]
  
  # Either classify all antibiograms that are not among the top n_antibiograms 
  # as "Other" or remove them from the dataset
  antibiogram_levels <- paste0("AB_", 1:n_antibiograms)
  if(include_other){
    df_micro_ab <- df_micro_ab %>% 
      mutate(antibiogram_name = ifelse(antibiogram_name%in% antibiogram_levels, antibiogram_name, "Other"))
    antibiogram_levels <- c(antibiogram_levels, "Other")
    df_micro_ab <- df_micro_ab %>% 
      mutate(antibiogram_name = factor(antibiogram_name, levels=antibiogram_levels))
  }else{
    df_micro_ab <- df_micro_ab %>% 
      dplyr::filter(antibiogram_name%in%antibiogram_levels)
  }
  
  # Add the proportion of patient days of antibiotics use to the microbiology data
  min_date <- min(df_micro_ab$specimentakendatetime)
  df_micro_abx <- df_micro_ab %>% 
    mutate(ab_name = ifelse(antibiogram_name=="Other", "Other", ab_name)) %>% 
    dplyr::filter(onset=="Hospital-onset") %>% 
    left_join(df_abx_prev %>% select(all_of(abx_colnames)), 
              by = c("sta6a"="Sta6a", "specimentakendatetime"="datevalue")) %>% 
    mutate(across(abx_prev_colnames, 
                  ~./prev_num_patient_days, 
                  .names = "prop_{.col}"), 
           across(abx_prev_colnames, 
                  ~.*1000/prev_num_patient_days, 
                  .names = "rate_{.col}"),
           t_i = as.numeric(specimentakendatetime-min_date))
  # Save input data for multinomial logistic regression
  save(df_micro_abx, ab_table, 
       file=paste0("data/", org_short, "/aim2-2/", data_id, "/mind_aim2-2_", data_id, "_df_micro_abx_HOI_patchwork.RData"))
  
  # Dates in the input data
  dates <- sort(unique(df_micro_abx$specimentakendatetime))
  
  # -------------------------------------------------------------------------- #
  # 2. Add importation term
  if(create_importation){
    fac <- df_micro_abx %>% select(sta6a) %>% unique() %>% unlist() %>% unname()
    time_period = unique(df$date_year)
    antibiogram_levels <- sort(unique(df_micro_ab$antibiogram_name))
    ab_levels <- unique(df_micro_abx$ab_name)
    df_micro_COI <- df_micro %>% filter(organismofinterest == org, 
                                           onset=="Community-onset", 
                                           sta6a%in%fac, 
                                           date_year%in%time_period) %>% 
      left_join(ab_table, 
                by = bug_drugs[[org]]) %>% 
      # Some antibiograms may not be represented in df_ab_table
      mutate(antibiogram_name = if_else(is.na(antibiogram_name), "Other", as.character(antibiogram_name)), 
             ab_name = if_else(is.na(ab_name), "Other", as.character(ab_name)))
    
    df_micro_imp <- create.importation.terms(df_1 = df_micro_abx, 
                                             df_2 = df_micro_COI, 
                                             org_short = org_short,
                                             start_date = start_date, 
                                             prev_time = imp_prev_days, 
                                             data_path = "P:/ORD_Samore_202109019D/Mui/data/", 
                                             data_id = data_id)
    
    
    df_micro_imp_wide <- pivot_wider(df_micro_imp %>% 
                                       select(date, sta6a, antibiogram_name, importation_rate, importation_perc), 
                                     names_from = antibiogram_name, 
                                     names_glue = "{antibiogram_name}_{.value}", 
                                     values_from = c(importation_rate, importation_perc))
    
    # Merge with micro data
    df_input <- left_join(df_micro_abx, 
                          df_micro_imp_wide,
                          by=c("specimentakendatetime" = "date",
                               "sta6a" = "sta6a")) %>% 
      filter(specimentakendatetime > as.Date(start_date, format="%Y-%m-%d")) %>% 
      left_join(df_micro_imp %>% select(antibiogram_name, sta6a, date, importation_perc, importation_rate), 
                by = c("specimentakendatetime" = "date",
                       "sta6a" = "sta6a", 
                       "antibiogram_name" = "antibiogram_name")) %>% 
      # Convert time variables to factor variables
      mutate(dayofweek = factor(dayofweek), 
             monthofyear = factor(monthofyear), 
             date_year = factor(date_year))
    # Save the result
    save(df_input, ab_table, file=paste0("data/", org_short, "/aim2-2/", data_id, "/mind_aim2-2_", data_id, "_df_input_", suffix, ".RData"))
    return(df_input)
  }else{
    
    return(df_micro_abx)
  }
}
