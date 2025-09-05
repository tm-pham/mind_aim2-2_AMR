# ============================================================================ #
# Project: MInD Aim 2.2
# FUNCTION
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# Title: Create importation terms for multinomial logistic regression
# ---------------------------------------------------------------------------- #
# For each calendar date (after a chosen start date) on which there are isolates 
# in df_1, it looks back prev_time days in df_2 within the same facilities and 
# computes, by facility × antibiogram, how many prior “importations” of that
# phenotype occurred, then turns those counts into proportions/rates. 
# It returns (and also saves) a long data frame with one row 
# per date × facility × antibiogram.
# ============================================================================ #

create.importation.terms <- function(df_1, # usually df_1
                                     df_2, # usually df_2 
                                     org_short, 
                                     start_date = "2007-01-31", 
                                     prev_time = 90, # time unit = days
                                     data_path, 
                                     data_id){
  # -------------------------------------------------------------------------- #
  # Data preparation
  # Set antibiograms levels for df_2 (equal to df_1)
  antibiogram_levels <- levels(df_1$antibiogram_name)
  if(is.null(antibiogram_levels)){
    antibiogram_levels <- sort(unique(df_1$antibiogram_name))
  }
  # Filter dayes after start date
  dates <- sort(unique(df_1$specimentakendatetime))
  dates <- as.Date(dates, format="%Y-%m-%d")
  dates <- dates[dates>as.Date(start_date, format="%Y-%m-%d")]
  
  # -------------------------------------------------------------------------- #
  # Loop through dates to calculate importation terms
  options(dplyr.summarise.inform = FALSE)
  df_micro_imp <- NULL
  for(i in 1:length(dates)){
    current_date <- dates[i]
    # Facilities with isolates on current date 
    current_fac <- df_1 %>% 
      filter(specimentakendatetime == current_date) %>% 
      select(sta6a) %>% unlist() %>% unique()
    temp <- df_2 %>% 
      mutate(antibiogram_name = factor(antibiogram_name, levels = antibiogram_levels)) %>% 
      # Filter data from the last prev_time days (usually 90)
      filter(specimentakendatetime < current_date, 
             specimentakendatetime >= current_date-prev_time, 
             sta6a %in% current_fac, 
             !is.na(antibiogram_name)) %>% 
      mutate(sta6a = factor(sta6a, levels = current_fac)) %>%
      ungroup() %>% 
      group_by(sta6a, antibiogram_name) %>% 
      # Count the number of importations of that specific antibiogram/phenotype per facility
      summarise(n_imp_prev=n(), 
                date = current_date) %>% 
      ungroup() %>% 
      complete(sta6a, antibiogram_name, fill = list(n_imp_prev=0, date = current_date)) %>%
      ungroup() %>% 
      group_by(sta6a) %>%
      # Proportion of antibiograms/phenotypes in previous prev_time days (usually 90)
      mutate(n_total_fac = sum(n_imp_prev, na.rm=T), 
             importation_perc = n_imp_prev/n_total_fac, 
             importation_rate = 1000*n_imp_prev/n_total_fac) %>%
      arrange(sta6a, date)
    df_micro_imp <- rbind(df_micro_imp, temp)
  }
  options(dplyr.summarise.inform = TRUE)
  
  # Save data 
  save(df_micro_imp, file=paste0("data/", org_short, "/aim2-2/", data_id, "/mind_aim2-2_", data_id, "_df_micro_imp.RData"))
  
  return(df_micro_imp)
}