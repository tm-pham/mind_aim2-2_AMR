# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Fit multinomial logistic regression to simulated data and summarize results
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ---------------------------------------------------------------------------- #
# Frequentist model using mblogit package
# This script runs a multinomial logistic regression model based on 
# simulation data created in the simulation_loop scripts and summarizea the 
# results in terms of bias, coverage, and confidence interval width/envelope
# ============================================================================ #
# Load simulation ids
load(here::here("results", "simulation", data_id, "sim_ids.RData"))

# ------------------------------------------------------------------------------
# Run frequentist model
model_fit <- list()
df_list <- list()
coef_table <- list()

for(i in 1:length(sim_ids)){
  cat("Loading ", i, "th simulation data for", sim_ids[i], "\n")
  load(here::here("results", "simulation", data_id, sim_ids[i], paste0(PREFIX, sim_ids[i], ".RData"))) # simReg, args_list
  abx_order <- names(simReg$abx_params)
  cat("Antibiotic classes:", abx_order, "\n")
  # Create the base formula string
  base_formula <- "outcome ~ 1 + year + day + month"
  # Dynamically append the antibiotic terms from abx_order
  for (antibiotic in abx_order) {
    base_formula <- paste(base_formula, "+", antibiotic)
  }
  
  # Convert the string to a formula object
  formula <- as.formula(base_formula)
  
  # Fit the formula using the formula
  cat("Fitting multinomial logistic regression model using the formula", base_formula, "\n")
  tryCatch({
    invisible(capture.output({
      model <- model_fit[[i]] <- mblogit.fit.model(data = simReg$data, formula = formula, random_effects =  ~1|facility)
    }))
  }, error = function(e) {
    cat("Error in fitting model", e$message, "\n")
  }) 
  
  tt <- broom::tidy(model, conf.int=TRUE)
  df_vars <- as.data.frame(do.call("rbind", stringr::str_split(tt$term, "~")))
  colnames(df_vars) <- c("antibiogram", "term")
  
  df_coef <- cbind(df_vars[,], tt[, -1]) %>% 
    mutate(term = stringr::str_remove_all(term, "^.*_prev_num_"), 
           term = stringr::str_remove_all(term, "class_")) %>% 
    filter(term%in%abx_order)
  
  positions <- paste0("pos", 1:length(abx_order))
  
  # Separate the antibiogram column into the dynamic positions
  df_list[[i]] <- df_coef %>%
    tidyr::separate(antibiogram, into = positions, sep = "-", remove = FALSE) %>%
    mutate(across(all_of(positions), ~ ifelse(. == "R", abx_order[as.numeric(sub("pos", "", cur_column()))], NA))) %>%
    rowwise() %>%
    mutate(
      flag = sum(across(all_of(positions), ~ ifelse(. == abx_order[as.numeric(sub("pos", "", cur_column()))] & term == abx_order[as.numeric(sub("pos", "", cur_column()))], 1, 0)), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(flag = ifelse(flag > 0, 1, 0))
  
  # Load the true coefficients
  coef_table[[i]] <- as.data.frame(cbind(true_coef = simReg$params$gamma, term = abx_order))
  coef_table[[i]]$true_coef <- as.numeric(coef_table[[i]]$true_coef)
}

# Combine the dataframes across all sim_ids 
df_combined <- bind_rows(df_list, .id = "sim_id") 
coef_combined <- bind_rows(coef_table, .id = "sim_id")