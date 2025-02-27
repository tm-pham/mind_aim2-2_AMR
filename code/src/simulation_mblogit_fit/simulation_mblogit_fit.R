# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Run multinomial logistic regression on simulated data
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ============================================================================ #
orderly2::orderly_dependency(
  "simulation_loop_scenario2", 
  "latest", 
  files = c("sim_ids.RData"))

# Print simulation ids
load("sim_ids.RData")

print(sim_ids)

################################################################################
# Run frequentist model
abx_order <- c("FQL", "ASBL", "MACR")
abx_labels <- c("Fluoroquinolones", "Antistaphylococcal Beta-lactams", "Macrolides")
model_fit <- list()
df_list <- list()
formula <- outcome ~ 1+ year + day + month + FQL + ASBL + MACR
for(i in 1:length(sim_ids)){
  load(paste0("/Users/tm-pham/academia/hsph/mind/publications/aim2-2/code/archive/simulation/", sim_ids[i], "/mind_aim2-2_simReg_", data_id, ".RData"))
  model <- model_fit[[i]] <- mblogit.fit.model(data = simReg$data, formula = formula, random_effects =  ~1|facility)
  tt <- broom::tidy(model, conf.int=TRUE)
  df_vars <- as.data.frame(do.call("rbind", stringr::str_split(tt$term, "~")))
  colnames(df_vars) <- c("antibiogram", "term")
  
  df_coef <- cbind(df_vars[ind,], tt[ind, -1]) %>% 
    mutate(term = stringr::str_remove_all(term, "^.*_prev_num_"), 
           term = stringr::str_remove_all(term, "class_")) %>% 
    filter(term%in%c("FQL", "ASBL", "MACR"))
  
  df_list[[i]] <- df_coef %>%
    tidyr::separate(antibiogram, into = c("pos1", "pos2", "pos3"), sep = "-", remove = FALSE) %>%
    mutate(
      pos1 = ifelse(pos1 == "R", "FQL", NA),
      pos2 = ifelse(pos2 == "R", "ASBL", NA),
      pos3 = ifelse(pos3 == "R", "MACR", NA), 
      flag = ifelse(pos1=="FQL" & term == "FQL", 1, 
                    ifelse(pos2=="ASBL" & term == "ASBL", 1, 
                           ifelse(pos3=="MACR" & term == "MACR", 1, 0))),
      flag = ifelse(is.na(flag), 0, flag)
    )
}
save(model_fit, file = "mblogit_fit.RData")

df_combined <- bind_rows(df_list, .id = "sim_id") 

coef_table <- as.data.frame(cbind(true_coef = params$gamma, term = abx_order))
coef_table$true_coef <- as.numeric(coef_table$true_coef)

bias_df <- df_combined %>%
  filter(flag == 1) %>% 
  left_join(coef_table, by = "term") %>%
  mutate(term = factor(term, levels = abx_order, labels = abx_labels)) %>% 
  group_by(sim_id, antibiogram, term) %>% 
  summarize(bias = estimate - true_coef)
print(bias_df)

coverage_df <- df_combined %>%
  filter(flag == 1) %>% 
  left_join(coef_table, by = "term") %>%
  mutate(term = factor(term, levels = abx_order, labels = abx_labels)) %>% 
  group_by(antibiogram, term) %>%
  summarize(
    coverage = sum(conf.low <= true_coef & conf.high >= true_coef, na.rm = TRUE)/length(sim_ids),
    .groups  = "drop"
  )
print(coverage_df)

# Calculate CI width and envelope metrics by grouping on antibiogram and term
ci_env_df <- df_combined %>%
  filter(flag == 1) %>% 
  mutate(term = factor(term, levels = abx_order, labels = abx_labels)) %>% 
  group_by(antibiogram, term) %>%
  summarise(
    # Average width of the confidence intervals
    avg_ci_width = mean(conf.high - conf.low, na.rm = TRUE),
    
    # Overall envelope: min of all lower bounds and max of all upper bounds
    overall_lower = min(conf.low, na.rm = TRUE),
    overall_upper = max(conf.high, na.rm = TRUE),
    envelope_width = overall_upper - overall_lower,
    
    # Alternative envelope metrics using quantiles (e.g., 5th and 95th percentiles)
    lower_5th = quantile(conf.low, probs = 0.05, na.rm = TRUE),
    upper_95th = quantile(conf.high, probs = 0.95, na.rm = TRUE)
  ) %>%
  ungroup()

save(bias_df, coverage_df, ci_env_df, file = "bias_coverage_df.RData")

