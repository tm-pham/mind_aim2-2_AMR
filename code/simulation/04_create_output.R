# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Create summary statistics and plots of simulation results
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ---------------------------------------------------------------------------- #
# Bias
# Coverage
# Confidence Interval Width and Envelope
# ============================================================================ #

# Bias -------------------------------------------------------------------------
# Defined as the difference between the estimated and true coefficients
# If bias is positive, 
bias_df <- df_combined %>%
  filter(flag == 1) %>% 
  left_join(coef_combined, by = c("sim_id","term")) %>%
  mutate(term = factor(term, levels = abx_order, labels = abx_labels)) %>% 
  group_by(sim_id, antibiogram, term) %>% 
  summarize(bias = estimate - true_coef)
print(bias_df)

# Coverage ---------------------------------------------------------------------
coverage_df <- df_combined %>%
  filter(flag == 1) %>% 
  left_join(coef_combined, by = c("sim_id","term")) %>%
  mutate(term = factor(term, levels = abx_order, labels = abx_labels)) %>% 
  group_by(antibiogram, term) %>%
  summarize(
    coverage = sum(conf.low <= true_coef & conf.high >= true_coef, na.rm = TRUE)/length(unique(df_combined$sim_id)),
    .groups  = "drop"
  )
print(coverage_df)

# Confidence Interval Width and Envelope --------------------------------------
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

# Plot -------------------------------------------------------------------------
mblogit.sim.coef.plot(model = model, 
                      true_coef = params$gamma,
                      abx_order = abx_order, 
                      abx_plot_names = abx_labels, 
                      figure_path = FIGURES, 
                      data_path ="", 
                      save_RData = F)