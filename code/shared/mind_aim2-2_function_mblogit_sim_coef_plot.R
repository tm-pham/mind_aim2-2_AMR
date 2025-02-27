# ============================================================================ #
# Function
# Project: MInD Aim 2.2
# Title: Plot of frequentist multinomial logistic regression (simulation study)
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ---------------------------------------------------------------------------- #
# This function plots the results from mblogit and saves the plots. 
# It is aimed at mblogit results from simulations studies. 
# ============================================================================ #
mblogit.sim.coef.plot <- function(model, 
                                  df_ab_table, 
                                  org_short = "SA",
                                  abx_order = c("beta_lactams", "fluoroquinolones", "lincosamides", "macrolides"),
                                  abx_plot_names = c("Anti-staphylococcal Beta-lactams", "Fluoroquinolones", 
                                                     "Lincosamides", "Macrolides"), 
                                  lookback = 14, 
                                  true_coef = 0.0025, 
                                  figure_path = "/Users/tm-pham/academia/hsph/mind/figures/", 
                                  data_path ="/Users/tm-pham/academia/hsph/mind/data/",
                                  imp_rate_flag = T, 
                                  suffix = "", 
                                  data_id = "", 
                                  plot_width = 16, 
                                  plot_height = 9, 
                                  plot_title = "", 
                                  y_axis_title = "terms"){
  library(ggplot2)
  source("/Users/tm-pham/academia/hsph/mind/publications/aim2-2/code/shared/plotting_template.R")
  tt <- broom::tidy(model, conf.int=TRUE)
  df_vars <- as.data.frame(do.call("rbind", stringr::str_split(tt$term, "~")))
  colnames(df_vars) <- c("antibiogram", "term")
  
  # Remove interaction terms
  ind <- which(!stringr::str_detect(df_vars$term, ":"))
  
  
  df_coef <- cbind(df_vars[ind,], tt[ind, -1]) %>% 
    mutate(term = stringr::str_remove_all(term, "^.*_prev_num_"), 
           term = stringr::str_remove_all(term, "class_"))
  
  
  # Variables for plotting 
  # All variables except time (i.e. including intercept and importation)
  vars_to_plot <- unique(df_coef$term[!stringr::str_detect(df_coef$term, "day|month|year|census_")])
  time_vars <- unique(df_coef$term[stringr::str_detect(df_coef$term, "day|month|year")])
  
  # Only variables related to past antibitic use
  vars_only_abx <- abx_order
  
  # Variables related to importation
  imp_var <- unique(vars_to_plot[stringr::str_detect(vars_to_plot, "imp")])
  antibiogram_match <- c(stringr::str_extract(imp_var, "AB_[0-9]"), stringr::str_extract(imp_var, "Other"))
  antibiogram_match <- antibiogram_match[!is.na(antibiogram_match)]
  
  
  imp_suffix <- ifelse(imp_rate_flag, "rate", "perc")
  if(length(imp_var)>1){
    if(stringr::str_detect(imp_var[1], "log")){
      imp_var_labels <- paste0("log(", ab_name_match, "_imp_", imp_suffix, ")")
    }else{
      imp_var_labels <- paste0(ab_name_match, "_imp_", imp_suffix)
    }
  }else{
    imp_var_labels <- paste0("Importation ", imp_suffix)
  }
  
  # Other covariates 
  other_covars <- unique(setdiff(df_coef$term, c(vars_only_abx, imp_var, time_vars)))
  
  # Data frames for the different plots
  df_coef_vars <- df_coef %>% 
    filter(term%in%vars_to_plot) %>% 
    mutate(term = factor(term, 
                         levels = c( "(Intercept)", imp_var, rev(abx_order)), 
                         labels = c("Intercept", imp_var, rev(abx_plot_names))))
  
  df_coef_time_vars <- df_coef %>% 
    filter(term%in%time_vars) 
  # Other covariates
  df_coef_other_covars <- df_coef %>% 
    filter(term%in%other_covars) 
  # Only antibiotic use
  df_coef_only_abx <- df_coef %>% 
    filter(term%in%vars_only_abx) %>% 
    mutate(term = factor(term, levels = rev(abx_order), labels=rev(abx_plot_names)))
  
  # Antibiotic use and importation terms
  if(length(imp_var)>0){
    df_coef_abx_imp <- df_coef %>% 
      filter(term%in%c(vars_only_abx, imp_var)) %>% 
      mutate(term = factor(term, levels = c(rev(abx_order), rev(imp_var)), labels= c(rev(abx_plot_names), rev(imp_var_labels))))
    # Only importation terms
    df_coef_imp <- df_coef %>% 
      filter(term%in%imp_var) %>% 
      mutate(term = factor(term, levels = c(rev(imp_var)), labels= c(rev(imp_var_labels))))
  }
  
  # All variables but no intercept
  df_coef_no_intercept <- df_coef %>% 
    filter(term%in%vars_to_plot, !term%in%c("(Intercept)")) %>% 
    mutate(term = factor(term, 
                         levels = c(rev(abx_order), setdiff(vars_to_plot, c(abx_order, "(Intercept)"))), 
                         labels = c(rev(abx_plot_names), setdiff(vars_to_plot, c(abx_order, "(Intercept)")))))
  
  
  coef_abx <- unique(df_vars$term[!stringr::str_detect(df_vars$term, "(Intercept)|importation|day|month|year")])
  # Determine the title for y-axis
  if(any(stringr::str_detect(coef_abx, "rate"))){
    y_axis_title <- paste0("Number of patient days receiving drug class\nin previous ", lookback, " days per 1,000 patient days")
  }
  
  if(any(stringr::str_detect(coef_abx, "prop"))){
    y_axis_title <- paste0("Proportion of patient days receiving drug class\nin previous ", lookback, " days")
  }
  
  # -------------------------------------------------------------------------- #
  # Intercept and abx plot
  coef_interval_plot_1 <- ggplot(df_coef_vars, 
                                 aes(x=estimate, y=term)) + 
    facet_wrap(~antibiogram, scales = "free_x") + 
    geom_vline(xintercept=0, linetype="dashed", color="gray") + 
    geom_pointrange(aes(xmin=conf.low, xmax = conf.high), color="blue", linewidth=1.5, size=1.5) + 
    labs(title = plot_title, x = "Log Odds Ratio (with respective to reference)", 
         y = y_axis_title) + 
    theme_template() + 
    theme(plot.title=element_text(hjust=-0.5))
  ggsave(coef_interval_plot_1, 
         file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_", suffix, ".pdf"), 
         width = plot_width, height = plot_height)
  
  # -------------------------------------------------------------------------- #
  # No intercept plot (but with importation terms)
  coef_interval_plot_2 <- ggplot(df_coef_no_intercept, 
                                 aes(x=estimate, y=term)) + 
    facet_wrap(~antibiogram, scales = "free_x") + 
    geom_vline(xintercept=0, linetype="dashed", color="gray") + 
    geom_pointrange(aes(xmin=conf.low, xmax = conf.high), color="blue", linewidth=1.5, size=1.5) + 
    labs(title = plot_title, 
         x = "Log Odds Ratio (with respective to reference)", 
         y = y_axis_title) + 
    theme_template() + 
    theme(plot.title=element_text(hjust=-0.5))
  ggsave(coef_interval_plot_2, 
         file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_noIntercept_", suffix, ".pdf"), 
         width = plot_width, height = plot_height)
  
  # -------------------------------------------------------------------------- #
  # Time variables plot
  coef_time_vars_plot <- ggplot(df_coef_time_vars, 
                                aes(x=estimate, y=term)) + 
    facet_wrap(~antibiogram, scales = "free_x") + 
    geom_vline(xintercept=0, linetype="dashed", color="gray") + 
    geom_pointrange(aes(xmin=conf.low, xmax = conf.high), color="blue", linewidth=1.5, size=1.5) + 
    labs(title = plot_title, 
         x = "Log Odds Ratio (with respective to reference)", 
         y = y_axis_title) + 
    theme_template() + 
    theme(plot.title=element_text(hjust=-0.5))
  ggsave(coef_time_vars_plot, 
         file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_time_vars_", suffix, ".pdf"), 
         width = plot_width, height = plot_height*1.5)
  
  # -------------------------------------------------------------------------- #
  # Other covariates plot
  coef_other_plot <- ggplot(df_coef_other_covars, 
                            aes(x=estimate, y=term)) + 
    facet_wrap(~antibiogram, scales = "free_x") + 
    geom_vline(xintercept=0, linetype="dashed", color="gray") + 
    geom_pointrange(aes(xmin=conf.low, xmax = conf.high), color="blue", linewidth=1.5, size=1.5) + 
    labs(title = plot_title, 
         x = "Log Odds Ratio (with respective to reference)", 
         y = y_axis_title) + 
    theme_template() + 
    theme(plot.title=element_text(hjust=-0.5))
  ggsave(coef_other_plot, 
         file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_other_covars_", suffix, ".pdf"), 
         width = plot_width, height = plot_height)
  
  # -------------------------------------------------------------------------- #
  # No intercept plot and no importation (only antibiotics)
  df_true_coef <- as.data.frame(cbind(terms = abx_plot_names, true_coef = as.numeric(true_coef))) %>% 
    mutate(true_coef = as.numeric(true_coef))
  
  # Apply the function to all antibiograms and bind them together
  df_true_coef <- lapply(unique(df_coef_only_abx$antibiogram), expand_antibiogram, abx_plot_names, true_coef) %>% 
    bind_rows(.id = "Antibiogram_ID")
  
  coef_interval_plot <- ggplot(df_coef_only_abx, 
                               aes(x=estimate, y=term)) + 
    facet_wrap(~antibiogram, scales = "free_x") + 
    geom_vline(xintercept=0, linetype="dashed", color="gray") + 
    geom_pointrange(aes(xmin=conf.low, xmax = conf.high), color="blue", linewidth=1.5, size=1.5) + 
    geom_point(data = df_true_coef, aes(x = true_coef, y = antibiotic), color = "red", shape = 18, size = 4) + 
    labs(title = plot_title, 
         x = "Log Odds Ratio (with respective to reference)", 
         y = y_axis_title) + 
    theme_template() + 
    theme(plot.title=element_text(hjust=-0.5))
  ggsave(coef_interval_plot, 
         file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_onlyAbx_", suffix, ".pdf"), 
         width = plot_width, height = plot_height)
  
  # -------------------------------------------------------------------------- #
  # Abx use and importation rate coefficient estimates
  if(length(imp_var)>0){
    coef_abx_imp_plot <- ggplot(df_coef_abx_imp, 
                                aes(x=estimate, y=term)) + 
      facet_wrap(~antibiogram, scales = "free_x") + 
      geom_vline(xintercept=0, linetype="dashed", color="gray") + 
      geom_pointrange(aes(xmin=conf.low, xmax = conf.high), color="blue", linewidth=1.5, size=1.5) + 
      labs(title = plot_title, 
           x = "Log Odds Ratio (with respective to reference)", 
           y = y_axis_title) + 
      theme_template() + 
      theme(plot.title=element_text(hjust=-0.5))
    ggsave(coef_abx_imp_plot, 
           file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_abx_imp_", suffix, ".pdf"), 
           width = plot_width, height = plot_height)
    # -------------------------------------------------------------------------- #
    # Only importation rate coefficient estimates
    coef_imp_plot <- ggplot(df_coef_imp, 
                            aes(x=estimate, y=term)) + 
      facet_wrap(~antibiogram, scales = "free_x") + 
      geom_vline(xintercept=0, linetype="dashed", color="gray") + 
      geom_pointrange(aes(xmin=conf.low, xmax = conf.high), color="blue", linewidth=1.5, size=1.5) + 
      labs(title = plot_title, 
           x = "Log Odds Ratio (with respective to reference)", 
           y = y_axis_title) + 
      theme_template() + 
      theme(plot.title=element_text(hjust=-0.5))
    ggsave(coef_imp_plot, 
           file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_imp_", suffix, ".pdf"), 
           width = plot_width, height = plot_height)
  }
  
  
  # -------------------------------------------------------------------------- #
  # Plot all parameters
  all_coef_interval_plot <- ggplot(df_coef, 
                                   aes(x=estimate, y=term)) + 
    facet_wrap(~antibiogram, scales = "free_x") + 
    geom_vline(xintercept=0, linetype="dashed", color="gray") + 
    geom_pointrange(aes(xmin=conf.low, xmax = conf.high), color="blue", linewidth=1.5, size=1.5) + 
    labs(title = plot_title, 
         x = "Log Odds Ratio (with respective to reference)", 
         y = y_axis_title) + 
    theme_template() + 
    theme(plot.title=element_text(hjust=-0.5))
  ggsave(all_coef_interval_plot, 
         file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_all_", suffix, ".pdf"), 
         width = plot_width, height = plot_height*3.5)
  
  # Save all plots in RData file
  save(coef_interval_plot, coef_interval_plot_1, coef_interval_plot_2, all_coef_interval_plot, 
       file = paste0(data_path, "mind_aim2-2_", data_id, "_mblogit_coef_plots_", suffix, ".RData"))
  
  # Return one of the plots
  return(coef_interval_plot)
}


# Create a function to expand each antibiogram into a table
expand_antibiogram <- function(antibiogram, antibiotics, true_coef, fill_susc = 0) {
  # Split the antibiogram into individual "S" or "R" values
  status <- strsplit(antibiogram, "-")[[1]]
  
  # Create a data frame with antibiotics, status, and true coefficients
  data.frame(
    antibiogram = antibiogram, 
    antibiotic = antibiotics,
    status = status,
    true_coef = ifelse(status == "R", true_coef, fill_susc)
  )
}


