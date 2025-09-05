# ============================================================================ #
# Project: MInD Aim 2.2
# FUNCTION
# Title: Plot of frequentist multinomial logistic regression for manuscript
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# ---------------------------------------------------------------------------- #
# This function plots the results from mblogit and saves the plots. 
# ============================================================================ #
mblogit.coef.plot <- function(model, 
                              org_short = "SA",
                              abx_order = c("antistaphbetalactams", 
                                            "flq", 
                                            "linc", 
                                            "macr"),
                              abx_plot_names = c("Anti-staphylococcal Beta-lactams", 
                                                 "Fluoroquinolones", 
                                                 "Lincosamides", 
                                                 "Macrolides"), 
                              lookback = 14, 
                              figure_path = "path/to/folder/figures/", 
                              data_path ="path/to/folder/results/",
                              df_perc = NULL,
                              ab_name_order = NULL,
                              suffix = "", 
                              exponentiate = T,
                              outcome_option = 1, # 1 = increase of 14 treatment days per 1,000 pt days, 2 = increase of treatment day per 100 pt-days
                              round_digits = 2, 
                              plot_width = 16, 
                              plot_height = 9, 
                              plot_title = "", 
                              y_title = "terms"){
  source("path_to_folder/code/mind_aim2-2_function_assign_codes.R")
  
  library(ggplot2)
  library(ggtext)
  
  # -------------------------------------------------------------------------- #
  # Function to create color coded facet labels
  # Color is hard-coded 
  label_fun <- function(ab1, ab2, ab3, ab4=NULL){
    res <- paste0(
      ifelse(ab1 == 1, "<span style='color:violetred4;'>R</span>", "S"), "-",
      ifelse(ab2 == 1, "<span style='color:violetred4;'>R</span>", "S"), "-",
      ifelse(ab3 == 1, "<span style='color:violetred4;'>R</span>", "S")
    )
    
    if(!is.null(ab4)){
      res <- paste0(res, ifelse(ab4 == 1, "-<span style='color:violetred4;'>R</span>", "-S"))
    }
    return(res)
  }
  
  # Table for phenotype names 
  df_ab_table = model$data %>% select(antibiogram_name, ab_name) %>% unique() 
  # Determine reference category for plot title
  antibiogram_reference <- df_ab_table[df_ab_table$antibiogram_name==levels(model$data$antibiogram_name)[1], "ab_name"]
  
  
  ##############################################################################
  # Data frame with estimate and confidence intervals
  # Transform to odds ratio (instead of log-odds ratio)
  df_confint <- cbind(summary(model)$coefficients[, "Estimate"], 
                      confint(model, level = 0.8), 
                      confint(model, level = 0.95))
  if(exponentiate){
    df_confint <- cbind(exp(df_confint), 
                        summary(model)$coefficients[, 4])
  }
  colnames(df_confint) <- c("estimate","conf10", "conf90", "conf2-5", "conf97-5", "p-value")
  # Make rownames to a column
  df_confint <- as.data.frame(df_confint) %>% 
    mutate(term = rownames(.), .before=1)
  rownames(df_confint) <- NULL

  # Table of variables for each antibiogram
  df_vars <- as.data.frame(do.call("rbind", stringr::str_split(df_confint$term, "~")))
  colnames(df_vars) <- c("antibiogram", "term")
  
  # Remove interaction terms
  ind <- which(!stringr::str_detect(df_vars$term, ":"))
  
  # Data frame for color coding 
  df_codes <- as.data.frame(do.call(rbind, lapply(df_ab_table$ab_name, assign_codes)))
  colnames(df_codes) <- abx_order[1:ncol(df_codes)]
  if(length(abx_order)>ncol(df_codes)){
    for(i in (ncol(df_codes)+1):length(abx_order)){
      df_codes[, abx_order[i]] <- rep(0, nrow(df_codes))
    }
  }
  df_codes <- cbind(ab_name = df_ab_table$ab_name, df_codes)
  rownames(df_codes) <- NULL

  # Add color codes to antibiogram name table
  df_ab_table <- left_join(df_ab_table, df_codes)
  
  # Make data frame for plotting
  df_coef <- cbind(df_vars[ind,], df_confint[ind, -1]) %>% 
    left_join(df_ab_table, 
              by = c("antibiogram" = "antibiogram_name")) %>% 
    mutate(term = stringr::str_remove_all(term, "^.*_prev_num_"), 
           term = stringr::str_remove_all(term, "class_|subclass_"))
  print(head(df_coef))
  
  
  df_coef$flag <- ifelse(df_coef$term==abx_order[1], df_coef[,abx_order[1]], 0)
  for(i in 1:(length(abx_order))){
    df_coef$flag <- ifelse(df_coef$term==abx_order[i], df_coef[,abx_order[i]], df_coef$flag)
  }  
  
  ##############################################################################
  # Variables for plotting 
  vars_to_plot <- unique(df_coef$term)
  # Only variables related to past antibitic use
  vars_only_abx <- unique(vars_to_plot[!stringr::str_detect(vars_to_plot, 
                                                            "(Intercept)|importation|facility|median|perc_|Score|complexitylevel|prop_|day|month|year|census|_rate")])
  # Variables related to importation
  imp_var <- unique(vars_to_plot[stringr::str_detect(vars_to_plot, "importation")])
  
  # Interaction variables
  varsX <- unique(vars_to_plot[stringr::str_detect(vars_to_plot, ":")])
  
  ##############################################################################
  # Data frames for the different plots
  df_imp <- df_coef %>% 
    filter(term%in%imp_var) %>% 
    mutate(term_ab = stringr::str_remove(term, "_importation_rate")) %>% 
    left_join(df_ab_table %>% select(antibiogram_name, ab_name) %>% rename(term_ab_name = ab_name), 
              by = c("term_ab" = "antibiogram_name")) %>% 
    mutate(flag = ifelse(ab_name == term_ab_name, 1, 0))
  write.csv(df_imp, 
            file = paste0(data_path, "/mind_aim2-2_", data_id, "_mblogit_table_importation_", suffix, "_odds.csv"), 
            row.names = F)
  
  if(outcome_option==1){
    df_coef_only_abx <- df_coef %>% 
      filter(term%in%vars_only_abx) %>% 
      mutate(term = factor(term, levels = rev(abx_order), labels=rev(abx_plot_names))) %>% 
      mutate(across(.cols=c("estimate", "conf10", "conf90", "conf2-5", "conf97-5"), 
                    .fns = ~round((.x^14-1),round_digits), 
                    .names = "{.col}_interp"))
    write.csv(df_coef_only_abx, 
              file = paste0(data_path, "/mind_aim2-2_", data_id, "_mblogit_table_onlyAbx_", suffix, "_14dodds.csv"), 
              row.names = F)
    
    df_imp <- df_imp %>% 
      mutate(across(.cols=c("estimate", "conf10", "conf90", "conf2-5", "conf97-5"), 
                    .fns = ~round((.x^14-1),round_digits), 
                    .names = "{.col}_interp"))
    
    write.csv(df_imp, 
              file = paste0(data_path, "/mind_aim2-2_", data_id, "_mblogit_table_importation_", suffix, "_14dodds.csv"), 
              row.names = F)
  }
  
  if(outcome_option==2){
    df_coef_only_abx <- df_coef %>% 
      filter(term%in%vars_only_abx) %>% 
      mutate(term = factor(term, levels = rev(abx_order), labels=rev(abx_plot_names))) %>% 
      mutate(across(.cols=c("estimate", "conf10", "conf90", "conf2-5", "conf97-5"), 
                    .fns = ~round((.x^10-1), round_digits), 
                    .names = "{.col}_interp"))
    write.csv(df_coef_only_abx %>% 
                mutate(across(.cols=c("estimate_interp", "conf10_interp", "conf90_interp", "conf2-5_interp", "conf97-5_interp"), 
                              .fns = ~round(100*.x, 1))), 
              file = paste0(data_path, "/mind_aim2-2_", data_id, "_mblogit_table_onlyAbx_", suffix, "oddsPer100.csv"), 
              row.names = F)
    
    df_imp <- df_imp %>% 
      mutate(across(.cols=c("estimate", "conf10", "conf90", "conf2-5", "conf97-5"), 
                    .fns = ~round(.x^10-1,round_digits), 
                    .names = "{.col}_interp"))
    
    write.csv(df_imp %>% 
                mutate(across(.cols=c("estimate_interp", "conf10_interp", "conf90_interp", "conf2-5_interp", "conf97-5_interp"), 
                              .fns = ~round(100*.x, 1))), 
              file = paste0(data_path, "/mind_aim2-2_", data_id, "_mblogit_table_importation_", suffix, "oddsPer100.csv"), 
              row.names = F)
  }
  
  
  if(!is.null(ab_name_order)){
    df_coef_only_abx <- df_coef_only_abx %>% 
      mutate(ab_name = factor(ab_name, levels = ab_name_order))
    
    df_imp <- df_imp %>% 
      mutate(ab_name = factor(ab_name, levels = ab_name_order))
  }
  
  ##############################################################################
  # Create color coded facet texts
  if(length(abx_order)==3){
    df_coef_only_abx$ab_label <- mapply(label_fun, df_coef_only_abx[,abx_order[1]], df_coef_only_abx[,abx_order[2]], df_coef_only_abx[, abx_order[3]])
    df_imp$ab_label <-  mapply(label_fun, df_imp[,abx_order[1]], df_imp[,abx_order[2]], df_imp[,abx_order[3]])
  }
  
  if(length(abx_order)==4){
    df_coef_only_abx$ab_label <- mapply(label_fun, df_coef_only_abx[,abx_order[1]], df_coef_only_abx[,abx_order[2]], df_coef_only_abx[, abx_order[3]], df_coef_only_abx[, abx_order[4]])
    df_imp$ab_label <-  mapply(label_fun, df_imp[,abx_order[1]], df_imp[,abx_order[2]], df_imp[,abx_order[3]], df_imp[,abx_order[4]])
  }

  ab_label_order <- df_coef_only_abx %>% select(ab_name, ab_label) %>% unique() %>% arrange(ab_name) %>% select(ab_label) %>% unlist() %>% unique()
  
  df_coef_only_abx$ab_label <- factor(df_coef_only_abx$ab_label, levels = ab_label_order)
  
  df_imp$ab_label <- factor(df_imp$ab_label, levels = ab_label_order)

  
  ##############################################################################
  # Titles for axes
  # y_axis_title <- paste0("Number of patient days receiving drug class\nin previous ", lookback, " days per 1,000 patient days")
  y_axis_title <- "Recent facility-level antibiotic prescribing"
  x_title <- paste0("Odds Ratio\n(compared to ", antibiogram_reference, ")")
  x_axis_title <- ifelse(!exponentiate, paste0("Log ", x_title), x_title)
  x_axis_title_interp <- "Increase in odds\n(compared to S-S-S)"
  
  plot_title <- ""
  
  ##############################################################################
  # Plot of only previous antibiotic use coefficients
  (plot1 <- ggplot(df_coef_only_abx, 
                   aes(x=estimate, y=term)) + 
     facet_wrap(~ab_label, labeller = label_value, scales = "free_y", ncol = 1) +
     geom_vline(xintercept=1, linetype="dashed", color="black", linewidth = 1.3) + 
     geom_pointrange(aes(xmin=`conf2-5`, xmax = `conf97-5`), linewidth=1.5, size=1.5, color ="black") +
     geom_pointrange(aes(xmin=conf10, xmax = conf90, color = factor(flag)), 
                     linewidth=2.2, size=1.5, 
                     fill = "white", shape = 21) +
     labs(title = plot_title, 
          x = x_axis_title, 
          y = y_axis_title) + 
     scale_color_manual(values = c("darkgrey", "violetred4")) +
     theme_template() + 
     theme(plot.title=element_text(hjust=-0.5), 
           strip.text = element_markdown(size=28, face = "bold"),
           axis.text.y = element_text(size = 24),
           axis.title.y = element_text(size = 26),
           axis.text.x = element_text(size = 22),
           axis.title.x = element_text(size = 24),
           legend.position = "none"))
  ggsave(plot1, 
         file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_onlyAbx_", suffix, "odds_long.pdf"), 
         width = plot_width, height = plot_height)
  
  ##############################################################################
  # Plot of only previous antibiotic use coefficients
  # Interpreted 
  (plot1.2 <- ggplot(df_coef_only_abx, 
                   aes(x=estimate_interp, y=term)) + 
     facet_wrap(~ab_label, labeller = label_value, scales = "free_y", ncol = 1) +
     geom_vline(xintercept=0, linetype="dashed", color="black", linewidth = 1.3) + 
     geom_pointrange(aes(xmin=`conf2-5_interp`, xmax = `conf97-5_interp`), linewidth=1.5, size=1.5, color ="black") +
     geom_pointrange(aes(xmin=conf10_interp, xmax = conf90_interp, color = factor(flag)), 
                     linewidth=2.2, size=1.5, 
                     fill = "white", shape = 21) +
     labs(title = plot_title, 
          x = x_axis_title_interp, 
          y = y_axis_title) + 
     scale_color_manual(values = c("darkgrey", "violetred4")) +
     scale_x_continuous(labels = scales::percent) + 
     theme_template() + 
     theme(plot.title=element_text(hjust=-0.5), 
           plot.margin = margin(0, 1, 0.1, 0.5, "cm"),
           strip.text = element_markdown(size=28, face = "bold"),
           axis.text.y = element_text(size = 24),
           axis.title.y = element_text(size = 26),
           axis.text.x = element_text(size = 22),
           axis.title.x = element_text(size = 24),
           legend.position = "none"))
  ggsave(plot1.2, 
         file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_onlyAbx_", suffix, "odds_long_interp.pdf"), 
         width = plot_width, height = plot_height)

  
  ##############################################################################
  # Plot of only previous antibiotic use coefficients 
  # 1% increase in antibiotic prescribing
  if(!is.null(df_perc)){
    df_coef_only_abx <- df_coef_only_abx %>% 
      left_join(df_perc, by = c("term"="variable")) %>% 
      mutate_at(c("estimate_interp", "conf10_interp", "conf90_interp", "conf2-5_interp", "conf97-5_interp"),
                ~./(100*perc)) %>% 
      mutate(term = factor(term, levels = rev(abx_plot_names)))
    
    (plot1.2 <- ggplot(df_coef_only_abx, 
                       aes(x=estimate_interp, y=term)) + 
        facet_wrap(~ab_label, labeller = label_value, scales = "free_y", ncol = 1) + 
        geom_vline(xintercept=0, linetype="dashed", color="black", linewidth = 1.3) + 
        geom_pointrange(aes(xmin=`conf2-5_interp`, xmax = `conf97-5_interp`), linewidth=1.5, size=1.5, color ="black") +
        geom_pointrange(aes(xmin=conf10_interp, xmax = conf90_interp, color = factor(flag)), 
                        linewidth=2.2, size=1.5, 
                        fill = "white", shape = 21) +
        labs(title = plot_title, 
             x = x_axis_title_interp, 
             y = y_axis_title) + 
        scale_color_manual(values = c("darkgrey", "violetred4")) +
        scale_x_continuous(labels = scales::percent) + 
        theme_template() + 
        theme(plot.title=element_text(hjust=-0.5), 
              plot.margin = margin(0, 1, 0.1, 0.5, "cm"),
              strip.text = element_markdown(size=28, face = "bold"),
              axis.text.y = element_text(size = 24),
              axis.title.y = element_text(size = 26),
              axis.text.x = element_text(size = 21),
              axis.title.x = element_text(size = 24),
              legend.position = "none"))
    ggsave(plot1.2, 
           file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_onlyAbx_perc_", suffix, "odds_long.pdf"), 
           width = plot_width+0.8, height = plot_height)
    
  }
  
  
  ##############################################################################
  # Plot of community importation coefficients
  (plot2 <- ggplot(df_imp[complete.cases(df_imp),], 
                   aes(x=estimate_interp, y=term_ab_name)) + 
     facet_wrap(~ab_label, labeller = label_value, scales = "free_y", ncol = 1) +
     geom_vline(xintercept=0, linetype="dashed", color="black", linewidth = 1.3) + 
     geom_pointrange(aes(xmin=`conf2-5_interp`, xmax = `conf97-5_interp`), linewidth=1.5, size=1.5, color ="black") +
     geom_pointrange(aes(xmin=conf10_interp, xmax = conf90_interp, color = factor(flag)), 
                     linewidth=2.2, size=1.5, 
                     fill = "white", shape = 21) +
     labs(title = plot_title, 
          x = x_axis_title_interp, 
          y = "Number of patient days with phenotype collected within 3 days of admission \nduring 90 days before isolate was collected (per 1,000 community-onset patient days)") + 
     scale_color_manual(values = c("darkgrey", "violetred4")) +
     scale_x_continuous(labels = scales::percent) + 
     theme_template() + 
     theme(plot.title=element_text(hjust=-0.5), 
           plot.margin = margin(0, 1, 0.1, 0.5, "cm"),
           strip.text = element_markdown(size=28, face = "bold"),
           axis.text.y = element_text(size = 24),
           axis.title.y = element_text(size = 26),
           axis.text.x = element_text(size = 22),
           axis.title.x = element_text(size = 24),
           legend.position = "none"))
  ggsave(plot2, 
         file = paste0(figure_path, "/mind_aim2-2_", data_id, "_mblogit_importation_terms_", suffix, "odds_long.pdf"), 
         width = 7, height = 18)
  

  return(plot1)
}
