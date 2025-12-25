# ============================================================================ #
# Project: MInD Aim 2.2
# PLOT
# Author: Thi Mui Pham, t.m.pham@hsph.harvard.edu
# Title: Figure 4B
# ------------------------------------------------------------------------------
# Description:
# FIGURE 4. Effect of recent antimicrobial prescribing on antimicrobial 
# resistance profiles in four target pathogens. 
# (B) Effects of recent 3/4GC (left) and BL/BLI (right) prescribing on E coli 
# (top) and K pneumoniae (bottom) AMR profiles.
# Patchwork data were used for all analyses, as described in the Appendix p. 41.
# AMR profiles used as outcome categories in the multilevel multinomial logistic 
# regression model are shown on the y-axis. The nested facet strip below the 
# pathogen header identifies the antimicrobial class used as the primary 
# exposure in the analysis. The x-axis shows the percent change in the odds of 
# each AMR profile associated with one additional treatment day per 100 
# patient-days of exposure to the specified antimicrobial class during the 
# preceding 14 days (expressed per 1000 patient-days), relative to isolates 
# susceptible to all key classes (S–S–S). Points indicate model estimates; 
# thick bars represent 80% confidence intervals and thin bars 95% intervals. 
# Grey backgrounds highlight pathogen–antimicrobial combinations for which an 
# association between prescribing and resistance was hypothesised. 
# ============================================================================ #
remove(list = ls())
# Load config file -------------------------------------------------------------
source(here::here("code", "00_config.R")) # phenotype_colors, bugs_ordered

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------
# Load data 
ec_data <- read.csv(paste0(RESULTS, "/figure3/mind_aim2-2_EC_6_mblogit_table_abx_results.csv"))
kp_data <- read.csv(paste0(RESULTS, "/figure3/mind_aim2-2_KP_5_mblogit_table_abx_results.csv"))

# Add pathogen identifier
ec_data$pathogen <- "Escherichia coli"
kp_data$pathogen <- "Klebsiella pneumoniae"

df <- rbind(ec_data, kp_data) %>% 
  filter(term %in% c("3/4GC", "BL/BLI")) %>% 
  left_join(facet_table, by = c("pathogen" = "organismofinterest")) %>% 
  mutate(org_ab = factor(paste0(pathogen, "__", ab_name), 
                         levels = c("Escherichia coli__R-R-R",
                                    "Escherichia coli__R-R-S",
                                    "Escherichia coli__S-R-R",
                                    "Escherichia coli__R-S-R",
                                    "Escherichia coli__R-S-S",
                                    "Escherichia coli__S-R-S",
                                    "Escherichia coli__S-S-R", 
                                    "Klebsiella pneumoniae__R-R-R",
                                    "Klebsiella pneumoniae__S-R-R",
                                    "Klebsiella pneumoniae__R-S-R",
                                    "Klebsiella pneumoniae__S-S-R")), 
         ab_name = factor(ab_name, levels = rev(GC_phenotype_order)), 
         term = factor(term, levels = c("3/4GC", "BL/BLI")), 
         organismofinterest = factor(pathogen, levels = c("Escherichia coli", "Klebsiella pneumoniae"))) %>% 
  group_by(pathogen) %>% 
  arrange(pathogen, term, org_ab) %>%
  mutate(y_pos = row_number()) %>% 
  ungroup()

EC_GC_shade_data <- df %>%
  filter(cph34 ==1, term =="3/4GC", pathogen == "Escherichia coli") %>% 
  mutate(y_pos = y_pos + 3)

EC_BLI_shade_data <- df %>%
  filter(b_lac ==1, term == "BL/BLI", pathogen == "Escherichia coli") %>% 
  mutate(y_pos = y_pos-6)

KP_GC_shade_data <- df %>%
  filter(cph34 ==1, term =="3/4GC", pathogen == "Klebsiella pneumoniae") %>% 
  mutate(y_pos = y_pos + 2)

KP_BLI_shade_data <- df %>%
  filter(b_lac ==1, term == "BL/BLI", pathogen == "Klebsiella pneumoniae") %>% 
  mutate(y_pos = y_pos-5)


# ------------------------------------------------------------------------------
# PLOT
# ------------------------------------------------------------------------------
(figure4B <- ggplot(df, aes(y = ab_name, x = estimate_interp)) +
   # Facet by pathogen
   facet_nested_wrap(~facet_labs + term, scales = "free") + 
   # Grey background for corresponding profiles (GC and BL/BLI-resistant)-------
   geom_rect(data = EC_GC_shade_data,
           aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5), xmin = -Inf, xmax = Inf,
           fill = "grey85", inherit.aes = FALSE) +
   geom_rect(data = EC_BLI_shade_data,
             aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5), xmin = -Inf, xmax = Inf,
             fill = "grey85", inherit.aes = FALSE) +
   geom_rect(data = KP_GC_shade_data,
             aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5), xmin = -Inf, xmax = Inf,
             fill = "grey85", inherit.aes = FALSE) +
   geom_rect(data = KP_BLI_shade_data,
             aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5),
             xmin = -Inf, xmax = Inf, fill = "grey85", inherit.aes = FALSE) +
   # Vertical reference line at 0 -----------------------------------------------
 geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
   # 95% CI - thin line ---------------------------------------------------------
 geom_segment(aes(x = conf2.5_interp, xend = conf97.5_interp, 
                  y = ab_name, 
                  yend = ab_name),
              linewidth = 0.7, color = "gray30") +
   # 80% CI - thick line --------------------------------------------------------
 geom_segment(aes(x = conf10_interp, xend = conf90_interp,
                  y = ab_name,
                  yend = ab_name, color = org_ab),
              linewidth = 2.8) +
   scale_color_manual(values = phenotype_colors[9:24]) +
   # Point estimate -------------------------------------------------------------
 geom_point(size = 3.5, shape = 21, fill = "white", color = "black") +
   # Add FQL label and resistance pattern on left -------------------------------
 geom_text(aes(x = -Inf, label = paste0("FQL\n", ab_name)),
           hjust = 1, vjust = 0.5, size = 2.8, lineheight = 0.85) +
   # Add estimate on right
   geom_text(aes(x = Inf, 
                 label = sprintf("%.1f%% (%.1f%%, %.1f%%)\n(compared to S-S-S)", 
                                 estimate_interp, conf2.5_interp, conf97.5_interp)),
             hjust = -0.05, vjust = 0.5, size = 2.5, lineheight = 0.85) +
   # Labels
   labs(
     x = "Increase in odds\n(%, compared to S-S-S)",
     y = "AMR profiles\n"
   ) +
   
   # Clean theme matching reference
   theme_minimal(base_size = 10) +
   theme(
     legend.position = "none",
     plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15)),
     axis.title.x = element_blank(), 
     axis.title.y = element_text(face = "plain", size = 16, margin = margin(b = 12)),
     axis.text.x = element_text(size = 14),
     axis.text.y = element_text(size = 14),
     axis.ticks.y = element_blank(),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     strip.text.x      = ggtext::element_markdown(size = 16), 
     strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.3),
     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
     panel.spacing.x = unit(1.2, "lines"),
     panel.spacing.y = unit(1, "lines"),
     plot.margin = margin(10, 20, 10, 10),
     plot.background = element_rect(fill = "transparent", color = NA), 
     panel.background = element_rect(fill = "transparent", color = "NA")
   ))

# Save figure ------------------------------------------------------------------
saveRDS(figure4B, paste0(RESULTS, "/figure4/figure4B_Enterobacterales.rds"))

