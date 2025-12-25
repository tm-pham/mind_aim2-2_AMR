# ============================================================================ #
# Project: MInD Aim 2.2
# PLOT
# Author: Thi Mui Pham, t.m.pham@hsph.harvard.edu
# Title: Figure 4C
# Description: multinomial logit model forest plot for P aeruginosa
# ============================================================================ #

################################################################################
# Load config file 
source(here::here("00_config.R")) # phenotype_colors, bugs_ordered

################################################################################
# Load data 
pa_data <- read.csv(paste0(RESULTS, "/figure3/mind_aim2-2_PA_5_mblogit_table_abx_results.csv"))

# Add pathogen identifier
pa_data$pathogen <- "Pseudomonas aeruginosa"

df_pa <- pa_data %>% 
  filter(term %in% c("BL/BLI", "CPM")) %>% 
  left_join(facet_table, by = c("pathogen" = "organismofinterest")) %>% 
  mutate(org_ab = factor(paste0(pathogen, "__", ab_name), 
                         levels = rev(paste0("Pseudomonas aeruginosa__", CPM_phenotype_order))), 
         ab_name = factor(ab_name, levels = rev(CPM_phenotype_order))) %>% 
  group_by(pathogen) %>% 
  arrange(pathogen, term, org_ab) %>%
  mutate(y_pos = row_number()) %>% 
  ungroup()

PA_BL_shade_data <- df_pa %>%
  filter(b_lac ==1, term =="BL/BLI", pathogen == "Pseudomonas aeruginosa") %>% 
  mutate(y_pos = y_pos)

PA_CPM_shade_data <- df_pa %>%
  filter(antipseudomonalcpm ==1, term =="CPM", pathogen == "Pseudomonas aeruginosa") %>% 
  mutate(y_pos = y_pos -6)


################################################################################
# FIGURE 4C
################################################################################
# Create plot
(figure4C <- ggplot(df_pa, aes(y = ab_name, x = estimate_interp)) +
   # Facet by pathogen
   facet_nested(~facet_labs + term, scales = "free_y") + 
   geom_rect(data = PA_BL_shade_data,
             aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5), xmin = -Inf, xmax = Inf,
             fill = "grey85", inherit.aes = FALSE) +
   geom_rect(data = PA_CPM_shade_data,
             aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5), xmin = -Inf, xmax = Inf,
             fill = "grey85", inherit.aes = FALSE) +
   # Vertical reference line at 0
   geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
   
   # 95% CI - thin line
   geom_segment(aes(x = conf2.5_interp, xend = conf97.5_interp, 
                    y = ab_name, 
                    yend = ab_name),
                linewidth = 0.8, color = "gray30") +
   
   # 80% CI - thick line
   geom_segment(aes(x = conf10_interp, xend = conf90_interp,
                    y = ab_name,
                    yend = ab_name, color = org_ab),
                linewidth = 2.5) +
   scale_color_manual(values = phenotype_colors[25:32]) +
   
   # Point estimate
   geom_point(size = 3.5, shape = 21, fill = "white", color = "black") +
   
   # Add FQL label and resistance pattern on left
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
     axis.title.x = element_text(face = "plain", size = 16, margin = margin(t = 8)),
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
     # plot.margin = margin(10, 20, 10, 10),
     plot.background = element_rect(fill = "transparent", color = NA), 
     panel.background = element_rect(color = "white")
   ))

################################################################################
# Combined figure 4
################################################################################
# Load figure 4A and 4B
figure4A <- readRDS(paste0(RESULTS, "/figure4/figure4A_SA.rds"))
figure4B <- readRDS(paste0(RESULTS, "/figure4/figure4B_Enterobacterales.rds"))

(figure4 <- (figure4A/figure4B/figure4C) + 
    plot_layout(nrow = 3, heights = c(0.7, 2, 0.9)) + 
    plot_annotation(tag_level = 'A') & 
   theme(plot.tag = element_text(face = "plain", size = 16)))
ggsave(figure4, file = paste0(FIGURES, "/figure4/figure4.pdf"), 
       width = 7, height = 13)



