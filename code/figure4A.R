# ============================================================================ #
# Project: MInD Aim 2.2
# PLOT
# Author: Thi Mui Pham, t.m.pham@hsph.harvard.edu
# Title: Figure 4A
# Dewscription: multinomial logit model forest plot for S aureus
# ============================================================================ #

################################################################################
# Load config file
source(here::here("00_config.R")) # phenotype_colors, bugs_ordered

################################################################################
# Load data 
sa_data <- read.csv(paste0(RESULTS, "/figure3/mind_aim2-2_SA_4_mblogit_table_abx_results.csv"))

# Add pathogen identifier
sa_data$pathogen <- "Staphylococcus aureus"

df_sa <- sa_data %>% 
  filter(term %in% c("ASBL", "MACR")) %>% 
  left_join(facet_table, by = c("pathogen" = "organismofinterest")) %>% 
  mutate(org_ab = factor(paste0(pathogen, "__", ab_name), 
                         levels = rev(paste0("Staphylococcus aureus__", ASBL_phenotype_order))), 
         ab_name = factor(ab_name, levels = rev(ASBL_phenotype_order))) %>% 
  group_by(pathogen) %>% 
  arrange(term, org_ab) %>%
  mutate(y_pos = row_number()) %>% 
  ungroup()

ASBL_shade_data <- df_sa %>%
  filter(ASBL ==1, term =="ASBL") %>% 
  mutate(y_pos = y_pos)

MACR_shade_data <- df_sa %>%
  filter(macr ==1, term == "MACR") %>% 
  mutate(y_pos = y_pos-5)

################################################################################
# PLOT
################################################################################
# Create plot
(figure4A <- ggplot(df_sa, aes(y = ab_name, x = estimate_interp)) +
   # Facet by pathogen (1 column)
   # facet_wrap(~term, ) + 
   facet_nested(~facet_labs + term, scales = "free_y") + 
   # Grey background for R-* phenotypes
   geom_rect(
     data = ASBL_shade_data,
     aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5),
     xmin = -Inf, xmax = Inf,
     fill = "grey85",
     inherit.aes = FALSE
   ) +
   geom_rect(
     data = MACR_shade_data,
     aes(ymin = y_pos - 0.5, ymax = y_pos + 0.5),
     xmin = -Inf, xmax = Inf,
     fill = "grey85",
     inherit.aes = FALSE
   ) +
   
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
   scale_color_manual(values = phenotype_colors[1:8]) +
   
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
     panel.background = element_rect(color = "white")
   ))
saveRDS(figure4A, paste0(RESULTS, "/figure4/figure4A_SA.rds"))


