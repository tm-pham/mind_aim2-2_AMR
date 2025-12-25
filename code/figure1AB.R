# ============================================================================ #
# Project: MInD Aim 2.2
# Title: Figure 1 for manuscript
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# ---------------------------------------------------------------------------- #
# A: Barplot of phenotype incidence 
# B: GEE time trend results for phenotype incidence
# ============================================================================ #
remove(list = ls())
################################################################################
# Load config file
source(here::here("00_config.R")) # phenotype_colors, bugs_ordered

################################################################################
# Figure A: Antibiogram incidence plot
load(here::here("results/figure1", "mind_aim2-2_phenotype_inc_data.RData"))
figure1A <- readRDS(paste0(RESULTS, "/figure1/figure1A_phenotype_trends.rds"))

################################################################################
# Figure B: Plot of GEE time trend results
# Load data --------------------------------------------------------------------
load(here::here("results", "figure1", "mind_aim2-2_phenotype_gee_time_trend_results.RData")) # df_time_trend

# Data preparation -------------------------------------------------------------
df_gee <- bind_rows(df_time_trend, .id = "organismofinterest") %>% 
  mutate(time_period = factor(time_period, levels = c("2020-2021", "2007-2019", "2007-2021")), 
         organismofinterest = factor(organismofinterest, levels = bugs_ordered), 
         org_ab = paste0(organismofinterest, "__", antibiogram)) %>% 
  left_join(df_ab_year %>% ungroup() %>% select(organismofinterest, ab_name, facet_lab) %>% unique(), 
            by = c("organismofinterest", "antibiogram"="ab_name")) %>%
  mutate(facet_lab = factor(facet_lab, levels = facet_labs), 
           trend_direction = case_when(
             lower.CL >= 0 | time.trend > 0 ~ "Increasing",
             upper.CL < 0 ~ "Decreasing",
             TRUE ~ "Stable"))

# Background data: one row per facet panel
facet_backgrounds <- df_gee %>%
  group_by(antibiogram, facet_lab) %>%
  filter(n() > 0) %>%     # ensures panel has data
  distinct(trend_direction, .keep_all = TRUE)

panel_info <- df_gee %>%
  count(antibiogram, facet_lab, name = "n_obs") %>%
  mutate(has_data = n_obs > 0)

################################################################################
# PLOT
(figure1B <- ggplot(df_gee, 
                   aes(x=as.numeric(time.trend), y = 1, color = org_ab)) + 
    # manual borders on data-containing panels only
    geom_rect(
      data = panel_info %>% filter(has_data),
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      fill = "white",
      color = "black",
      linewidth = 0.7,
      inherit.aes = FALSE
    ) + 
    ggh4x::facet_nested(rows = vars(antibiogram), 
                        cols = vars(facet_lab), 
                        scales = "free", 
                        render_empty = FALSE, 
                        switch = "y") + 
   # Background rectangles
   geom_rect(
     data = facet_backgrounds,
     aes(fill = trend_direction),
     xmin = -Inf, xmax = Inf,
     ymin = -Inf, ymax = Inf,
     alpha = 0.3,
     inherit.aes = FALSE
   ) +
   scale_fill_manual(values = trend_colors, guide = "none") +
    geom_vline(xintercept = 0, linetype = 1, linewidth = 0.8, color = "grey20") + 
    geom_errorbar(aes(xmin=as.numeric(lower.CL), xmax=as.numeric(upper.CL)), width = 0., linewidth = 1.5) + 
    geom_point(fill = "white", stroke = 2, size = 5, shape =21) + 
    labs(x = "Average annual percentage change (%)") + 
    scale_color_manual(values = phenotype_colors) +
    theme_template() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=20, face = "plain"),
          axis.text.x = element_text(size=20),
          legend.position = "none", 
          plot.tag = element_text(size=20),
          strip.text        = element_text(face = "plain"),
          strip.text.x = ggtext::element_markdown(size = 21), 
          strip.text.y = element_text(face = "bold", size = 21),
          panel.spacing.x = unit(1.3, "cm"),
          strip.background  = element_rect(fill = "lightgrey", color = "black", linewidth = 1.2),
          panel.spacing.y   = unit(0.5, "cm"),
          panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank(), 
          panel.background   = element_blank(),
          panel.border = element_blank()
          ))


################################################################################
# Combine Figure A and B
(figure <- (figure1A/figure1B) +  
   plot_layout(nrow = 2, heights = c(1.3, 1)) + 
   plot_annotation(tag_level = 'A') &
   theme(plot.tag = element_text(size = 24, face = "plain")))

ggsave(figure, file = paste0(FIGURES, "/figure1/figure1AB.pdf"), 
       width = 19, height = 26)
