# ============================================================================ #
# Project: MInD Aim 2.2
# PLOT
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# Title: Figure 1A in manuscript (AMR profile incidence trends)
# ---------------------------------------------------------------------------- #
# Description: 
# Trends in antimicrobial resistance profiles among hospital-onset isolates of 
# four target pathogens in the U.S. Veterans Affairs Healthcare Administration, 
# February 1, 2007–December 31, 2021. 
# Patchwork data were used for analyses, as described in the Appendix (p.41). 
# Each column represents one of four target pathogens; 
# each row represents an antimicrobial resistance (AMR) profile defined by 
# susceptibility to three key antimicrobial classes (listed below organism names). 
# (A) Incidence of hospital-onset isolates per 1,000 admissions stratified by 
# organism and by AMR profile. 
# Points represent observed yearly incidence of hospital-onset isolates; 
# lines show linear regression fits with color-shaded 95% confidence intervals. 
# Grey backgrounds highlight AMR profiles where estimates indicated increasing trends. 
# ============================================================================ #
remove(list = ls())

# Load config file -------------------------------------------------------------
source(here::here("code", "00_config.R")) # phenotype_colors, bugs_ordered

# Load data --------------------------------------------------------------------
load(here::here("results/figure1", "mind_aim2-2_phenotype_inc_data.RData"))
load(here::here("results/figure1", "mind_aim2-2_phenotype_gee_time_trend_results.RData")) # df_time_trend

# ------------------------------------------------------------------------------
# Data prep
# ------------------------------------------------------------------------------
df_gee <- bind_rows(df_time_trend, .id = "organismofinterest") %>% 
  mutate(time_period = factor(time_period, levels = c("2020-2021", "2007-2019", "2007-2021")), 
         organismofinterest = factor(organismofinterest, levels = bugs_ordered)) %>% 
  mutate(
    trend_direction = case_when(
      lower.CL >= 0 | time.trend > 0 ~ "Increasing",
      upper.CL < 0 ~ "Decreasing",
      TRUE ~ "Stable")
  )

df_plot <- df_ab_year %>% 
  left_join(
    df_gee %>% 
      select(organismofinterest, antibiogram, trend_direction), 
    by = c("organismofinterest", "ab_name" = "antibiogram")
  ) %>% 
  mutate(org_ab = paste0(organismofinterest, "__", ab_name), 
         facet_lab = factor(facet_lab, levels = facet_labs))

# ------------------------------------------------------------------------------
# Pathogen-specific Phenotype color definitions
# ------------------------------------------------------------------------------
# Background data: one row per facet panel
facet_backgrounds <- df_plot %>%
  group_by(ab_name, facet_lab) %>%
  filter(n() > 0) %>%     # ensures panel has data
  distinct(trend_direction, .keep_all = TRUE)

panel_info <- df_plot %>%
  count(ab_name, facet_lab, name = "n_obs") %>%
  mutate(has_data = n_obs > 0)


# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
figure1A <- ggplot(df_plot, aes(x = ymd(date_year, truncated = 2), y = inc, color = org_ab)) +
  # 2) Faceting
  ggh4x::facet_nested(
    rows = vars(ab_name), 
    cols = vars(facet_lab),
    scales = "free_y", 
    independent = "y", 
    switch = "y"
  ) +
  # manual borders on data-containing panels only
  geom_rect(
    data = panel_info %>% filter(has_data),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    fill = "white",
    color = "black",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) + 
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
  # Tell ggplot: from now on, we’re starting a NEW fill scale
  ggnewscale::new_scale_fill() +
  # Smoothed trend line + ribbon, fill & color by org_ab
  stat_smooth(
    method  = "lm", 
    formula = y ~ x, 
    linewidth = 1.5,
    aes(fill = org_ab, color = org_ab),
    alpha = 0.3) +
  geom_point(size = 3.5, shape = 19) +
  
  scale_color_manual(values = phenotype_colors) +
  scale_fill_manual(values = phenotype_colors) +
  scale_y_continuous(limits = c(-0.2, NA), expand = c(0, 0), n.breaks = 4, position = "right") +
  coord_cartesian(ylim = c(0, NA)) +
  scale_x_date(date_label = "%Y", breaks = seq(as.Date("2007-01-01"), as.Date("2022-01-01"), by = "4 years")) +
  labs(x = "Year", y = "Incident isolates per 1000 admissions, n", color = "Phenotype", fill  = "Phenotype") +
  theme_template_time() +
  theme(
    legend.position   = "none", 
    strip.text        = element_text(face = "plain"),
    strip.text.x      = ggtext::element_markdown(size = 21), 
    strip.text.y      = element_text(face = "bold", size = 21), 
    axis.text.y       = element_text(size = 18),
    axis.title        = element_text(size = 20, face = "plain"),
    axis.text.x       = element_text(size = 20), 
    strip.background  = element_rect(fill = "lightgrey", color = "black", linewidth = 1.2),
    panel.spacing.y   = unit(0.6, "cm"),
    panel.spacing.x   = unit(0.3, "cm"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(), 
    panel.background   = element_blank(),
    panel.border = element_blank()
  )
print(figure1A)

# Save plot --------------------------------------------------------------------
saveRDS(figure1A, file = paste0(RESULTS, "/figure1/figure1A_phenotype_trends.rds"))

ggsave(filename = "figure1A_phenotype_trends.pdf",
       plot = figure1A,
       path = paste0(FIGURES, "/figure1/"),
       width = 21,
       height = 15)
