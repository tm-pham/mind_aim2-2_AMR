# ============================================================================ #
# Project: MInD Aim 2.2
# PLOT
# Author: Thi Mui Pham, t.m.pham@hsph.harvard.edu
# Title: Figure 3 - multinomial logit model forest plot of fluoroquinolone 
# effects on antibiotic resistance outcomes across four pathogens
# ============================================================================ #
################################################################################
# Load config file
source(here::here("00_config.R")) # phenotype_colors, bugs_ordered

################################################################################
# Load data 
sa_data <- read.csv(paste0(RESULTS, "/figure3/mind_aim2-2_SA_4_mblogit_table_abx_results.csv"))
ec_data <- read.csv(paste0(RESULTS, "/figure3/mind_aim2-2_EC_6_mblogit_table_abx_results.csv"))
kp_data <- read.csv(paste0(RESULTS, "/figure3/mind_aim2-2_KP_5_mblogit_table_abx_results.csv"))
pa_data <- read.csv(paste0(RESULTS, "/figure3/mind_aim2-2_PA_5_mblogit_table_abx_results.csv"))

################################################################################
# Data preparation
################################################################################
# Add pathogen identifier
sa_data$pathogen <- "Staphylococcus aureus"
ec_data$pathogen <- "Escherichia coli"
kp_data$pathogen <- "Klebsiella pneumoniae"
pa_data$pathogen <- "Pseudomonas aeruginosa"

# Combine all datasets
all_data <- bind_rows(sa_data, ec_data, kp_data, pa_data)

# Filter for FQL term where flq == 1
# Filter for FQL term
fql_data <- all_data %>%
  filter(term == "FQL") %>%
  select(pathogen, ab_name, term, flq, estimate_interp, 
         conf10_interp, conf90_interp,
         conf2.5_interp, conf97.5_interp,
         p.value) %>%
  mutate(
    pathogen = factor(pathogen, 
                      levels = bugs_ordered),
    # Add antibiogram label
    label = paste0(term, "\n", ab_name),
    # Check if susceptible to FQL (starts with "S")
    fql_susceptible = (flq == 1)
  ) %>%
  arrange(pathogen, estimate_interp) %>%
  group_by(pathogen) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>% 
  mutate(org_ab = factor(paste0(pathogen, "__", ab_name), levels = org_ab_order), 
         ab_name = factor(ab_name, levels = rev(FQL_phenotype_order))) %>% 
  left_join(facet_table, by = c("pathogen" = "organismofinterest")) %>% 
  group_by(pathogen) %>% 
  arrange(estimate_interp) %>%
  mutate(y_pos = row_number()) %>% 
  ungroup()


# Create shade data - keep the SAME structure as fql_data
shade_data <- fql_data %>%
  filter(substr(as.character(ab_name), 1, 1) == "R")

################################################################################
# PLOT
################################################################################
# Create plot
(figure3 <- ggplot(fql_data, aes(y = ab_name, x = estimate_interp)) +
  # Grey background for R-* phenotypes
  geom_rect(
    data = shade_data,
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
               linewidth = 0.7, color = "gray30") +
  # 80% CI - thick line
  geom_segment(aes(x = conf10_interp, xend = conf90_interp,
                   y = ab_name,
                   yend = ab_name, color = org_ab),
               linewidth = 2.8) +
  scale_color_manual(values = phenotype_colors) +
  scale_x_continuous(limits =c(min(fql_data$conf2.5_interp)-0.2, max(fql_data$conf97.5_interp)+0.2), 
                     n.breaks = 10) + 
  # Point estimate
  geom_point(size = 3.5, shape = 21, fill = "white", color = "black") +
  # Facet by pathogen (1 column)
  facet_wrap(~ fct_reorder(facet_labs, pathogen), ncol = 1, scales = "free_y") +
  
  # Labels
  labs(x = "Effect of FQL prescribing\nIncrease in odds (%, compared to S-S-S)",
       y = "AMR profiles (outcome)") +
  
  # Clean theme matching reference
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "plain", size = 16, margin = margin(t = 5)),
    axis.title.y = element_text(face = "plain", size = 16, margin = margin(r = 12)),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.ticks.y = element_blank(),
    strip.text.x      = ggtext::element_markdown(size = 16), 
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(1, "lines"),
    plot.margin = margin(10, 30, 10, 10),
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(), 
    panel.background = element_rect(color = "white")
  ))

# Save the plot at high resolution for poster
ggsave(plot = figure3, paste0(FIGURES, "/figure3/figure3.pdf"), 
       width = 5, height = 11)

