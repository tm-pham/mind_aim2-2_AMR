# ============================================================================ #
# Project: MInD Aim 2.2
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# Title: Figure 2 for manuscript
# ---------------------------------------------------------------------------- #
# Figure 2. Trends in antimicrobial prescribing in the Veterans Affairs 
# Healthcare Administration, February 1, 2007–March 31, 2022.
# (A) Bar plots show overall antimicrobial prescribing rates, expressed as days 
# of therapy per 1,000 patient-days, in 138 VA medical centres for key 
# antimicrobial classes. Blue-shaded regions indicate time periods corresponding 
# to those analysed in panel B.
# (B) GEE time trend results for antibiotic use
# Average annual percentage change (AAPC) estimates from time-trend analyses 
# using generalized estimating equations for the antimicrobial classes in panel 
# A. Positive AAPC values indicate increasing trends; negative values indicate 
# decreasing trends. For anti-staphylococcal β-lactams, macrolides, 
# third-generation cephalosporins, and fluoroquinolones, 
# interactions with the COVID-19 period were evaluated, and separate trend 
# estimates are presented if the interaction term was statistically significant 
# (overall p < 0.05, Bonferroni-corrected). 
# For fluoroquinolones, the evaluated interaction term was not statistically 
# significant and thus time trends were reported for the whole study period. 
# For β-lactam/β-lactamase inhibitors and carbapenems, interactions with 2015 
# drug shortages and 2011 antimicrobial stewardship initiatives, respectively, 
# were similarly evaluated. “Before” and “after” denote the periods preceding 
# and following these breakpoints, as shown on the y-axis. 
# Further details are provided in the Appendix (pp. 36–37).
# ============================================================================ #
remove(list = ls())
# Load config file -------------------------------------------------------------
source(here::here("code", "00_config.R")) # phenotype_colors, bugs_ordered

# ------------------------------------------------------------------------------
# Data preparation
# Load antibiotic use data 
load(here::here("results/figure2/", "mind_aim2-2_overall_abx_comb_CPH34_sum_year.RData")) # DT_abx_overall_year

# Create table with abx use variable names and corresponding labels for plotting
abx_subset <- c("class_ANTISTAPHBETALACTAMS", 
                "class_MACR", 
                "class_B_LAC", 
                "subclass_CPH34",
                "class_CPM",
                "class_FLQ")
abx_labels <- c("Anti-staphylococcal beta-lactams", 
                "Macrolides", 
                "Beta-lactam/Beta-lactmase inhibitor", 
                "3rd/4th gen cephalosporins",
                "Carbapenems", 
                "Fluoroquinolones")
abx_table <- data.frame(cbind(variable = abx_subset, label = abx_labels))

# ------------------------------------------------------------------------------
# Figure 1A: Plot of overall antibiotic use overall
# ------------------------------------------------------------------------------
# Define time periods for background shading 
df_shade <- as.data.frame(cbind(label = abx_labels))
df_shade$xmin1 <- c(as.Date("2006-07-01", format = "%Y-%m-%d"), 
                    as.Date("2006-07-01", format = "%Y-%m-%d"), 
                    as.Date("2006-07-01", format = "%Y-%m-%d"), 
                    as.Date("2006-07-01", format = "%Y-%m-%d"), 
                    as.Date("2006-07-01", format = "%Y-%m-%d"), 
                    NA)
df_shade$xmax1 <- c(as.Date("2019-06-30", format = "%Y-%m-%d"),
                    as.Date("2019-06-30", format = "%Y-%m-%d"), 
                    as.Date("2015-06-30", format = "%Y-%m-%d"),
                    as.Date("2019-06-30", format = "%Y-%m-%d"),
                    as.Date("2011-06-30", format = "%Y-%m-%d"),
                    NA)
df_shade$xmin2 <- c(as.Date("2019-06-30", format = "%Y-%m-%d"), 
                    as.Date("2019-06-30", format = "%Y-%m-%d"), 
                    as.Date("2015-06-30", format = "%Y-%m-%d"), 
                    as.Date("2019-06-30", format = "%Y-%m-%d"),
                    as.Date("2011-06-30", format = "%Y-%m-%d"), 
                    NA)
df_shade$xmax2 <- rep(as.Date("2022-08-01", format = "%Y-%m-%d"), 6)

data <- DT_abx_overall_year %>% 
  filter(variable%in%abx_subset) %>%
  left_join(abx_table) %>% 
  as.data.frame() %>% 
  left_join(df_shade, by = "label") %>% 
  mutate(label = factor(label, levels = abx_labels))

# PLOT -------------------------------------------------------------------------
(figure2A <- ggplot(data, 
                    aes(x = ymd(date_year, truncated=2), 
                        y = rate, group = date_year)) +
    facet_wrap(.~label, scales = "free_y", nrow = 7) + 
    geom_rect(aes(xmin = xmin1, xmax = xmax1, 
                  ymin=-Inf, ymax=Inf), fill = "lightblue", alpha = 0.5) + 
    geom_rect(aes(xmin = xmin2, xmax = xmax2, 
                  ymin=-Inf, ymax=Inf), fill = "midnightblue", alpha = 0.5) + 
    geom_bar(stat = "identity", color="black", fill = "darkgrey", linewidth=2) +
    labs(x = "Year", 
         y = "Days of therapy per 1,000 inpatient days") + 
    scale_x_date(date_labels = "%Y", breaks = seq(as.Date("2007-01-01", format = "%Y-%m-%d"), 
                                                  as.Date("2022-12-31", format = "%Y-%m-%d"), by = "2 years")) + 
    scale_y_continuous(expand= expansion(mult=c(0.1,0.2))) + 
    theme_template() + 
    theme(axis.text.x = element_text(size = 13, angle = 0, hjust=0.5, vjust=0.5), 
          axis.title = element_text(size = 16, face = "plain"),
          strip.background = element_rect(fill="white", color = "black", linewidth = 1.2),
          strip.text = element_text(size=17, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA), 
          panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank()))

ggsave(figure2A, file = paste0(FIGURES, "/figure2/figure2A.pdf"), 
       width = 5, height = 14)

# ------------------------------------------------------------------------------
# Figure 2B: Plot of GEE time trend results
# ------------------------------------------------------------------------------
# Load GEE time trend results
load(paste0(RESULTS, "/figure2/mind_aim2-2_abx_use_gee_time_trend_results.RData")) # df_time_trend

# ------------------------------------------------------------------------------
# Data preparation 
df_gee <- df_time_trend %>% 
  filter(antibiotic%in%abx_subset) %>% 
  left_join(abx_table, by = c("antibiotic" = "variable"))%>% 
  mutate(label = factor(label, levels = abx_labels), 
         time_period = factor(time_period, 
                              levels = c("2020-2022", 
                                         "2007-2019", 
                                         "2007-2022", 
                                         "2012-2022", 
                                         "2007-2011", 
                                         "2016-2022", 
                                         "2007-2015"),
                              labels = c("2020-2022\n(COVID-19 period)", 
                                         "2007-2019\n(pre-COVID-19 period)", 
                                         "2007-2022", 
                                         "2012-2022\n(post-stewardship initiative)", 
                                         "2007-2011\n(pre-stewardship initiative)", 
                                         "2016-2022\n(post-pip/tazo shortage)", 
                                         "2007-2015\n(pre-pip/tazo shortage)")))
df_gee$period <- c(rep(c("before intervention/event", "after intervention/event"), 4), 
                   c("before intervention/event", "after intervention/event"), "overall")
df_gee$period <- factor(df_gee$period, 
                        levels = c("before intervention/event", "after intervention/event", "overall"))

# PLOT -------------------------------------------------------------------------
(figure2B <- ggplot(df_gee, 
                    aes(x=time.trend, y = time_period)) + 
    facet_wrap(.~label, ncol = 1, scales = "free_y") + 
    geom_vline(xintercept = 0, linetype = 1, linewidth = 0.5, color = "black") + 
    geom_errorbar(aes(xmin=lower.CL, xmax=upper.CL, color = period), width = 0., linewidth = 1.4) + 
    # geom_point(aes(fill = "white", color = "white"), stroke = 2, size = 5) + 
    geom_point(aes(color = period), fill = "white", stroke = 1.4, size = 3.5, shape =21) + 
    labs(x = "Average annual percentage change (%)", color = "Time period", fill = "Time period") + 
    scale_x_continuous(limits = c(min(df_gee$lower.CL), max(df_gee$upper.CL)), breaks = seq(min(df_gee$lower.CL),  max(df_gee$upper.CL), by = 5))+
    # scale_fill_manual(values = rev(c("black", "midnightblue", "deepskyblue3"))) +
    scale_color_manual(values = rev(c("black", "midnightblue", "deepskyblue3"))) +
    guides(color = guide_legend(nrow=3, byrow=T)) + 
    theme_template() + 
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_text(size=16, face = "plain"),
          strip.background = element_rect(fill="white", color = "black", linewidth=1.2), 
          strip.text = element_text(size=17, face = "bold"),
          legend.position = "right", 
          legend.title = element_text(size=17, hjust = 0.5), 
          plot.background = element_rect(fill = "transparent", color = NA), 
          panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank()))


# ------------------------------------------------------------------------------
# Figure 2: Combine Figure 2A and B
combined_figure <- figure2A + plot_spacer() +  figure2B + plot_layout(guides = "collect", widths = c(1,0.1,1)) & theme(legend.position = "bottom")
(figure2 <- combined_figure + plot_annotation(tag_level = 'A') & theme(plot.tag = element_text(size = 18, face = "plain")))

# Save figure ------------------------------------------------------------------
ggsave(figure2, file = paste0(FIGURES, "/figure2/figure2AB.pdf"), 
       width = 13.5, height = 14.5)
