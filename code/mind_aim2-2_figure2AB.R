# ============================================================================ #
# Project: MInD Aim 2.2
# Title: Figure 2 for manuscript
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# ---------------------------------------------------------------------------- #
# 2A: Barplot of antibiotic use
# 2B: GEE time trend results for antibiotic use
# ============================================================================ #
remove(list = ls())
setwd("path/to/folder")

# Load packages
library(lubridate)
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)

# Load files
source("code/plotting_template.R")
load("data/mind_hospitallevel_abx_use.RData")

################################################################################
# Data preparation
# Load abx use data 
DT_abx_overall_year <- read.csv("data/mind_aim2-2_overall_abx_comb_sum_year")

# Table with abx use labels
abx_subset <- c("class_ANTISTAPHBETALACTAMS", "class_MACR", 
                "class_B_LAC", 
                "subclass_CPH34",
                "class_CPM",
                "class_FLQ")
abx_labels <- c("Anti-Staphylococcal Beta-lactams", "Macrolides", 
                "Beta-lactam/Beta-lactmase inhibitor", 
                "3rd/4th gen Cephalosporins",
                "Carbapenems", 
                "Fluoroquinolones")
abx_table <- data.frame(cbind(variable = abx_subset, label = abx_labels))

################################################################################
# Figure 1A: Plot of overall antibiotic use overall
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

# ---------------------------------------------------------------------------- #
# PLOT
(figure2A <- ggplot(data, 
                    aes(x = ymd(date_year, truncated=2), 
                        y = rate, group = date_year)) +
    facet_wrap(.~label, scales = "free_y", nrow = 7) + 
    geom_rect(aes(xmin = xmin1, xmax = xmax1, 
                  ymin=-Inf, ymax=Inf), fill = "lightblue", alpha = 0.5) + 
    geom_rect(aes(xmin = xmin2, xmax = xmax2, 
                  ymin=-Inf, ymax=Inf), fill = "midnightblue", alpha = 0.5) + 
    geom_bar(stat = "identity", color="black", fill = "darkgrey", linewidth=2) +
    labs(x = "Calendar year", 
         y = "Days of therapy per 1,000 inpatient days") + 
    scale_x_date(date_labels = "%Y", breaks = seq(as.Date("2007-01-01", format = "%Y-%m-%d"), 
                                                  as.Date("2022-12-31", format = "%Y-%m-%d"), by = "2 years")) + 
    scale_y_continuous(expand= expansion(mult=c(0.1,0.2))) + 
    theme_template_white() + 
    theme(axis.text.x = element_text(size = 13, angle = 0, hjust=0.5, vjust=0.5), 
          strip.background = element_rect(fill="white", color = "black", linewidth = 1.2),
          axis.title = element_text(size = 16)))
ggsave(figure2A, file = "figures/abx_use/mind_aim2-2_figure2A_v4.pdf", 
       width = 5, height = 14)

################################################################################
# Figure 2B: Plot of GEE time trend results
# Load GEE time trend results
df_time_trend <- read.csv("results/mind_aim2-2_abx_use_gee_time_trend_results.csv")

# Data preparation
df_gee <- df_time_trend %>% 
  filter(antibiotic%in%abx_subset) %>% 
  left_join(abx_table, by = c("antibiotic" = "variable"))%>% 
  mutate(label = factor(label, levels = abx_labels), 
         time_period = factor(time_period, 
                              levels = c("2020-2022", "2007-2019", "2007-2022", "2012-2022", "2007-2011", "2016-2022", "2007-2015")))
df_gee$period <- c(rep(c("before intervention/event", "after intervention/event"), 4), 
                   c("before intervention/event", "after intervention/event"), "overall")
df_gee$period <- factor(df_gee$period, 
                        levels = c("before intervention/event", "after intervention/event", "overall"))

# ---------------------------------------------------------------------------- #
# PLOT
(figure2B <- ggplot(df_gee, 
                    aes(x=time.trend, y = time_period)) + 
    facet_wrap(.~label, ncol = 1, scales = "free_y") + 
    geom_vline(xintercept = 0, linetype = 2, linewidth = 1.3, color = "darkred") + 
    geom_errorbar(aes(xmin=lower.CL, xmax=upper.CL, color = period), width = 0., linewidth = 1.8) + 
    # geom_point(aes(fill = "white", color = "white"), stroke = 2, size = 5) + 
    geom_point(aes(color = period), fill = "white", stroke = 2, size = 3, shape =21) + 
    labs(x = "Average annual percentage change (%)", color = "Time period", fill = "Time period") + 
    scale_x_continuous(limits = c(min(df_gee$lower.CL), max(df_gee$upper.CL)), breaks = seq(min(df_gee$lower.CL),  max(df_gee$upper.CL), by = 5))+
    # scale_fill_manual(values = rev(c("black", "midnightblue", "deepskyblue3"))) +
    scale_color_manual(values = rev(c("black", "midnightblue", "deepskyblue3"))) +
    guides(color = guide_legend(nrow=3, byrow=T)) + 
    theme_template_white() + 
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_text(size=16),
          strip.background = element_rect(fill="white", color = "black", linewidth=1.2), 
          legend.position = "right", 
          legend.title = element_text(size=18, hjust = 0.5)))


################################################################################
# Figure 1: Combine Figure 1A and B
combined_figure <- figure2A + plot_spacer() +  figure2B + plot_layout(guides = "collect", widths = c(1,0.1,1)) & theme(legend.position = "bottom")
(figure2 <- combined_figure + plot_annotation(tag_level = 'A'))

ggsave(figure2, file = "figures/abx_use/mind_aim2-2_figure2AB.pdf", 
       width = 12, height = 14)
