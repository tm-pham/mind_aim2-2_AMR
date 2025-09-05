# ============================================================================ #
# Project: MInD Aim 2.2
# Title: Figure 1 for manuscript
# Author: Thi Mui Pham, mui.k.pham@gmail.com
# ---------------------------------------------------------------------------- #
# A: Barplot of phenotype incidence 
# B: GEE time trend results for phenotype incidence
# ============================================================================ #
remove(list = ls())
setwd("path/to/folder")

# Load packages
library(lubridate)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load files
source("code/plotting_template.R")

bugs_ordered <- c("Staphylococcus aureus", 
                  "Escherichia coli", 
                  "Klebsiella pneumoniae",
                  "Pseudomonas aeruginosa")

################################################################################
# Figure A: Antibiogram incidence plot
load("P:/ORD_Samore_202109019D/Mui/results/aim2-2/mind_aim2-2_figureA_phenotype_trend.RData") # figureA


################################################################################
# Figure B: Plot of GEE time trend results
df_time_trend <- read.csv("results/mind_aim2-2_phenotype_gee_time_trend_results.csv", header = T)

df_gee <- bind_rows(df_time_trend, .id = "organismofinterest") %>% 
  mutate(time_period = factor(time_period, levels = c("2020-2021", "2007-2019", "2007-2021")), 
         organismofinterest = factor(organismofinterest, levels = bugs_ordered))

(figureB <- ggplot(df_gee, 
                   aes(x=as.numeric(time.trend), y = 1)) + 
    ggh4x::facet_nested(rows = vars(antibiogram), 
                        cols = vars(organismofinterest), 
                        scales = "free", render_empty = FALSE) + 
    geom_vline(xintercept = 0, linetype = 2, linewidth = 1.3, color = "darkred") + 
    geom_errorbar(aes(xmin=as.numeric(lower.CL), xmax=as.numeric(upper.CL), color = time_period), width = 0., linewidth = 1.5) + 
    geom_point(aes(color = time_period), fill = "white", stroke = 2, size = 5, shape =21) + 
    labs(x = "Average annual percentage change (%)") + 
    scale_color_manual(values = rev(c("black"))) +
    theme_template() + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=20),
          axis.text.x = element_text(size=18),
          legend.position = "none", 
          plot.tag = element_text(size=20),
          strip.text.x = element_text(face = "bold.italic"), 
          strip.text.y = element_text(face = "bold"),
          panel.spacing.x = unit(1, "cm"),
          strip.background = element_rect(fill = "lightgrey")))


################################################################################
# Figure: Combine Figure A and B
(figure <- (figureA/figureB) +  
   plot_layout(nrow = 2, heights = c(1.4, 1)) + 
   plot_annotation(tag_level = 'A') +
   theme(plot.tag = element_text(size=18)))

ggsave(figure, file = "figures/mind_aim2-2_figure1AB.pdf", 
       width = 16, height = 26)
