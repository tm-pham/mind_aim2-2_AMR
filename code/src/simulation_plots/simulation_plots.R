# ============================================================================ #
# Project: MInD Aim 2.2
# RUNFILE
# Title: Plot bias and coverage of the siimulation analyses 
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# ============================================================================ #
orderly2::orderly_dependency(
  "simulation_mblogit_fit", 
  "latest", 
  files = c("bias_coverage_df.RData"))

load("bias_coverage_df.RData")

library(patchwork)

antibiogram_colors <- c("#E91E63", "#E64A19", "#FFC107", "#1E88E5", "#388E3C")

################################################################################
# Bias
# Average of estimates minus the true value 
(bias_plot <- ggplot(bias_df, aes(x = antibiogram, y = bias, color = antibiogram, fill = antibiogram)) +
   facet_wrap(.~term, scales = "fixed") + 
   geom_boxplot(, color = "black") + 
   geom_jitter(width = 0.1, shape = 21, color = "white") +
   geom_hline(yintercept = 0, linetype = "dashed") +
   scale_color_manual(values = antibiogram_colors) +
   scale_fill_manual(values = antibiogram_colors) +
   labs(x = "Antibiogram", 
        y = "Bias") +
   theme_template() + 
   theme(axis.title.x = element_blank(), 
         axis.title.y = element_text(face = "plain")))

################################################################################
# Coverage 
# Ratio of times the confidence intervals overlaps the true mean 
(coverage_plot <- ggplot(coverage_df, aes(x = antibiogram, y = coverage, fill = antibiogram)) +
   facet_wrap(.~term, scales = "fixed") +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_hline(yintercept = 0.95, linetype = "dashed") +
   scale_y_continuous(labels = scales::percent_format(scale = 100)) +
   scale_fill_manual(values = antibiogram_colors) +
   labs(y = "Coverage") +
   theme_template() + 
   theme(axis.title.x = element_blank(), 
         axis.title.y = element_text(face = "plain")))

################################################################################
# Confidence interval width and envelope 
# Define a scaling factor for avg_ci_width so it can be overlaid with CI bounds.
# Adjust scale_factor so that the bar heights align visually with the envelope.
scale_factor <- 1  # Modify as needed

(ci_env_plot <- ggplot(ci_env_df, aes(x = antibiogram)) +
    facet_wrap(~term, scales = "fixed") +
    # Plot the envelope as error bars and points for each antibiogram
    geom_errorbar(aes(ymin = overall_lower, ymax = overall_upper, color = antibiogram),
                  width = 0.3, position = position_dodge(width = 0.6)) +
    geom_point(aes(y = overall_lower, color = antibiogram),
               position = position_dodge(width = 0.8)) +
    geom_point(aes(y = overall_upper, color = antibiogram),
               position = position_dodge(width = 0.8)) +
    # Plot the average CI width as bars, scaled appropriately
    geom_bar(aes(y = avg_ci_width * scale_factor, fill = antibiogram),
             stat = "identity", position = position_dodge(width = 0.5), alpha = 0.6, color = "black") +
    scale_color_manual(values = antibiogram_colors) +
    scale_fill_manual(values = antibiogram_colors) +
    labs(x = "Antibiogram") +
    # Set up the primary y-axis for the CI bounds and a secondary y-axis for avg_ci_width
    scale_y_continuous(
      name = "Envelope",
      sec.axis = sec_axis(~ . / scale_factor, name = "Average Confidence Interval Width")) +
    theme_template() + 
    theme(axis.title.y = element_text(face = "plain"))
)


# Combine the plots and add annotations
figure <- (bias_plot/coverage_plot/ci_env_plot) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(figure, 
       file = "simulation_analysis_plots.pdf", 
       width = 18, height = 15)


