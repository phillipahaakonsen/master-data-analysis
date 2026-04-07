# Growth trajectory plot (discrete day axis)
# Prerequisite:
# source("scripts/02_growth_helpers.R")

library(dplyr)
library(ggplot2)

if (!exists("rep_means2") || !exists("emm_df2")) {
  stop("Run scripts/02_growth_helpers.R first")
}

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

base_family <- "serif"
base_size <- 18
axis_title_size <- 22
axis_text_size <- 18
legend_title_size <- 22
legend_text_size <- 18
probiotic_rep_color <- "#6F8FF0"
control_rep_color <- "#F4B15A"
probiotic_mean_color <- "#4367E0"
control_mean_color <- "#F09A2E"

rep_means2_plot <- rep_means2 %>%
  mutate(color_key = if_else(treatment == "probiotic", "probiotic_rep", "control_rep"))

emm_df2_plot <- emm_df2 %>%
  mutate(color_key = if_else(treatment == "probiotic", "probiotic_mean", "control_mean"))

p_growth_discrete <- ggplot() +
  geom_line(
    data = rep_means2_plot,
    aes(x = day, y = mean_length, group = interaction(treatment, unit), color = color_key),
    linewidth = 1,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_point(
    data = rep_means2_plot,
    aes(x = day, y = mean_length, group = interaction(treatment, unit), color = color_key),
    size = 2,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_line(
    data = emm_df2_plot,
    aes(x = day, y = emmean, group = treatment, color = color_key),
    linewidth = 2.2
  ) +
  geom_point(
    data = emm_df2_plot,
    aes(x = day, y = emmean, color = color_key),
    size = 4
  ) +
  geom_errorbar(
    data = emm_df2_plot,
    aes(x = day, ymin = lower.CL, ymax = upper.CL, color = color_key),
    width = 0.10,
    linewidth = 1
  ) +
  scale_color_manual(
    values = c(
      "control_rep" = control_rep_color,
      "control_mean" = control_mean_color,
      "probiotic_rep" = probiotic_rep_color,
      "probiotic_mean" = probiotic_mean_color
    ),
    breaks = c("probiotic_mean", "control_mean"),
    labels = c("probiotic_mean" = "Probiotic-treated", "control_mean" = "Control"),
    na.translate = FALSE
  ) +
  scale_x_discrete(labels = day_labels) +
  labs(
    x = "Sampling day",
    y = "Length (um)",
    color = "Treatment"
  ) +
  theme_classic(base_size = base_size, base_family = base_family) +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    legend.text = element_text(size = legend_text_size),
    legend.position = "right"
  ) +
  coord_cartesian(clip = "off")

print(p_growth_discrete)

growth_discrete_file <- file.path(results_dir, "growth_trajectory_discrete.png")
ggsave(
  filename = growth_discrete_file,
  plot = p_growth_discrete,
  width = 12,
  height = 8,
  dpi = 300
)
cat("Saved plot:", growth_discrete_file, "\n")
