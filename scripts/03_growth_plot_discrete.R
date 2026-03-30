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

growth_label <- sprintf(
  "Growth rate (Oct 20-23)\nControl: %.1f um/day\nProbiotic: %.1f um/day",
  rate_control,
  rate_probiotic
)

p_growth_discrete <- ggplot() +
  geom_line(
    data = rep_means2,
    aes(x = day, y = mean_length, group = interaction(treatment, unit), color = treatment),
    linewidth = 1,
    alpha = 0.25
  ) +
  geom_point(
    data = rep_means2,
    aes(x = day, y = mean_length, group = interaction(treatment, unit), color = treatment),
    size = 2,
    alpha = 0.25
  ) +
  geom_line(
    data = emm_df2,
    aes(x = day, y = emmean, group = treatment, color = treatment),
    linewidth = 2.2
  ) +
  geom_point(
    data = emm_df2,
    aes(x = day, y = emmean, color = treatment),
    size = 4
  ) +
  geom_errorbar(
    data = emm_df2,
    aes(x = day, ymin = lower.CL, ymax = upper.CL, color = treatment),
    width = 0.10,
    linewidth = 1
  ) +
  scale_color_manual(
    values = c("control" = "darkgreen", "probiotic" = "purple4"),
    labels = c("control" = "Control", "probiotic" = "Probiotic"),
    na.translate = FALSE
  ) +
  scale_x_discrete(labels = day_labels) +
  labs(
    title = "Copepod growth trajectories across developmental sampling days",
    x = "Sampling day",
    y = "Length (um)",
    color = "Treatment"
  ) +
  theme_classic(base_size = 18, base_family = base_family) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18)
  ) +
  annotate(
    "label",
    x = "Oct23",
    y = 845,
    label = growth_label,
    hjust = 1,
    vjust = 1,
    size = 5.2,
    family = base_family,
    linewidth = 0.3,
    fill = "white"
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
