# Growth trajectory plot (numeric x-axis + right-side labels)
# Prerequisite:
# source("scripts/02_growth_helpers.R")

library(dplyr)
library(ggplot2)
library(tibble)

if (!exists("rep_means2") || !exists("emm_df2")) {
  stop("Run scripts/02_growth_helpers.R first")
}

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

label_df <- tibble(
  treatment = c("probiotic", "control"),
  x = c(22.6, 22.6),
  y = c(830, 780),
  lab = c(
    sprintf("Rate (Oct. 20th-23rd): %.1f um/day", rate_probiotic),
    sprintf("Rate (Oct. 20th-23rd): %.1f um/day", rate_control)
  )
)

p_growth_numeric <- ggplot() +
  geom_line(
    data = rep_means2,
    aes(x = x, y = mean_length, group = interaction(treatment, unit), color = treatment),
    linewidth = 1,
    alpha = 0.25
  ) +
  geom_point(
    data = rep_means2,
    aes(x = x, y = mean_length, group = interaction(treatment, unit), color = treatment),
    size = 2,
    alpha = 0.25
  ) +
  geom_line(
    data = emm_df2,
    aes(x = x, y = emmean, group = treatment, color = treatment),
    linewidth = 2.2
  ) +
  geom_point(
    data = emm_df2,
    aes(x = x, y = emmean, color = treatment),
    size = 4
  ) +
  geom_errorbar(
    data = emm_df2,
    aes(x = x, ymin = lower.CL, ymax = upper.CL, color = treatment),
    width = 0.6,
    linewidth = 1
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = lab, color = treatment),
    hjust = 0,
    size = 5.5,
    family = "serif"
  ) +
  scale_color_manual(
    values = c("control" = "darkgreen", "probiotic" = "purple4"),
    labels = c("control" = "Control", "probiotic" = "Probiotic"),
    breaks = c("control", "probiotic"),
    na.translate = FALSE
  ) +
  scale_x_continuous(
    breaks = c(15, 20, 23),
    labels = c("Oct. 15th", "Oct. 20th", "Oct. 23rd"),
    expand = expansion(mult = c(0.08, 0.35))
  ) +
  coord_cartesian(ylim = c(60, 860), clip = "off") +
  labs(
    title = "Copepod growth trajectories across developmental sampling days",
    x = "Sampling day",
    y = "Length (um)",
    color = "Treatment"
  ) +
  theme_classic(base_size = 18, base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    plot.margin = margin(10, 40, 10, 10)
  )

print(p_growth_numeric)

growth_numeric_file <- file.path(results_dir, "growth_trajectory_numeric.png")
ggsave(
  filename = growth_numeric_file,
  plot = p_growth_numeric,
  width = 12,
  height = 8,
  dpi = 300
)
cat("Saved plot:", growth_numeric_file, "\n")
