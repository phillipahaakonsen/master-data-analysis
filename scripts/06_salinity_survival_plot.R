# Salinity survival plot
# Prerequisite:
# source("scripts/05_salinity_stress_stats.R")

library(dplyr)
library(ggplot2)

if (!exists("sst_main")) {
  stop("Run scripts/05_salinity_stress_stats.R first")
}

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

sst_plot <- sst_main %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Control", "Probiotic"),
      labels = c("Control", "Probiotic-treated")
    )
  )

sum_plot <- sst_plot %>%
  group_by(Treatment) %>%
  summarise(
    mean_survival = mean(pct_alive, na.rm = TRUE),
    sd_survival = sd(pct_alive, na.rm = TRUE),
    .groups = "drop"
  )

p_survival <- ggplot(sst_plot, aes(x = Treatment, y = pct_alive, color = Treatment)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    linewidth = 0.8,
    alpha = 0.18
  ) +
  geom_jitter(width = 0.06, size = 2.8, alpha = 0.9) +
  geom_errorbar(
    data = sum_plot,
    aes(
      x = Treatment,
      ymin = mean_survival - sd_survival,
      ymax = mean_survival + sd_survival,
      color = Treatment
    ),
    width = 0.08,
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = sum_plot,
    aes(x = Treatment, y = mean_survival),
    shape = 18,
    size = 4,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c("Control" = "darkgreen", "Probiotic-treated" = "purple4")
  ) +
  labs(
    x = "Treatment group",
    y = "Survival (%)"
  ) +
  theme_classic(base_family = "serif") +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15, color = "black"),
    legend.position = "none"
  )

print(p_survival)

survival_plot_file <- file.path(results_dir, "salinity_survival_plot.png")
ggsave(
  filename = survival_plot_file,
  plot = p_survival,
  width = 8,
  height = 6,
  dpi = 300
)
cat("Saved plot:", survival_plot_file, "\n")
