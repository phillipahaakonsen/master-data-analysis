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
base_family <- "serif"
base_size <- 18
axis_title_size <- 22
axis_text_size <- 18

sst_plot <- sst_main %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Probiotic", "Control"),
      labels = c("Probiotic-treated", "Control")
    )
  )

sum_plot <- sst_plot %>%
  group_by(Treatment) %>%
  summarise(
    mean_survival = mean(pct_alive, na.rm = TRUE),
    se_survival = sd(pct_alive, na.rm = TRUE) / sqrt(sum(!is.na(pct_alive))),
    .groups = "drop"
  )

p_survival <- ggplot(sst_plot, aes(x = Treatment, y = pct_alive)) +
  geom_jitter(
    aes(fill = Treatment),
    width = 0.08,
    size = 3.2,
    shape = 21,
    stroke = 0.5,
    color = "grey25",
    alpha = 0.95
  ) +
  geom_errorbar(
    data = sum_plot,
    aes(
      x = Treatment,
      ymin = mean_survival - se_survival,
      ymax = mean_survival + se_survival
    ),
    width = 0.06,
    linewidth = 0.6,
    color = "grey30",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = sum_plot,
    aes(x = Treatment, y = mean_survival, fill = Treatment),
    shape = 23,
    size = 4.4,
    stroke = 0.7,
    color = "grey20",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = sum_plot,
    aes(x = Treatment, y = mean_survival),
    shape = 4,
    stroke = 0.8,
    size = 3.4,
    color = "grey30",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("Probiotic-treated" = "#4367E0", "Control" = "#F09A2E")
  ) +
  labs(
    x = NULL,
    y = "Survival (%)"
  ) +
  theme_classic(base_size = base_size, base_family = base_family) +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size, color = "black"),
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
