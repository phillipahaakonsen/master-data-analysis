# PCoA D21 plot from PAST scores export
# Input: PCoA_D21_scores.csv and PCoA_D21_summary.csv in project root
# Output: pcoa_d21_plot.png in results directory

source("scripts/08_pcoa_plot_helper.R", local = TRUE)
make_pcoa_plot(
  "PCoA_D21_scores.csv",
  "PCoA_D21_summary.csv",
  "D21",
  "pcoa_d21_plot.png",
  probiotic_color = "#4169E1",
  control_color = "#FF9500",
  x_limits = c(-0.300, 0.375),
  y_limits = c(-0.20, 0.25),
  x_breaks = seq(-0.300, 0.375, by = 0.075),
  y_breaks = seq(-0.20, 0.25, by = 0.05),
  show_legend = FALSE
)
