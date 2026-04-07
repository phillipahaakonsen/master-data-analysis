# PCoA D16 plot from PAST scores export
# Input: PCoA_D16_scores.csv and PCoA_D16_summary.csv in project root
# Output: pcoa_d16_plot.png in results directory

source("scripts/08_pcoa_plot_helper.R", local = TRUE)
make_pcoa_plot(
  "PCoA_D16_scores.csv",
  "PCoA_D16_summary.csv",
  "D16",
  "pcoa_d16_plot.png",
  show_legend = FALSE
)
