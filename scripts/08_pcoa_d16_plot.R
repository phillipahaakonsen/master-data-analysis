# PCoA D16 plot from PAST scores export
# Input: PCoA_D16_scores.csv and PCoA_D16_summary.csv in project root
# Output: pcoa_d16_plot.png in results directory

source("scripts/08_pcoa_plot_helper.R", local = TRUE)
data_dir <- getOption("project_data_dir", getwd())
make_pcoa_plot(
  file.path(data_dir, "PCoA_D16_scores.csv"),
  file.path(data_dir, "PCoA_D16_summary.csv"),
  "D16",
  "pcoa_d16_plot.png",
  show_legend = FALSE
)
