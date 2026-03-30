# Script index for this project
#
# One-command run:
# source("run_all.R")
#
# Run these scripts by task:
# 1) D16 similarity boxplot
#    source("scripts/01_d16_similarity_boxplot.R")
#
# 2) Copepod growth trajectories (requires rep_means and emm_df in memory)
#    source("scripts/07_growth_inputs_from_length_measurements.R")
#    source("scripts/02_growth_helpers.R")
#    source("scripts/03_growth_plot_discrete.R")
#    source("scripts/04_growth_plot_numeric.R")
#
# 3) Salinity stress analysis
#    source("scripts/05_salinity_stress_stats.R")
#    source("scripts/06_salinity_survival_plot.R")
#
# Notes:
# - 03 and 04 produce alternative visual styles of the same growth trajectory.
# - 06 expects that 05 has already created sst_main and summary tables.
