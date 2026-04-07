# Run all project analyses in order and save artifacts to results/

results_dir <- file.path(getwd(), "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
options(project_results_dir = results_dir)

scripts_to_run <- c(
  "scripts/01_d16_similarity_boxplot.R",
  "scripts/01b_d21_similarity_boxplot.R",
  "scripts/01c_d24_similarity_boxplot.R",
  "scripts/01d_all_days_similarity_boxplot.R",
  "scripts/01e_copepod_stage_similarity_boxplots.R",
  "scripts/02_growth_helpers.R",
  "scripts/03_growth_plot_discrete.R",
  "scripts/04_growth_plot_numeric.R",
  "scripts/05_salinity_stress_stats.R",
  "scripts/06_salinity_survival_plot.R",
  "scripts/08_pcoa_d16_plot.R",
  "scripts/08b_pcoa_d21_plot.R",
  "scripts/08c_pcoa_d24_plot.R",
  "scripts/09_pcoa_all_plot.R",
  "scripts/09b_pcoa_all_copepods_plot.R",
  "scripts/10_alpha_diversity_boxplots.R",
  "scripts/10b_alpha_diversity_stage_boxplots.R",
  "scripts/11_asv_relative_abundance_scatterplots.R",
  "scripts/16_asv_cop_ctr_vs_cop_pro_scatterplots.R"
)

cat("Results directory:", results_dir, "\n")

for (script_path in scripts_to_run) {
  cat("\n--- Running", script_path, "---\n")
  source(script_path, local = .GlobalEnv)
}

manifest_file <- file.path(results_dir, "results_manifest.txt")
result_files <- list.files(results_dir, full.names = FALSE)
writeLines(c("Generated result files:", sort(result_files)), con = manifest_file)

cat("\nAll scripts completed successfully.\n")
cat("Manifest:", manifest_file, "\n")
