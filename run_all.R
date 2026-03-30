# Run all project analyses in order and save artifacts to results/

results_dir <- file.path(getwd(), "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
options(project_results_dir = results_dir)

scripts_to_run <- c(
  "scripts/01_d16_similarity_boxplot.R",
  "scripts/02_growth_helpers.R",
  "scripts/03_growth_plot_discrete.R",
  "scripts/04_growth_plot_numeric.R",
  "scripts/05_salinity_stress_stats.R",
  "scripts/06_salinity_survival_plot.R"
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
