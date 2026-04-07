# Salinity stress statistical analysis
# Input: Salinity_stress_test.csv
# Output: objects in global environment + console summaries

library(dplyr)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
data_dir <- getOption("project_data_dir", getwd())
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

csv_candidates <- c(
  file.path(data_dir, "Salinity_stress_test.csv"),
  "~/Documents/NTNU/Masteroppgave/R_Lengde/Salinity_stress_test.csv"
)

csv_path <- csv_candidates[file.exists(path.expand(csv_candidates))][1]
if (is.na(csv_path)) {
  stop("Could not find Salinity_stress_test.csv in expected locations")
}

Salinity_stress_test <- read.csv(
  path.expand(csv_path),
  sep = ";",
  stringsAsFactors = FALSE,
  strip.white = TRUE,
  check.names = FALSE
)

required_cols <- c("Treatment", "Cone", "FS_alive", "PS_NS_dead", "Total")
missing_cols <- setdiff(required_cols, names(Salinity_stress_test))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "Missing required columns in Salinity_stress_test.csv: %s",
    paste(missing_cols, collapse = ", ")
  ))
}

Salinity_stress_test <- Salinity_stress_test[, required_cols]

sst <- Salinity_stress_test %>%
  mutate(
    Dead = PS_NS_dead,
    FS = FS_alive,
    p_alive = FS / Total,
    pct_alive = 100 * p_alive
  )

sst_main <- sst %>%
  filter(Treatment %in% c("Probiotic", "Control"))

desc <- sst_main %>%
  group_by(Treatment) %>%
  summarise(
    n_cones = n(),
    mean_survival = mean(pct_alive, na.rm = TRUE),
    sd_survival = sd(pct_alive, na.rm = TRUE),
    .groups = "drop"
  )

glm_fit <- glm(cbind(FS, Dead) ~ Treatment, data = sst_main, family = binomial())

ttest_result <- t.test(pct_alive ~ Treatment, data = sst_main)
odds_ci <- exp(confint(glm_fit))
odds_ratio <- exp(coef(glm_fit))

control_check <- sst %>%
  filter(Treatment == "Control_check") %>%
  select(Cone, FS, Dead, Total, pct_alive)

cat("\nDescriptive summary:\n")
print(desc)

cat("\nGLM summary:\n")
print(summary(glm_fit))

cat("\nOdds ratios (exp(coef)):\n")
print(odds_ratio)

cat("\n95% CI for odds ratios (exp(confint)):\n")
print(odds_ci)

cat("\nT-test:\n")
print(ttest_result)

cat("\nControl_check rows:\n")
print(control_check)

odds_table <- data.frame(
  term = names(odds_ratio),
  odds_ratio = as.numeric(odds_ratio),
  ci_lower = as.numeric(odds_ci[, 1]),
  ci_upper = as.numeric(odds_ci[, 2]),
  row.names = NULL
)

write.csv(desc, file.path(results_dir, "salinity_descriptive_summary.csv"), row.names = FALSE)
write.csv(control_check, file.path(results_dir, "salinity_control_check.csv"), row.names = FALSE)
write.csv(odds_table, file.path(results_dir, "salinity_glm_odds_ratios.csv"), row.names = FALSE)

report_file <- file.path(results_dir, "salinity_stats_report.txt")
report_lines <- c(
  "Salinity stress analysis report",
  paste("Source file:", path.expand(csv_path)),
  "",
  "Descriptive summary:",
  capture.output(print(desc)),
  "",
  "GLM summary:",
  capture.output(print(summary(glm_fit))),
  "",
  "Odds ratios (with 95% CI):",
  capture.output(print(odds_table)),
  "",
  "Welch t-test:",
  capture.output(print(ttest_result)),
  "",
  "Control_check rows:",
  capture.output(print(control_check))
)
writeLines(report_lines, con = report_file)

cat("Saved:", file.path(results_dir, "salinity_descriptive_summary.csv"), "\n")
cat("Saved:", file.path(results_dir, "salinity_control_check.csv"), "\n")
cat("Saved:", file.path(results_dir, "salinity_glm_odds_ratios.csv"), "\n")
cat("Saved:", report_file, "\n")
