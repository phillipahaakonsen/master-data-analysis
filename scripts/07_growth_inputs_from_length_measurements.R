# Build growth-model inputs and appendix-ready outputs from length_measurements.csv
# Output objects in the global environment:
# - length_raw
# - growth_counts
# - rep_means
# - growth_model
# - growth_anova
# - growth_fixed_effects
# - growth_random_effects
# - emm_df
# - growth_pairwise
# - growth_rate_by_unit
# - growth_rate_summary
# - growth_rate_ttest_20_23

library(dplyr)
library(emmeans)
library(lmerTest)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

input_file <- "length_measurements.csv"
if (!file.exists(input_file)) {
  stop("Missing input file: length_measurements.csv in project root")
}

length_raw <- read.csv(input_file, sep = ";", stringsAsFactors = FALSE) %>%
  mutate(
    treatment = trimws(treatment),
    day = trimws(day),
    unit = trimws(unit)
  ) %>%
  filter(
    !is.na(length),
    !is.na(treatment), treatment != "",
    !is.na(day), day != "",
    !is.na(unit), unit != ""
  ) %>%
  mutate(
    treatment = case_when(
      tolower(treatment) %in% c("control", "ctr") ~ "control",
      tolower(treatment) %in% c("probiotic", "pb", "pro") ~ "probiotic",
      TRUE ~ NA_character_
    ),
    day = case_when(
      day %in% c("Oct15", "Oct20", "Oct23") ~ day,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(treatment), !is.na(day)) %>%
  mutate(
    treatment = factor(treatment, levels = c("control", "probiotic")),
    day = factor(day, levels = c("Oct15", "Oct20", "Oct23")),
    unit = factor(unit)
  )

growth_counts <- length_raw %>%
  count(treatment, day, unit, name = "n_individuals")

growth_sampling_summary <- growth_counts %>%
  group_by(treatment, day) %>%
  summarise(
    n_units = n(),
    n_individuals = sum(n_individuals),
    .groups = "drop"
  )

rep_means <- length_raw %>%
  group_by(treatment, unit, day) %>%
  summarise(mean_length = mean(length, na.rm = TRUE), .groups = "drop")

calc_growth_interval <- function(data, start_day, end_day, interval_days) {
  start_df <- data %>%
    filter(day == start_day) %>%
    transmute(treatment, unit, len_start = mean_length)

  end_df <- data %>%
    filter(day == end_day) %>%
    transmute(treatment, unit, len_end = mean_length)

  inner_join(start_df, end_df, by = c("treatment", "unit")) %>%
    mutate(
      day_start = start_day,
      day_end = end_day,
      interval_days = interval_days,
      interval_label = paste0(start_day, "_", end_day),
      rate_um_per_day = (len_end - len_start) / interval_days
    ) %>%
    select(
      treatment,
      unit,
      day_start,
      day_end,
      interval_days,
      interval_label,
      len_start,
      len_end,
      rate_um_per_day
    )
}

growth_rate_by_unit <- bind_rows(
  calc_growth_interval(rep_means, "Oct15", "Oct20", 5),
  calc_growth_interval(rep_means, "Oct20", "Oct23", 3)
)

growth_rate_summary <- growth_rate_by_unit %>%
  group_by(interval_label, day_start, day_end, interval_days, treatment) %>%
  summarise(
    n_units = n(),
    mean_rate = mean(rate_um_per_day, na.rm = TRUE),
    sd_rate = if (n() > 1) sd(rate_um_per_day, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

growth_rate_20_23 <- growth_rate_by_unit %>%
  filter(interval_label == "Oct20_Oct23")

growth_rate_ttest_obj <- t.test(rate_um_per_day ~ treatment, data = growth_rate_20_23, var.equal = FALSE)

growth_rate_ttest_20_23 <- data.frame(
  interval_label = "Oct20_Oct23",
  method = unname(growth_rate_ttest_obj$method),
  statistic = unname(growth_rate_ttest_obj$statistic),
  parameter = unname(growth_rate_ttest_obj$parameter),
  p_value = growth_rate_ttest_obj$p.value,
  conf_low = growth_rate_ttest_obj$conf.int[1],
  conf_high = growth_rate_ttest_obj$conf.int[2],
  mean_control = unname(growth_rate_ttest_obj$estimate["mean in group control"]),
  mean_probiotic = unname(growth_rate_ttest_obj$estimate["mean in group probiotic"]),
  stringsAsFactors = FALSE
)

old_contrasts <- getOption("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))
on.exit(options(contrasts = old_contrasts), add = TRUE)

emm_options(lmer.df = "satterthwaite")

# Model individuals while accounting for clustering within replicate unit/day samples.
growth_model <- lmer(length ~ treatment * day + (1 | unit/day), data = length_raw, REML = TRUE)

growth_anova <- as.data.frame(anova(growth_model, type = 3)) %>%
  mutate(term = row.names(.), .before = 1) %>%
  rename(
    sum_sq = `Sum Sq`,
    mean_sq = `Mean Sq`,
    num_df = NumDF,
    den_df = DenDF,
    f_value = `F value`,
    p_value = `Pr(>F)`
  ) %>%
  mutate(row.names = NULL)

growth_fixed_effects <- as.data.frame(coef(summary(growth_model))) %>%
  mutate(term = row.names(.), .before = 1) %>%
  rename(
    estimate = Estimate,
    std_error = `Std. Error`,
    df = df,
    statistic = `t value`,
    p_value = `Pr(>|t|)`
  ) %>%
  mutate(row.names = NULL)

growth_random_effects <- as.data.frame(VarCorr(growth_model)) %>%
  select(grp, var1, var2, vcov, sdcor) %>%
  rename(
    group = grp,
    term = var1,
    term2 = var2,
    variance = vcov,
    std_dev = sdcor
  )

growth_emmeans <- emmeans(growth_model, ~ treatment | day)

emm_df <- as.data.frame(summary(growth_emmeans, infer = c(TRUE, FALSE))) %>%
  mutate(
    treatment = as.character(treatment),
    day = as.character(day)
  )

growth_pairwise <- as.data.frame(
  summary(pairs(growth_emmeans, adjust = "holm"), infer = c(TRUE, TRUE))
) %>%
  rename(
    std_error = SE,
    statistic = t.ratio,
    p_value = p.value,
    conf_low = lower.CL,
    conf_high = upper.CL
  )

write.csv(growth_counts, file.path(results_dir, "growth_length_counts.csv"), row.names = FALSE)
write.csv(growth_sampling_summary, file.path(results_dir, "growth_sampling_summary.csv"), row.names = FALSE)
write.csv(rep_means, file.path(results_dir, "growth_rep_means.csv"), row.names = FALSE)
write.csv(growth_anova, file.path(results_dir, "growth_type3_anova.csv"), row.names = FALSE)
write.csv(growth_fixed_effects, file.path(results_dir, "growth_fixed_effects.csv"), row.names = FALSE)
write.csv(growth_random_effects, file.path(results_dir, "growth_random_effects.csv"), row.names = FALSE)
write.csv(emm_df, file.path(results_dir, "growth_emm_df.csv"), row.names = FALSE)
write.csv(growth_pairwise, file.path(results_dir, "growth_pairwise_by_day.csv"), row.names = FALSE)
write.csv(growth_rate_by_unit, file.path(results_dir, "growth_rate_by_unit.csv"), row.names = FALSE)
write.csv(growth_rate_summary, file.path(results_dir, "growth_rate_summary.csv"), row.names = FALSE)
write.csv(growth_rate_ttest_20_23, file.path(results_dir, "growth_rate_welch_ttest_20_23.csv"), row.names = FALSE)

report_file <- file.path(results_dir, "growth_stats_report.txt")
report_lines <- c(
  "Growth analysis report",
  paste("Source file:", normalizePath(input_file)),
  "Primary model: lmer(length ~ treatment * day + (1 | unit/day), REML = TRUE)",
  paste("Singular fit:", lme4::isSingular(growth_model)),
  "Design note: October 15 includes one replicate culture unit per treatment.",
  "Interpret treatment comparisons for that day cautiously.",
  "",
  "Sampling summary by treatment/day:",
  capture.output(print(growth_sampling_summary)),
  "",
  "Counts by treatment/day/unit:",
  capture.output(print(growth_counts)),
  "",
  "Replicate-level means used for trajectory plots:",
  capture.output(print(rep_means)),
  "",
  "Mixed-model summary:",
  capture.output(print(summary(growth_model))),
  "",
  "Type III ANOVA:",
  capture.output(print(growth_anova)),
  "",
  "Fixed-effect estimates:",
  capture.output(print(growth_fixed_effects)),
  "",
  "Random-effect estimates:",
  capture.output(print(growth_random_effects)),
  "",
  "Estimated marginal means:",
  capture.output(print(emm_df)),
  "",
  "Pairwise treatment contrasts within day:",
  capture.output(print(growth_pairwise)),
  "",
  "Descriptive growth rate by unit:",
  capture.output(print(growth_rate_by_unit)),
  "",
  "Descriptive growth rate summary:",
  capture.output(print(growth_rate_summary)),
  "",
  "Welch t-test for Oct20-Oct23 growth rate:",
  capture.output(print(growth_rate_ttest_20_23)),
  "",
  "Full Welch t-test output for Oct20-Oct23 growth rate:",
  capture.output(print(growth_rate_ttest_obj))
)
writeLines(report_lines, con = report_file)

cat("Created objects: growth_counts, rep_means, growth_model, growth_anova, growth_fixed_effects, growth_random_effects, emm_df, growth_pairwise, growth_rate_by_unit, growth_rate_summary, growth_rate_ttest_20_23\n")
cat("length_raw rows:", nrow(length_raw), "\n")
cat("growth_counts rows:", nrow(growth_counts), "\n")
cat("growth_sampling_summary rows:", nrow(growth_sampling_summary), "\n")
cat("rep_means rows:", nrow(rep_means), "\n")
cat("emm_df rows:", nrow(emm_df), "\n")
cat("growth_rate_by_unit rows:", nrow(growth_rate_by_unit), "\n")
cat("Saved:", file.path(results_dir, "growth_length_counts.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_sampling_summary.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_rep_means.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_type3_anova.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_fixed_effects.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_random_effects.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_emm_df.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_pairwise_by_day.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_rate_by_unit.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_rate_summary.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_rate_welch_ttest_20_23.csv"), "\n")
cat("Saved:", report_file, "\n")
