# Build rep_means and emm_df from length_measurements.csv
# Output objects in the global environment:
# - length_raw
# - rep_means
# - growth_model
# - emm_df

library(dplyr)
library(emmeans)

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

rep_means <- length_raw %>%
  group_by(treatment, unit, day) %>%
  summarise(mean_length = mean(length, na.rm = TRUE), .groups = "drop")

growth_model <- lm(mean_length ~ treatment * day, data = rep_means)

emm_df <- as.data.frame(emmeans(growth_model, ~ treatment | day)) %>%
  mutate(
    treatment = as.character(treatment),
    day = as.character(day)
  )

write.csv(rep_means, file.path(results_dir, "growth_rep_means.csv"), row.names = FALSE)
write.csv(emm_df, file.path(results_dir, "growth_emm_df.csv"), row.names = FALSE)

cat("Created objects: rep_means and emm_df\n")
cat("rep_means rows:", nrow(rep_means), "\n")
cat("emm_df rows:", nrow(emm_df), "\n")
cat("Saved:", file.path(results_dir, "growth_rep_means.csv"), "\n")
cat("Saved:", file.path(results_dir, "growth_emm_df.csv"), "\n")
