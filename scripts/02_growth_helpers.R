# Growth helpers and shared objects
# Requires these data frames already in memory:
# - rep_means: treatment, unit, day, mean_length
# - emm_df: day, treatment, emmean, lower.CL, upper.CL

library(dplyr)

if (!exists("rep_means") || !exists("emm_df")) {
  source("scripts/07_growth_inputs_from_length_measurements.R")
}

if (!exists("rep_means")) stop("rep_means is missing from environment")
if (!exists("emm_df")) stop("emm_df is missing from environment")

day_to_x <- function(d) dplyr::case_when(
  d == "Oct15" ~ 15,
  d == "Oct20" ~ 20,
  d == "Oct23" ~ 23,
  TRUE ~ NA_real_
)

rep_means2 <- rep_means %>%
  filter(!is.na(treatment), treatment != "") %>%
  mutate(
    day = factor(as.character(day), levels = c("Oct15", "Oct20", "Oct23")),
    x = day_to_x(as.character(day))
  ) %>%
  arrange(treatment, unit, x)

emm_df2 <- emm_df %>%
  mutate(
    day = factor(as.character(day), levels = c("Oct15", "Oct20", "Oct23")),
    x = day_to_x(as.character(day))
  )

growth_20_23 <- rep_means2 %>%
  filter(x %in% c(20, 23)) %>%
  group_by(treatment, unit) %>%
  summarise(
    len_20 = mean(mean_length[x == 20], na.rm = TRUE),
    len_23 = mean(mean_length[x == 23], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(rate_um_per_day = (len_23 - len_20) / 3) %>%
  group_by(treatment) %>%
  summarise(rate = mean(rate_um_per_day, na.rm = TRUE), .groups = "drop")

rate_control <- growth_20_23$rate[growth_20_23$treatment == "control"]
rate_probiotic <- growth_20_23$rate[growth_20_23$treatment == "probiotic"]

day_labels <- c("Oct15" = "Oct. 15th", "Oct20" = "Oct. 20th", "Oct23" = "Oct. 23rd")
