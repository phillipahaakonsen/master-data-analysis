# Alpha-diversity boxplots for Cop_pro vs Cop_ctr across D16, D21, and D24
# Input: AlphaDiversityCopProVsCopCtr.csv in project root
# Output: alpha_diversity_boxplots.png in results directory

library(dplyr)
library(ggplot2)
library(patchwork)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
data_dir <- getOption("project_data_dir", getwd())
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

base_family <- "serif"
base_size <- 18
axis_text_size <- 15
legend_text_size <- 15
axis_title_size <- 20
axis_line_width <- 0.45
axis_tick_width <- 0.45
x_axis_text_size <- 20

alpha_path <- file.path(data_dir, "AlphaDiversityCopProVsCopCtr.csv")

alpha_raw <- read.delim(
  alpha_path,
  sep = ";",
  row.names = 1,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

required_metrics <- c("Taxa_S", "Shannon_H", "Evenness_e^H/S")
missing_metrics <- setdiff(required_metrics, rownames(alpha_raw))

if (length(missing_metrics) > 0) {
  stop(
    "Missing required alpha-diversity rows in ",
    alpha_path,
    ": ",
    paste(missing_metrics, collapse = ", ")
  )
}

alpha_samples <- as.data.frame(t(alpha_raw[required_metrics, , drop = FALSE]))
alpha_samples[] <- lapply(alpha_samples, as.numeric)
alpha_samples <- data.frame(
  sample_id = rownames(alpha_samples),
  alpha_samples,
  row.names = NULL,
  check.names = FALSE
)
alpha_samples$exp_shannon <- exp(alpha_samples$Shannon_H)

alpha_plot_df <- bind_rows(
  data.frame(
    sample_id = alpha_samples$sample_id,
    metric = "Observed ASV Richness",
    value = alpha_samples$Taxa_S,
    stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = alpha_samples$sample_id,
    metric = "Exponential Shannon's diversity",
    value = exp(alpha_samples$Shannon_H),
    stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = alpha_samples$sample_id,
    metric = "Evenness",
    value = alpha_samples[["Evenness_e^H/S"]],
    stringsAsFactors = FALSE
  )
) %>%
  mutate(
    treatment = case_when(
      grepl("^Cop_pro_", sample_id) ~ "Cop_pro",
      grepl("^Cop_ctr_", sample_id) ~ "Cop_ctr",
      TRUE ~ NA_character_
    ),
    day = sub("^(Cop_(?:pro|ctr))_(D[0-9]+)_.+$", "\\2", sample_id, perl = TRUE),
    treatment = factor(treatment, levels = c("Cop_pro", "Cop_ctr")),
    day = factor(day, levels = c("D16", "D21", "D24")),
    metric = factor(
      metric,
      levels = c(
        "Observed ASV Richness",
        "Exponential Shannon's diversity",
        "Evenness"
      )
    )
  )

if (any(is.na(alpha_plot_df$treatment)) || any(is.na(alpha_plot_df$day))) {
  stop("Could not parse treatment/day from one or more sample IDs in ", alpha_path, ".")
}

metric_column_map <- c(
  "Observed ASV Richness" = "Taxa_S",
  "Exponential Shannon's diversity" = "exp_shannon",
  "Evenness" = "Evenness_e^H/S"
)

build_welch_stats_row <- function(day_label, metric_label, value_column) {
  day_df <- alpha_samples[
    grepl(paste0("^Cop_(?:pro|ctr)_", day_label, "_"), alpha_samples$sample_id),
    ,
    drop = FALSE
  ]

  cop_pro_values <- day_df[grepl("^Cop_pro_", day_df$sample_id), value_column]
  cop_ctr_values <- day_df[grepl("^Cop_ctr_", day_df$sample_id), value_column]
  welch_result <- t.test(cop_pro_values, cop_ctr_values, var.equal = FALSE)

  data.frame(
    day = day_label,
    metric = metric_label,
    n_Cop_pro = length(cop_pro_values),
    n_Cop_ctr = length(cop_ctr_values),
    mean_Cop_pro = round(mean(cop_pro_values), 4),
    mean_Cop_ctr = round(mean(cop_ctr_values), 4),
    median_Cop_pro = round(median(cop_pro_values), 4),
    median_Cop_ctr = round(median(cop_ctr_values), 4),
    p_value_welch = round(welch_result$p.value, 4),
    stringsAsFactors = FALSE
  )
}

alpha_welch_stats <- do.call(
  rbind,
  lapply(levels(alpha_plot_df$day), function(day_label) {
    do.call(
      rbind,
      lapply(names(metric_column_map), function(metric_label) {
        build_welch_stats_row(
          day_label = day_label,
          metric_label = metric_label,
          value_column = unname(metric_column_map[[metric_label]])
        )
      })
    )
  })
)

cat("\nAlpha-diversity day-by-day Welch t-tests:\n")
print(alpha_welch_stats, row.names = FALSE)
cat("\n")

alpha_stats_file <- file.path(results_dir, "alpha_diversity_welch_stats.csv")
write.csv(alpha_welch_stats, file = alpha_stats_file, row.names = FALSE)
cat("Saved Welch stats:", alpha_stats_file, "\n")
cat("\n")

day_colors <- c(
  "D16" = "#FFD84A",
  "D21" = "#90C95A",
  "D24" = "#E58FA2"
)

outline_colors <- c(
  "D16" = "#C29400",
  "D21" = "#5E9F2E",
  "D24" = "#B25E78"
)

base_plot_theme <- theme_classic(base_size = base_size, base_family = base_family) +
  theme(
    axis.text = element_text(size = axis_text_size, color = "black"),
    axis.text.x = element_text(size = x_axis_text_size, color = "black", lineheight = 0.95),
    axis.title = element_text(size = axis_title_size, color = "black", face = "plain"),
    axis.line = element_line(linewidth = axis_line_width, color = "grey30"),
    axis.ticks = element_line(linewidth = axis_tick_width, color = "grey30"),
    legend.text = element_text(size = legend_text_size),
    legend.position = "top",
    panel.spacing = grid::unit(1.2, "lines")
  )

make_metric_plot <- function(metric_name, y_label) {
  ggplot(
    alpha_plot_df %>% filter(metric == metric_name),
    aes(x = treatment, y = value, fill = day, color = day)
  ) +
    geom_boxplot(
      aes(group = interaction(treatment, day)),
      position = position_dodge(width = 0.72),
      width = 0.62,
      linewidth = 0.9,
      alpha = 0.45,
      outlier.shape = NA
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72),
      shape = 21,
      fill = "white",
      stroke = 1,
      size = 2.3,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    stat_summary(
      aes(group = interaction(treatment, day)),
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.72),
      shape = 4,
      size = 3.1,
      stroke = 1.15,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = day_colors) +
    scale_color_manual(values = outline_colors, guide = "none") +
    scale_x_discrete(
      labels = c(
        "Cop_pro" = "Cop_pro",
        "Cop_ctr" = "Cop_ctr"
      )
    ) +
    labs(
      x = NULL,
      y = y_label,
      fill = NULL
    ) +
    base_plot_theme
}

p_alpha <- (
  make_metric_plot("Observed ASV Richness", "Observed ASV Richness")
) +
  make_metric_plot("Exponential Shannon's diversity", "Exponential Shannon's diversity") +
  make_metric_plot("Evenness", "Evenness") +
  plot_layout(nrow = 1, guides = "collect") &
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  )

output_file <- file.path(results_dir, "alpha_diversity_boxplots.png")
ggsave(
  filename = output_file,
  plot = p_alpha,
  width = 15,
  height = 6.5,
  dpi = 300
)

cat("Saved plot:", output_file, "\n")
