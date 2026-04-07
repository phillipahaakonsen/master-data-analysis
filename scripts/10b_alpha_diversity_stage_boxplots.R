# Alpha-diversity boxplots by developmental stage with treatment groups
# Input: AlphaDiversityCopProVsCopCtr.csv in project root
# Output: alpha_diversity_stage_boxplots.png in results directory

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
x_axis_text_size <- 18

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
      grepl("^Cop_pro_", sample_id) ~ "Probiotic-treated",
      grepl("^Cop_ctr_", sample_id) ~ "Untreated control",
      TRUE ~ NA_character_
    ),
    day = sub("^(Cop_(?:pro|ctr))_(D[0-9]+)_.+$", "\\2", sample_id, perl = TRUE),
    stage = case_when(
      day == "D16" ~ "Nauplii",
      day == "D21" ~ "Copepodites",
      day == "D24" ~ "Adults",
      TRUE ~ NA_character_
    ),
    treatment = factor(
      treatment,
      levels = c("Probiotic-treated", "Untreated control")
    ),
    day = factor(day, levels = c("D16", "D21", "D24")),
    stage = factor(stage, levels = c("Nauplii", "Copepodites", "Adults")),
    metric = factor(
      metric,
      levels = c(
        "Observed ASV Richness",
        "Exponential Shannon's diversity",
        "Evenness"
      )
    )
  )

if (any(is.na(alpha_plot_df$treatment)) ||
    any(is.na(alpha_plot_df$day)) ||
    any(is.na(alpha_plot_df$stage))) {
  stop("Could not parse treatment/day/stage from one or more sample IDs in ", alpha_path, ".")
}

treatment_fill_colors <- c(
  "Probiotic-treated" = "#20B7EA",
  "Untreated control" = "#FDC300"
)

treatment_outline_colors <- c(
  "Probiotic-treated" = "#148CB6",
  "Untreated control" = "#C89700"
)

stage_axis_labels <- c(
  "Nauplii" = "Nauplii\n(D16)",
  "Copepodites" = "Copepodites\n(D21)",
  "Adults" = "Adults\n(D24)"
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
    aes(x = stage, y = value, fill = treatment, color = treatment)
  ) +
    geom_boxplot(
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
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.72),
      shape = 4,
      size = 3.1,
      stroke = 1.15,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = treatment_fill_colors) +
    scale_color_manual(values = treatment_outline_colors) +
    scale_x_discrete(labels = stage_axis_labels) +
    labs(
      x = NULL,
      y = y_label,
      fill = NULL,
      color = NULL
    ) +
    base_plot_theme
}

p_alpha_stage <- (
  make_metric_plot("Observed ASV Richness", "Observed ASV Richness")
) +
  make_metric_plot("Exponential Shannon's diversity", "Exponential Shannon's diversity") +
  make_metric_plot("Evenness", "Evenness") +
  plot_layout(nrow = 1, guides = "collect") &
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  )

output_file <- file.path(results_dir, "alpha_diversity_stage_boxplots.png")
ggsave(
  filename = output_file,
  plot = p_alpha_stage,
  width = 16,
  height = 6.5,
  dpi = 300
)

cat("Saved plot:", output_file, "\n")
