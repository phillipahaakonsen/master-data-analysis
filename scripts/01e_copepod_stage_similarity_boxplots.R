# Mean Bray-Curtis similarity boxplot for developmental-stage shifts in copepod microbiota
# Input: BrayCurtisSimilarities_copepods.csv in project root
# Output:
#   - copepod_stage_similarity_boxplot_combined.png in results/
#   - copepod_stage_similarity_boxplot_values.csv in results/

library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(readr)

source("scripts/01_similarity_boxplot_helper.R", local = TRUE)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

base_family <- "serif"
base_size <- 18
axis_title_size <- 20
axis_text_size <- 14
x_axis_text_size <- 13
legend_text_size <- 15

comparison_specs <- tibble::tribble(
  ~comparison_id, ~day_1, ~day_2, ~comparison_label,
  "within_d16", "D16", "D16", "Within nauplii\n(D16)",
  "within_d21", "D21", "D21", "Within copepodites\n(D21)",
  "within_d24", "D24", "D24", "Within adults\n(D24)",
  "d16_vs_d21", "D16", "D21", "Nauplii vs. copepodites\n(D16 vs. D21)",
  "d21_vs_d24", "D21", "D24", "Copepodites vs. adults\n(D21 vs. D24)",
  "d16_vs_d24", "D16", "D24", "Nauplii vs. adults\n(D16 vs. D24)"
)

treatment_palette <- tibble::tribble(
  ~treatment, ~fill_color, ~outline_color,
  "Cop_pro", "#20B7EA", "#148CB6",
  "Cop_ctr", "#FDC300", "#C89700"
)

build_stage_mean_group <- function(similarity_matrix, sample_meta, treatment_value, comparison_id, day_1, day_2, comparison_label) {
  day_1_samples <- sample_meta %>%
    filter(treatment == treatment_value, day == day_1) %>%
    pull(sample_id)

  day_2_samples <- sample_meta %>%
    filter(treatment == treatment_value, day == day_2) %>%
    pull(sample_id)

  if (day_1 == day_2) {
    # Within-stage groups average each sample's similarity to the other samples from the same stage.
    mean_data <- compute_mean_similarity_values(
      similarity_matrix = similarity_matrix,
      focal_samples = day_1_samples,
      target_samples = day_1_samples,
      exclude_self = TRUE
    )
  } else {
    # Between-stage groups average each sample's similarity to all samples from the paired stage.
    mean_data <- bind_rows(
      compute_mean_similarity_values(
        similarity_matrix = similarity_matrix,
        focal_samples = day_1_samples,
        target_samples = day_2_samples,
        exclude_self = FALSE
      ),
      compute_mean_similarity_values(
        similarity_matrix = similarity_matrix,
        focal_samples = day_2_samples,
        target_samples = day_1_samples,
        exclude_self = FALSE
      )
    )
  }

  mean_data %>%
    mutate(
      treatment = treatment_value,
      comparison_id = comparison_id,
      comparison_label = comparison_label
    )
}

plot_similarity_boxplot <- function(plot_data) {
  dodge_width <- 0.72

  plot_data <- plot_data %>%
    mutate(
      comparison_label = factor(
        comparison_label,
        levels = comparison_specs$comparison_label
      ),
      treatment = factor(treatment, levels = c("Cop_pro", "Cop_ctr"))
    )

  p_similarity_stage <- ggplot(
    plot_data,
    aes(x = comparison_label, y = similarity, fill = treatment, color = treatment)
  ) +
    geom_boxplot(
      position = position_dodge(width = dodge_width),
      width = 0.34,
      linewidth = 0.8,
      outlier.shape = NA,
      alpha = 0.45
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.08, dodge.width = dodge_width),
      size = 1.9,
      alpha = 0.8,
      shape = 21,
      stroke = 0.25
    ) +
    scale_fill_manual(
      values = setNames(treatment_palette$fill_color, treatment_palette$treatment),
      labels = c("Probiotic-treated", "Untreated control"),
      name = NULL
    ) +
    scale_color_manual(
      values = setNames(treatment_palette$outline_color, treatment_palette$treatment),
      labels = c("Probiotic-treated", "Untreated control"),
      name = NULL
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0.02)
    ) +
    labs(
      x = NULL,
      y = "Mean Bray-Curtis similarity"
    ) +
    theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.title = element_text(size = axis_title_size, color = "black"),
      axis.text = element_text(size = axis_text_size, color = "black"),
      axis.text.x = element_text(size = x_axis_text_size, color = "black", lineheight = 0.95),
      legend.text = element_text(size = legend_text_size),
      legend.position = "top",
      legend.direction = "horizontal",
      plot.margin = margin(12, 16, 10, 10)
    )

  output_file <- file.path(results_dir, "copepod_stage_similarity_boxplot_combined.png")
  ggsave(
    filename = output_file,
    plot = p_similarity_stage,
    width = 12,
    height = 7.5,
    dpi = 300
  )

  cat("Saved combined plot:", output_file, "\n")
}

similarity_matrix <- read_similarity_matrix()
sample_meta <- build_similarity_sample_meta(similarity_matrix)

stage_similarity_data <- pmap_dfr(
  comparison_specs,
  function(comparison_id, day_1, day_2, comparison_label) {
    map_dfr(
      treatment_palette$treatment,
      function(treatment_value) {
        build_stage_mean_group(
          similarity_matrix = similarity_matrix,
          sample_meta = sample_meta,
          treatment_value = treatment_value,
          comparison_id = comparison_id,
          day_1 = day_1,
          day_2 = day_2,
          comparison_label = comparison_label
        )
      }
    )
  }
)

if (nrow(stage_similarity_data) == 0) {
  stop("No developmental-stage mean similarity values were produced from BrayCurtisSimilarities_copepods.csv.")
}

values_output_file <- file.path(results_dir, "copepod_stage_similarity_boxplot_values.csv")
readr::write_csv(
  stage_similarity_data %>%
    mutate(
      treatment_label = recode(
        treatment,
        "Cop_pro" = "Probiotic-treated",
        "Cop_ctr" = "Untreated control"
      )
    ) %>%
    select(
      treatment,
      treatment_label,
      comparison_id,
      comparison_label,
      sample_id,
      n_compared,
      mean_similarity = similarity
    ),
  values_output_file
)
cat("Saved values:", values_output_file, "\n")

plot_similarity_boxplot(stage_similarity_data)
