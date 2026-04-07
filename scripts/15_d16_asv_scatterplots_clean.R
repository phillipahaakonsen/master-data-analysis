# D16 ASV relative abundance scatterplots with clean manual labels
# Input:
#   - results/asv_relative_abundance_summary_D16.csv
# Output:
#   - asv_relative_abundance_scatterplots_D16_clean.png in results/

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(patchwork)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

summary_path <- file.path(results_dir, "asv_relative_abundance_summary_D16.csv")
output_path <- file.path(results_dir, "asv_relative_abundance_scatterplots_D16_clean.png")

if (!file.exists(summary_path)) {
  stop("Missing ", summary_path, ". Run scripts/11_asv_relative_abundance_scatterplots.R first.")
}

base_family <- "Times New Roman"
x_axis_limit_percent <- 15

sample_type_colors <- c(
  "Cop_pro_D16" = "#20B7EA",
  "Cop_ctr_D16" = "#FDC300"
)

make_pretty_axis_limit <- function(max_value) {
  if (!is.finite(max_value) || max_value <= 0) {
    return(1)
  }

  if (max_value <= 10) {
    step <- 2
  } else if (max_value <= 25) {
    step <- 5
  } else if (max_value <= 50) {
    step <- 10
  } else {
    step <- 20
  }

  ceiling(max_value / step) * step
}

make_axis_breaks <- function(axis_limit) {
  if (!is.finite(axis_limit) || axis_limit <= 0) {
    return(c(0, 1))
  }

  if (axis_limit <= 10) {
    step <- 2
  } else if (axis_limit <= 25) {
    step <- 5
  } else if (axis_limit <= 50) {
    step <- 10
  } else {
    step <- 20
  }

  seq(0, axis_limit, by = step)
}

extract_genus <- function(taxonomy) {
  taxonomy_chr <- as.character(taxonomy)
  genus <- str_match(taxonomy_chr, "(^|[,;])\\s*g:([^,;]+)")[, 3]
  ifelse(is.na(genus) | genus == "", NA_character_, genus)
}

format_asv_label <- function(asv_id, taxonomy) {
  display_asv_id <- str_replace(as.character(asv_id), regex("^zotu", ignore_case = TRUE), "ASV")
  genus <- extract_genus(taxonomy)
  ifelse(is.na(genus), display_asv_id, paste(display_asv_id, genus, sep = " - "))
}

build_label_df <- function(plot_df, placements) {
  plot_df %>%
    filter(ASV_ID %in% placements$ASV_ID) %>%
    transmute(
      ASV_ID,
      x_value,
      y_value,
      plot_label
    ) %>%
    left_join(placements, by = "ASV_ID")
}

make_manual_plot <- function(plot_df, label_df, x_label, y_label, point_color, y_axis_limit) {
  x_breaks <- make_axis_breaks(x_axis_limit_percent)
  y_breaks <- make_axis_breaks(y_axis_limit)

  ggplot(plot_df, aes(x = x_value, y = y_value)) +
    geom_point(
      shape = 21,
      size = 4.2,
      stroke = 0.6,
      fill = point_color,
      color = point_color,
      alpha = 0.85
    ) +
    geom_segment(
      data = label_df,
      aes(x = x_value, y = y_value, xend = line_x, yend = line_y),
      inherit.aes = FALSE,
      color = "grey35",
      linewidth = 0.28,
      na.rm = TRUE
    ) +
    geom_text(
      data = label_df,
      aes(x = label_x, y = label_y, label = plot_label),
      inherit.aes = FALSE,
      family = base_family,
      size = 6.2,
      hjust = label_df$label_hjust,
      vjust = label_df$label_vjust
    ) +
    coord_cartesian(
      xlim = c(0, x_axis_limit_percent),
      ylim = c(-1.8, y_axis_limit),
      expand = TRUE
    ) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(
      title = NULL,
      x = "Mean relative abundance (%) in Prob",
      y = y_label
    ) +
    theme_classic(base_size = 16, base_family = base_family) +
    theme(
      axis.title = element_text(color = "black", size = 22),
      axis.text = element_text(color = "black", size = 22)
    )
}

summary_tbl <- read_csv(summary_path, show_col_types = FALSE)

y_axis_limit <- make_pretty_axis_limit(
  max(c(summary_tbl$mean_Cop_pro, summary_tbl$mean_Cop_ctr), na.rm = TRUE)
)

plot_prob_vs_cop_pro <- summary_tbl %>%
  transmute(
    ASV_ID,
    taxonomy,
    x_value = mean_Prob,
    y_value = mean_Cop_pro,
    plot_label = format_asv_label(ASV_ID, taxonomy)
  )

plot_prob_vs_cop_ctr <- summary_tbl %>%
  transmute(
    ASV_ID,
    taxonomy,
    x_value = mean_Prob,
    y_value = mean_Cop_ctr,
    plot_label = format_asv_label(ASV_ID, taxonomy)
  )

# Left panel: remove ASV1 and ASV6, but keep the two x-axis Lactobacillus points.
cop_pro_label_positions <- tibble(
  ASV_ID = c("Zotu2", "Zotu8", "Zotu12", "Zotu70", "Zotu60", "Zotu29", "Zotu17", "Zotu13"),
  label_x = c(0.78, 0.78, 0.62, 2.10, 3.80, 5.10, 7.55, 10.95),
  label_y = c(17.6, 15.0, 7.65, -0.20, -1.20, 0.62, 1.42, 1.00),
  line_x = c(0.70, 0.62, 0.30, 2.95, 3.92, 5.20, 8.38, 10.28),
  line_y = c(17.3, 14.8, 7.35, -0.02, -0.55, 0.28, 0.55, 0.36),
  label_hjust = c(0, 0, 0, 0, 0, 0, 0, 0),
  label_vjust = c(0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5)
)

cop_ctr_label_positions <- tibble(
  ASV_ID = c("Zotu5", "Zotu9", "Zotu2", "Zotu21", "Zotu60", "Zotu29", "Zotu17", "Zotu13"),
  label_x = c(0.38, 0.38, 0.38, 0.38, 2.30, 5.55, 7.95, 9.80),
  label_y = c(28.0, 15.4, 10.8, 5.7, 1.95, 0.82, -0.62, 2.80),
  line_x = c(0.16, 0.10, 0.16, 0.09, 2.62, 5.42, 8.30, 9.58),
  line_y = c(27.6, 15.0, 10.45, 5.25, 1.15, 0.30, -0.22, 2.36),
  label_hjust = c(0, 0, 0, 0, 0, 0, 0.5, 0),
  label_vjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5)
)

cop_pro_labels <- build_label_df(plot_prob_vs_cop_pro, cop_pro_label_positions)
cop_ctr_labels <- build_label_df(plot_prob_vs_cop_ctr, cop_ctr_label_positions)

p_cop_pro <- make_manual_plot(
  plot_df = plot_prob_vs_cop_pro,
  label_df = cop_pro_labels,
  x_label = "Mean relative abundance (%) in Prob",
  y_label = "Mean relative abundance (%) in Cop_pro",
  point_color = sample_type_colors[["Cop_pro_D16"]],
  y_axis_limit = y_axis_limit
)

p_cop_ctr <- make_manual_plot(
  plot_df = plot_prob_vs_cop_ctr,
  label_df = cop_ctr_labels,
  x_label = "Mean relative abundance (%) in Prob",
  y_label = "Mean relative abundance (%) in Cop_ctr",
  point_color = sample_type_colors[["Cop_ctr_D16"]],
  y_axis_limit = y_axis_limit
)

combined_plot <- p_cop_pro + p_cop_ctr + plot_layout(ncol = 2)

ggsave(
  filename = output_path,
  plot = combined_plot,
  width = 16,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("Saved clean D16 plot:", output_path, "\n")
