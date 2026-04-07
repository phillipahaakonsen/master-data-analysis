# PCoA plot helper for PAST score and summary exports

library(dplyr)
library(ggplot2)

build_pcoa_plot <- function(
  scores_file,
  summary_file,
  day_label,
  probiotic_color = "#20B7EA",
  control_color = "#FDC300",
  x_limits = NULL,
  y_limits = NULL,
  x_breaks = NULL,
  y_breaks = NULL,
  x_label_accuracy = NULL,
  y_label_accuracy = NULL,
  show_legend = TRUE
) {
  results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  base_family <- "serif"
  base_size <- 18
  axis_title_size <- 22
  axis_text_size <- 18
  legend_text_size <- 18
  pcoa_raw <- read.csv(
    scores_file,
    sep = ";",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(pcoa_raw) >= 1) {
    names(pcoa_raw)[1] <- "sample_id"
  }

  pcoa_summary <- read.csv(
    summary_file,
    sep = ";",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(pcoa_summary) >= 3) {
    names(pcoa_summary)[1:3] <- c("Axis", "Eigenvalue", "Percent")
  }

  coord1_percent <- pcoa_summary %>%
    filter(Axis == 1) %>%
    pull(Percent)

  coord2_percent <- pcoa_summary %>%
    filter(Axis == 2) %>%
    pull(Percent)

  coord1_eigenvalue <- pcoa_summary %>%
    filter(Axis == 1) %>%
    pull(Eigenvalue)

  coord2_eigenvalue <- pcoa_summary %>%
    filter(Axis == 2) %>%
    pull(Eigenvalue)

  coord1_label <- if (length(coord1_percent) == 1) {
    sprintf("Coordinate 1 (%.3f%%)", coord1_percent)
  } else {
    "Coordinate 1"
  }

  coord2_label <- if (length(coord2_percent) == 1) {
    sprintf("Coordinate 2 (%.3f%%)", coord2_percent)
  } else {
    "Coordinate 2"
  }

  coord1_scale <- if (length(coord1_eigenvalue) == 1 && coord1_eigenvalue > 0) sqrt(coord1_eigenvalue) else 1
  coord2_scale <- if (length(coord2_eigenvalue) == 1 && coord2_eigenvalue > 0) sqrt(coord2_eigenvalue) else 1

  pro_pattern <- sprintf("^Cop_pro_%s", day_label)
  ctr_pattern <- sprintf("^Cop_ctr_%s", day_label)
  replicate_pattern <- sprintf("^(?:Cop_pro_%s|Cop_ctr_%s)_([1-4])_.*$", day_label, day_label)

  pcoa_plot <- pcoa_raw %>%
    transmute(
      sample_id = sample_id,
      coord1 = as.numeric(`Coord 1`) * coord1_scale,
      coord2 = as.numeric(`Coord 2`) * coord2_scale,
      sample_type = case_when(
        grepl(pro_pattern, sample_id) ~ sprintf("Cop_pro_%s", day_label),
        grepl(ctr_pattern, sample_id) ~ sprintf("Cop_ctr_%s", day_label),
        TRUE ~ NA_character_
      ),
      replicate_unit = sub(replicate_pattern, "\\1", sample_id)
    ) %>%
    filter(!is.na(sample_type)) %>%
    mutate(replicate_unit = factor(replicate_unit, levels = c("1", "2", "3", "4")))

  hull_plot <- pcoa_plot %>%
    group_by(sample_type) %>%
    slice(chull(coord1, coord2)) %>%
    ungroup()

  if (is.null(x_limits)) {
    x_limits <- c(floor(min(pcoa_plot$coord1) * 10) / 10, ceiling(max(pcoa_plot$coord1) * 10) / 10)
  }
  if (is.null(y_limits)) {
    y_limits <- c(floor(min(pcoa_plot$coord2) * 10) / 10, ceiling(max(pcoa_plot$coord2) * 10) / 10)
  }
  if (is.null(x_breaks)) {
    x_breaks <- seq(x_limits[1], x_limits[2], by = 0.1)
  }
  if (is.null(y_breaks)) {
    y_breaks <- seq(y_limits[1], y_limits[2], by = 0.1)
  }

  x_labels <- waiver()
  y_labels <- waiver()
  if (!is.null(x_label_accuracy)) {
    x_labels <- scales::label_number(accuracy = x_label_accuracy)
  } else {
    x_labels <- scales::label_number(accuracy = 0.01)
  }
  if (!is.null(y_label_accuracy)) {
    y_labels <- scales::label_number(accuracy = y_label_accuracy)
  } else {
    y_labels <- scales::label_number(accuracy = 0.01)
  }

  pro_label <- sprintf("Cop_pro_%s", day_label)
  ctr_label <- sprintf("Cop_ctr_%s", day_label)

  ggplot(
    pcoa_plot,
    aes(x = coord1, y = coord2, color = sample_type, fill = sample_type, shape = replicate_unit, size = replicate_unit)
  ) +
    geom_polygon(
      data = hull_plot,
      aes(group = sample_type),
      alpha = 0.25,
      linewidth = 1
    ) +
    geom_point(stroke = 1.4) +
    scale_color_manual(
      values = setNames(c(probiotic_color, control_color), c(pro_label, ctr_label)),
      breaks = c(pro_label, ctr_label),
      labels = setNames(c(pro_label, ctr_label), c(pro_label, ctr_label))
    ) +
    scale_fill_manual(
      values = setNames(c(probiotic_color, control_color), c(pro_label, ctr_label)),
      breaks = c(pro_label, ctr_label),
      labels = setNames(c(pro_label, ctr_label), c(pro_label, ctr_label))
    ) +
    scale_shape_manual(
      values = c("1" = 1, "2" = 8, "3" = 18, "4" = 0),
      guide = "none"
    ) +
    scale_size_manual(
      values = c("1" = 6.1, "2" = 6.8, "3" = 6.8, "4" = 6.1),
      guide = "none"
    ) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks, labels = y_labels, expand = c(0, 0)) +
    labs(
      x = coord1_label,
      y = coord2_label,
      color = NULL,
      fill = NULL
    ) +
    theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size, color = "black"),
      legend.text = element_text(size = legend_text_size),
      legend.position = if (show_legend) "right" else "none"
    )
}

make_pcoa_plot <- function(
  scores_file,
  summary_file,
  day_label,
  output_name,
  probiotic_color = "#20B7EA",
  control_color = "#FDC300",
  x_limits = NULL,
  y_limits = NULL,
  x_breaks = NULL,
  y_breaks = NULL,
  x_label_accuracy = NULL,
  y_label_accuracy = NULL,
  show_legend = TRUE
) {
  results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  p_pcoa <- build_pcoa_plot(
    scores_file = scores_file,
    summary_file = summary_file,
    day_label = day_label,
    probiotic_color = probiotic_color,
    control_color = control_color,
    x_limits = x_limits,
    y_limits = y_limits,
    x_breaks = x_breaks,
    y_breaks = y_breaks,
    x_label_accuracy = x_label_accuracy,
    y_label_accuracy = y_label_accuracy,
    show_legend = show_legend
  )

  print(p_pcoa)

  output_file <- file.path(results_dir, output_name)
  ggsave(
    filename = output_file,
    plot = p_pcoa,
    width = 10,
    height = 8,
    dpi = 300
  )
  cat("Saved plot:", output_file, "\n")
}
