# D24 ASV relative abundance scatterplots with ASV29 placed below its point
# Input:
#   - results/asv_relative_abundance_summary_D24.csv
# Output:
#   - asv_relative_abundance_scatterplots_D24_asv29_below.png in results/

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(patchwork)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

summary_path <- file.path(results_dir, "asv_relative_abundance_summary_D24.csv")
output_path <- file.path(results_dir, "asv_relative_abundance_scatterplots_D24_asv29_below.png")

if (!file.exists(summary_path)) {
  stop(
    "Missing ", summary_path, ". Run scripts/11_asv_relative_abundance_scatterplots.R first."
  )
}

base_family <- "Times New Roman"
label_max_n <- 8
label_outlier_percent <- 2
x_axis_limit_percent <- 15

sample_type_colors <- c(
  "Cop_pro_D24" = "#1624F2",
  "Cop_ctr_D24" = "#FF6200",
  "Prob_D24" = "#7FBF00"
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

make_comparison_plot <- function(summary_tbl, x_col, y_col, x_label, y_label, title_text, point_color, y_axis_limit) {
  plot_df <- summary_tbl %>%
    mutate(
      x_value = .data[[x_col]],
      y_value = .data[[y_col]],
      combined_abundance = x_value + y_value,
      plot_label = format_asv_label(ASV_ID, taxonomy)
    )

  label_df <- plot_df %>%
    filter(pmax(x_value, y_value) >= label_outlier_percent) %>%
    arrange(desc(combined_abundance), ASV_ID) %>%
    slice_head(n = label_max_n)

  if (nrow(label_df) == 0) {
    label_df <- plot_df %>%
      arrange(desc(combined_abundance), ASV_ID) %>%
      slice_head(n = min(3, n()))
  }

  label_df <- label_df %>%
    filter(!ASV_ID %in% c("Zotu18", "Zotu24"))

  manual_label_df <- label_df %>%
    slice(0) %>%
    mutate(
      label_x = numeric(),
      label_y = numeric(),
      label_hjust = numeric(),
      label_vjust = numeric()
    )

  if (title_text == "D24: Prob vs Cop_pro") {
    manual_label_df <- label_df %>%
      filter(ASV_ID %in% c("Zotu29", "Zotu2")) %>%
      mutate(
        label_x = case_when(
          ASV_ID == "Zotu29" ~ x_value,
          ASV_ID == "Zotu2" ~ x_value + 0.35,
          TRUE ~ x_value
        ),
        label_y = case_when(
          ASV_ID == "Zotu29" ~ -0.75,
          ASV_ID == "Zotu2" ~ y_value - 0.45,
          TRUE ~ y_value
        ),
        label_hjust = case_when(
          ASV_ID == "Zotu29" ~ 0.5,
          ASV_ID == "Zotu2" ~ 0.0,
          TRUE ~ 0.5
        ),
        label_vjust = case_when(
          ASV_ID == "Zotu29" ~ 1.0,
          ASV_ID == "Zotu2" ~ 1.0,
          TRUE ~ 1.0
        )
      )

    label_df <- label_df %>%
      filter(!ASV_ID %in% c("Zotu29", "Zotu2"))
  }

  x_breaks <- make_axis_breaks(x_axis_limit_percent)
  y_breaks <- make_axis_breaks(y_axis_limit)
  label_nudge_x <- x_axis_limit_percent * 0.03
  label_nudge_y <- y_axis_limit * 0.03

  label_df <- label_df %>%
    mutate(
      label_nudge_x = label_nudge_x,
      label_nudge_y = label_nudge_y
    )

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
      data = manual_label_df,
      aes(x = x_value, y = y_value, xend = label_x, yend = label_y),
      inherit.aes = FALSE,
      color = "grey35",
      linewidth = 0.35
    ) +
    geom_text(
      data = manual_label_df,
      aes(x = label_x, y = label_y, label = plot_label),
      inherit.aes = FALSE,
      family = base_family,
      size = 6.2,
      hjust = manual_label_df$label_hjust,
      vjust = manual_label_df$label_vjust
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = plot_label),
      size = 6.2,
      family = base_family,
      max.overlaps = Inf,
      seed = 42,
      nudge_x = label_df$label_nudge_x,
      nudge_y = label_df$label_nudge_y,
      box.padding = 0.6,
      point.padding = 0.5,
      force = 2.5,
      force_pull = 0.5,
      max.time = 3,
      max.iter = 20000,
      segment.color = "grey35",
      segment.size = 0.35,
      segment.alpha = 1,
      min.segment.length = 0
    ) +
    coord_cartesian(xlim = c(0, x_axis_limit_percent), ylim = c(-1, y_axis_limit), expand = TRUE) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(
      title = NULL,
      x = x_label,
      y = y_label
    ) +
    theme_classic(base_size = 16, base_family = base_family) +
    theme(
      axis.title = element_text(color = "black", size = 22),
      axis.text = element_text(color = "black", size = 22)
    )
}

summary_tbl <- read_csv(summary_path, show_col_types = FALSE)
y_axis_limit <- make_pretty_axis_limit(max(c(summary_tbl$mean_Cop_pro, summary_tbl$mean_Cop_ctr), na.rm = TRUE))

p_cop_pro <- make_comparison_plot(
  summary_tbl = summary_tbl,
  x_col = "mean_Prob",
  y_col = "mean_Cop_pro",
  x_label = "Mean relative abundance (%) in Prob",
  y_label = "Mean relative abundance (%) in Cop_pro",
  title_text = "D24: Prob vs Cop_pro",
  point_color = sample_type_colors[["Cop_pro_D24"]],
  y_axis_limit = y_axis_limit
)

p_cop_ctr <- make_comparison_plot(
  summary_tbl = summary_tbl,
  x_col = "mean_Prob",
  y_col = "mean_Cop_ctr",
  x_label = "Mean relative abundance (%) in Prob",
  y_label = "Mean relative abundance (%) in Cop_ctr",
  title_text = "D24: Prob vs Cop_ctr",
  point_color = sample_type_colors[["Cop_ctr_D24"]],
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

cat("Saved D24 ASV29-below plot:", output_path, "\n")
