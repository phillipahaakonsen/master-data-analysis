# PCoA plot for all copepod samples from PAST score and summary exports
# Input: PCoA_AllCopepods_Scores.csv and PCoA_AllCopepods_summary.csv in project root
# Output: pcoa_all_copepods_plot.png in results directory

library(dplyr)
library(ggplot2)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
data_dir <- getOption("project_data_dir", getwd())
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

base_family <- "serif"
base_size <- 18
axis_title_size <- 22
axis_text_size <- 18
legend_text_size <- 18

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

treatment_shapes <- c(
  "Cop_pro" = 21,
  "Cop_ctr" = 24
)

pcoa_raw <- read.csv(
  file.path(data_dir, "PCoA_AllCopepods_Scores.csv"),
  sep = ";",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)
if (ncol(pcoa_raw) >= 1) {
  names(pcoa_raw)[1] <- "sample_id"
}

pcoa_summary <- read.csv(
  file.path(data_dir, "PCoA_AllCopepods_summary.csv"),
  sep = ";",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
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

coord1_scale <- if (length(coord1_eigenvalue) == 1 && coord1_eigenvalue > 0) sqrt(coord1_eigenvalue) else 1
coord2_scale <- if (length(coord2_eigenvalue) == 1 && coord2_eigenvalue > 0) sqrt(coord2_eigenvalue) else 1

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

pcoa_plot <- pcoa_raw %>%
  transmute(
    sample_id = sample_id,
    coord1 = as.numeric(`Coord 1`) * coord1_scale,
    coord2 = as.numeric(`Coord 2`) * coord2_scale,
    treatment = case_when(
      grepl("^Cop_pro_", sample_id) ~ "Cop_pro",
      grepl("^Cop_ctr_", sample_id) ~ "Cop_ctr",
      TRUE ~ NA_character_
    ),
    day = sub("^(Cop_(?:pro|ctr))_(D[0-9]+)_.+$", "\\2", sample_id, perl = TRUE)
  ) %>%
  filter(!is.na(treatment)) %>%
  mutate(
    treatment = factor(treatment, levels = c("Cop_pro", "Cop_ctr")),
    day = factor(day, levels = c("D16", "D21", "D24"))
  )

if (any(is.na(pcoa_plot$day))) {
  stop("Could not parse one or more sampling days from PCoA_AllCopepods_Scores.csv.")
}

x_limits <- c(floor(min(pcoa_plot$coord1) * 10) / 10, ceiling(max(pcoa_plot$coord1) * 10) / 10)
y_limits <- c(floor(min(pcoa_plot$coord2) * 10) / 10, ceiling(max(pcoa_plot$coord2) * 10) / 10)
x_breaks <- pretty(x_limits, n = 7)
y_breaks <- pretty(y_limits, n = 7)
y_step <- if (length(y_breaks) > 1) diff(y_breaks)[1] else 0.1
x_limits <- c(min(x_breaks), 0.60)
y_limits <- c(min(y_breaks), max(y_breaks) + y_step)
x_breaks <- pretty(x_limits, n = 7)
y_breaks <- pretty(y_limits, n = 7)
x_breaks <- sort(unique(c(x_breaks[x_breaks <= x_limits[2]], x_limits[2])))
axis_cap_x <- 0.025 * diff(x_limits)
axis_cap_y <- 0.025 * diff(y_limits)

p_pcoa_all_copepods <- ggplot(
  pcoa_plot,
  aes(x = coord1, y = coord2, fill = day, color = day, shape = treatment)
) +
  geom_point(
    size = 6.6,
    stroke = 1.2
  ) +
  annotate(
    "segment",
    x = x_limits[1],
    xend = x_limits[1] + axis_cap_x,
    y = y_limits[2],
    yend = y_limits[2],
    linewidth = 0.45,
    color = "grey30"
  ) +
  annotate(
    "segment",
    x = x_limits[2],
    xend = x_limits[2],
    y = y_limits[1],
    yend = y_limits[1] + axis_cap_y,
    linewidth = 0.45,
    color = "grey30"
  ) +
  scale_fill_manual(
    values = day_colors,
    breaks = c("D16", "D21", "D24"),
    labels = c("D16", "D21", "D24"),
    name = NULL
  ) +
  scale_color_manual(
    values = outline_colors,
    guide = "none"
  ) +
  scale_shape_manual(
    values = treatment_shapes,
    breaks = c("Cop_pro", "Cop_ctr"),
    labels = c("Cop_pro", "Cop_ctr"),
    name = NULL
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 22,
        size = 5.0,
        stroke = 1.0,
        color = unname(outline_colors[c("D16", "D21", "D24")])
      )
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(fill = "white", color = "black", size = 5.0, stroke = 1.2)
    )
  ) +
  scale_x_continuous(
    limits = x_limits,
    breaks = x_breaks,
    labels = scales::label_number(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    labels = scales::label_number(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  labs(
    x = coord1_label,
    y = coord2_label
  ) +
  theme_classic(base_size = base_size, base_family = base_family) +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size, color = "black"),
    legend.text = element_text(size = legend_text_size),
    legend.position = "right",
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.height = grid::unit(0.85, "cm"),
    legend.key.width = grid::unit(0.85, "cm"),
    legend.spacing.y = grid::unit(0.15, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(10, 14, 10, 10)
  )

output_file <- file.path(results_dir, "pcoa_all_copepods_plot.png")
ggsave(
  filename = output_file,
  plot = p_pcoa_all_copepods,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Saved plot:", output_file, "\n")
