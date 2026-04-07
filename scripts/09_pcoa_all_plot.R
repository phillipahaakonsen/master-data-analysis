# PCoA plot for all samples from PAST score and summary exports
# Input: PCoA_All_scores.csv and PCoA_All_summary.csv in project root
# Output: pcoa_all_plot.png in results directory

library(dplyr)
library(ggplot2)
library(grid)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
data_dir <- getOption("project_data_dir", getwd())
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

base_family <- "serif"
base_size <- 18
axis_title_size <- 22
axis_text_size <- 18
legend_text_size <- 18

group_order <- c(
  "Cop_pro_D16", "Cop_pro_D21", "Cop_pro_D24",
  "Cop_ctr_D16", "Cop_ctr_D21", "Cop_ctr_D24",
  "Alg_D16", "Alg_D21", "Alg_D24",
  "Prob_D16", "Prob_D21", "Prob_D24"
)

group_colors <- c(
  "Cop_pro_D16" = "#20B7EA",
  "Cop_ctr_D16" = "#FDC300",
  "Cop_pro_D21" = "#4169E1",
  "Cop_ctr_D21" = "#FF9500",
  "Cop_pro_D24" = "#1624F2",
  "Cop_ctr_D24" = "#FF6200",
  "Alg_D16" = "#FF8CCB",
  "Alg_D21" = "#FF78B8",
  "Alg_D24" = "#FF6FB3",
  "Prob_D16" = "#9ACD32",
  "Prob_D21" = "#86C400",
  "Prob_D24" = "#7FBF00"
)

pcoa_raw <- read.csv(
  file.path(data_dir, "PCoA_All_scores.csv"),
  sep = ";",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (ncol(pcoa_raw) >= 1) {
  names(pcoa_raw)[1] <- "sample_id"
}

pcoa_summary <- read.csv(
  file.path(data_dir, "PCoA_All_summary.csv"),
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
    sample_type = case_when(
      grepl("^Cop_pro_D16", sample_id) ~ "Cop_pro_D16",
      grepl("^Cop_ctr_D16", sample_id) ~ "Cop_ctr_D16",
      grepl("^Cop_pro_D21", sample_id) ~ "Cop_pro_D21",
      grepl("^Cop_ctr_D21", sample_id) ~ "Cop_ctr_D21",
      grepl("^Cop_pro_D24", sample_id) ~ "Cop_pro_D24",
      grepl("^Cop_ctr_D24", sample_id) ~ "Cop_ctr_D24",
      grepl("^Alg_D16", sample_id) ~ "Alg_D16",
      grepl("^Alg_D21", sample_id) ~ "Alg_D21",
      grepl("^Alg_D24", sample_id) ~ "Alg_D24",
      grepl("^Prob_D16", sample_id) ~ "Prob_D16",
      grepl("^Prob_D21", sample_id) ~ "Prob_D21",
      grepl("^Prob_D24", sample_id) ~ "Prob_D24",
      TRUE ~ NA_character_
    ),
    replicate_unit = case_when(
      grepl("^Cop_", sample_id) ~ sub("^Cop_[^_]+_[^_]+_([1-4])_.*$", "\\1", sample_id),
      TRUE ~ NA_character_
    ),
    day_group = case_when(
      grepl("D16", sample_id) ~ "D16",
      grepl("D21", sample_id) ~ "D21",
      grepl("D24", sample_id) ~ "D24",
      TRUE ~ NA_character_
    ),
    sample_class = case_when(
      grepl("^Cop_", sample_id) ~ "copepod",
      grepl("^Alg_", sample_id) ~ "microalgae",
      grepl("^Prob_", sample_id) ~ "probiotic_solution",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sample_type)) %>%
  mutate(sample_type = factor(sample_type, levels = group_order))

hull_plot <- pcoa_plot %>%
  group_by(sample_type) %>%
  filter(n() >= 3) %>%
  slice(chull(coord1, coord2)) %>%
  ungroup()

alg_hull_plot <- pcoa_plot %>%
  filter(sample_class == "microalgae") %>%
  slice(chull(coord1, coord2)) %>%
  mutate(hull_group = "Alg_all")

prob_hull_plot <- pcoa_plot %>%
  filter(sample_class == "probiotic_solution") %>%
  slice(chull(coord1, coord2)) %>%
  mutate(hull_group = "Prob_all")

copepod_plot <- pcoa_plot %>%
  filter(sample_class == "copepod") %>%
  mutate(replicate_unit = factor(replicate_unit, levels = c("1", "2", "3", "4")))

alg_plot <- pcoa_plot %>%
  filter(sample_class == "microalgae") %>%
  mutate(day_group = factor(day_group, levels = c("D16", "D21", "D24")))

prob_plot <- pcoa_plot %>%
  filter(sample_class == "probiotic_solution") %>%
  mutate(day_group = factor(day_group, levels = c("D16", "D21", "D24")))

x_limits <- c(floor(min(pcoa_plot$coord1) * 10) / 10, ceiling(max(pcoa_plot$coord1) * 10) / 10)
y_limits <- c(floor(min(pcoa_plot$coord2) * 10) / 10, ceiling(max(pcoa_plot$coord2) * 10) / 10)
x_breaks <- seq(x_limits[1], x_limits[2], by = 0.1)
y_breaks <- seq(y_limits[1], y_limits[2], by = 0.1)

p_pcoa_all <- ggplot() +
  geom_polygon(
    data = hull_plot,
    aes(x = coord1, y = coord2, group = sample_type, color = sample_type, fill = sample_type),
    alpha = 0.25,
    linewidth = 1,
    show.legend = FALSE
  ) +
  geom_polygon(
    data = alg_hull_plot,
    aes(x = coord1, y = coord2, group = hull_group),
    color = "#FF78B8",
    fill = "#FF78B8",
    alpha = 0.18,
    linewidth = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_polygon(
    data = prob_hull_plot,
    aes(x = coord1, y = coord2, group = hull_group),
    color = "#86C400",
    fill = "#86C400",
    alpha = 0.18,
    linewidth = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(
    data = copepod_plot,
    aes(x = coord1, y = coord2, color = sample_type, fill = sample_type, shape = replicate_unit, size = replicate_unit),
    stroke = 1.4,
    show.legend = FALSE
  ) +
  geom_point(
    data = alg_plot,
    aes(x = coord1, y = coord2, color = sample_type, shape = day_group),
    size = 5.2,
    stroke = 1.2,
    show.legend = FALSE
  ) +
  geom_point(
    data = prob_plot,
    aes(x = coord1, y = coord2, color = sample_type, shape = day_group),
    size = 5.6,
    stroke = 1.2,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = group_colors,
    breaks = group_order,
    labels = group_order
  ) +
  scale_fill_manual(
    values = group_colors,
    breaks = group_order,
    labels = group_order
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(
      "1" = 1,
      "2" = 8,
      "3" = 18,
      "4" = 0,
      "D16" = 2,
      "D21" = 17,
      "D24" = 4
    ),
    breaks = group_order,
    labels = group_order
  ) +
  scale_size_manual(
    name = NULL,
    values = c("1" = 4.9, "2" = 5.4, "3" = 5.4, "4" = 4.9),
    guide = "none"
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
    legend.position = "none"
  )

print(p_pcoa_all)

pcoa_all_file <- file.path(results_dir, "pcoa_all_plot.png")
png(filename = pcoa_all_file, width = 2400, height = 1500, res = 150)
grid.newpage()

legend_labels <- c(
  "Cop_pro_D16",
  "Cop_pro_D21",
  "Cop_pro_D24",
  "Cop_ctr_D16",
  "Cop_ctr_D21",
  "Cop_ctr_D24",
  "Alg",
  "Prob"
)

legend_colors <- c(
  "Cop_pro_D16" = unname(group_colors["Cop_pro_D16"]),
  "Cop_pro_D21" = unname(group_colors["Cop_pro_D21"]),
  "Cop_pro_D24" = unname(group_colors["Cop_pro_D24"]),
  "Cop_ctr_D16" = unname(group_colors["Cop_ctr_D16"]),
  "Cop_ctr_D21" = unname(group_colors["Cop_ctr_D21"]),
  "Cop_ctr_D24" = unname(group_colors["Cop_ctr_D24"]),
  "Alg" = "#FF78B8",
  "Prob" = "#86C400"
)

layout <- grid.layout(
  nrow = 1,
  ncol = 2,
  widths = unit(c(1.12, 0.34), "null")
)
pushViewport(viewport(layout = layout))

print(p_pcoa_all, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
for (i in seq_along(legend_labels)) {
  y_pos <- 0.88 - (i - 1) * 0.058
  key_x <- 0.13
  grid.rect(
    x = unit(key_x, "npc"),
    y = unit(y_pos, "npc"),
    width = unit(5.4, "mm"),
    height = unit(5.4, "mm"),
    gp = gpar(
      col = legend_colors[legend_labels[i]],
      fill = grDevices::adjustcolor(legend_colors[legend_labels[i]], alpha.f = 0.4),
      lwd = 1.2
    )
  )
  grid.text(
    label = legend_labels[i],
    x = unit(0.19, "npc"),
    y = unit(y_pos, "npc"),
    just = "left",
    gp = gpar(fontfamily = base_family, fontsize = legend_text_size, col = "black")
  )
}
popViewport()
popViewport()

dev.off()
cat("Saved plot:", pcoa_all_file, "\n")
