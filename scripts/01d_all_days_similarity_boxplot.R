# Combined mean-similarity boxplot across D16, D21, and D24
# Input: BrayCurtisSimilarities_copepods.csv in project root
# Output: all_days_similarity_boxplot.png in results directory

library(dplyr)
library(ggplot2)

source("scripts/01_similarity_boxplot_helper.R", local = TRUE)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

base_family <- "serif"
base_size <- 18
axis_title_size <- 20
axis_text_size <- 15
legend_text_size <- 15
axis_line_width <- 0.45
axis_tick_width <- 0.45
x_axis_text_size <- 20

similarity_matrix <- read_similarity_matrix()
sample_meta <- build_similarity_sample_meta(similarity_matrix)

similarity_all <- bind_rows(
  build_day_mean_similarity_data("D16", similarity_matrix = similarity_matrix, sample_meta = sample_meta),
  build_day_mean_similarity_data("D21", similarity_matrix = similarity_matrix, sample_meta = sample_meta),
  build_day_mean_similarity_data("D24", similarity_matrix = similarity_matrix, sample_meta = sample_meta)
) %>%
  mutate(
    group_key = tolower(group),
    group_label = recode(
      group_key,
      "within_pro" = "Cop_pro",
      "within_ctr" = "Cop_ctr",
      "between" = "Cop_pro vs. Cop_ctr"
    ),
    group_label = factor(
      group_label,
      levels = c(
        "Cop_pro",
        "Cop_ctr",
        "Cop_pro vs. Cop_ctr"
      )
    ),
    day = factor(day, levels = c("D16", "D21", "D24"))
  )

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

p_similarity_all <- ggplot(
  similarity_all,
  aes(x = group_label, y = similarity, fill = day, color = day)
) +
  geom_boxplot(
    position = position_dodge(width = 0.72),
    width = 0.62,
    linewidth = 0.9,
    alpha = 0.45,
    staplewidth = 0.5
  ) +
  scale_fill_manual(values = day_colors) +
  scale_color_manual(values = outline_colors, guide = "none") +
  scale_x_discrete(
    labels = c(
      "Cop_pro" = "Cop_pro",
      "Cop_ctr" = "Cop_ctr",
      "Cop_pro vs. Cop_ctr" = "Cop_pro vs. Cop_ctr"
    )
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Mean Bray-Curtis similarity",
    fill = NULL
  ) +
  theme_classic(base_size = base_size, base_family = base_family) +
  theme(
    axis.title = element_text(size = axis_title_size, color = "black", face = "plain"),
    axis.text = element_text(size = axis_text_size, color = "black"),
    axis.text.x = element_text(size = x_axis_text_size, color = "black", lineheight = 0.95),
    axis.line = element_line(linewidth = axis_line_width, color = "grey30"),
    axis.ticks = element_line(linewidth = axis_tick_width, color = "grey30"),
    legend.text = element_text(size = legend_text_size),
    legend.position = "top",
    legend.direction = "horizontal"
  )

print(p_similarity_all)

output_file <- file.path(results_dir, "all_days_similarity_boxplot.png")
ggsave(
  filename = output_file,
  plot = p_similarity_all,
  width = 11,
  height = 8,
  dpi = 300
)
cat("Saved plot:", output_file, "\n")
