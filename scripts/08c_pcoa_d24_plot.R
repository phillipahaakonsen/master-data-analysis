# PCoA D24 plot from PAST scores export
# Input: PCoA_D24_scores.csv and PCoA_D24_summary.csv in project root
# Output: pcoa_d24_plot.png in results directory

library(ggplot2)
library(grid)
source("scripts/08_pcoa_plot_helper.R", local = TRUE)

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
data_dir <- getOption("project_data_dir", getwd())
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

base_family <- "serif"
legend_text_size <- 18

p_d24 <- build_pcoa_plot(
  file.path(data_dir, "PCoA_D24_scores.csv"),
  file.path(data_dir, "PCoA_D24_summary.csv"),
  "D24",
  probiotic_color = "#1624F2",
  control_color = "#FF6200"
) +
  theme(legend.position = "none")

legend_order <- c(
  "Cop_pro_D16", "Cop_pro_D21", "Cop_pro_D24",
  "Cop_ctr_D16", "Cop_ctr_D21", "Cop_ctr_D24"
)

legend_colors <- c(
  "Cop_pro_D16" = "#20B7EA",
  "Cop_pro_D21" = "#4169E1",
  "Cop_pro_D24" = "#1624F2",
  "Cop_ctr_D16" = "#FDC300",
  "Cop_ctr_D21" = "#FF9500",
  "Cop_ctr_D24" = "#FF6200"
)

output_file <- file.path(results_dir, "pcoa_d24_plot.png")
png(filename = output_file, width = 3900, height = 2400, res = 300)
grid.newpage()

layout <- grid.layout(
  nrow = 1,
  ncol = 2,
  widths = unit(c(1.12, 0.34), "null")
)
pushViewport(viewport(layout = layout))

print(p_d24, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
for (i in seq_along(legend_order)) {
  y_pos <- 0.84 - (i - 1) * 0.072
  key_x <- 0.13
  grid.rect(
    x = unit(key_x, "npc"),
    y = unit(y_pos, "npc"),
    width = unit(5.4, "mm"),
    height = unit(5.4, "mm"),
    gp = gpar(
      col = legend_colors[legend_order[i]],
      fill = grDevices::adjustcolor(legend_colors[legend_order[i]], alpha.f = 0.4),
      lwd = 1.2
    )
  )
  grid.text(
    label = legend_order[i],
    x = unit(0.19, "npc"),
    y = unit(y_pos, "npc"),
    just = "left",
    gp = gpar(fontfamily = base_family, fontsize = legend_text_size, col = "black")
  )
}
popViewport()

dev.off()
cat("Saved plot:", output_file, "\n")
