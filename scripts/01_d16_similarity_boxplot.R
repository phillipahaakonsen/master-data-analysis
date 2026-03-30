# D16 similarity boxplot
# Input: D16.csv in project root
# Output: d16_similarity_boxplot.png in results directory

results_dir <- getOption("project_results_dir", file.path(getwd(), "results"))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv("D16.csv")

# Rename group labels for plotting
label_map <- c(
  "within_pro" = "Probiotic-treated",
  "within_ctr" = "Untreated control",
  "between" = "Probiotic-treated vs. Untreated control"
)

data$group_key <- tolower(data$group)
data$group_label <- unname(label_map[data$group_key])

if (any(is.na(data$group_label))) {
  missing_levels <- unique(data$group[is.na(data$group_label)])
  stop(sprintf("Unmapped group labels found: %s", paste(missing_levels, collapse = ", ")))
}

data$group_label <- factor(
  data$group_label,
  levels = c(
    "Probiotic-treated",
    "Untreated control",
    "Probiotic-treated vs. Untreated control"
  )
)

output_file <- file.path(results_dir, "d16_similarity_boxplot.png")
png(filename = output_file, width = 1000, height = 700, res = 120)

boxplot(
  similarity ~ group_label,
  data = data,
  main = "Similarity Scores by Group",
  xlab = "Group",
  ylab = "Similarity",
  col = c("lightblue", "lightgreen", "lightpink")
)

dev.off()
cat("Saved boxplot to:", output_file, "\n")
