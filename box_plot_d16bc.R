# 1. Import the CSV file
data <- read.csv("D16.csv")

# 1b. Rename group labels for plotting
data$group_label <- tolower(data$group)
data$group_label[data$group_label == "within_pro"] <- "Probiotic-treated"
data$group_label[data$group_label == "within_ctr"] <- "Untreated control"
data$group_label[data$group_label == "between"] <- "Probiotic-treated vs. Untreated control"
data$group_label <- factor(
  data$group_label,
  levels = c(
    "Probiotic-treated",
    "Untreated control",
    "Probiotic-treated vs. Untreated control"
  )
)

# 2. Save the box plot in the project root folder
output_file <- file.path(getwd(), "boxplot.png")
png(filename = output_file, width = 1000, height = 700, res = 120)

# The '~' tells R to plot 'similarity' broken down by 'group'
boxplot(similarity ~ group_label, data = data,
        main = "Similarity Scores by Group",  # Title of the plot
        xlab = "Group",                       # X-axis label
        ylab = "Similarity",                  # Y-axis label
        col = c("lightblue", "lightgreen", "lightpink")) # Box colors

dev.off()
cat("Saved boxplot to:", output_file, "\n")
