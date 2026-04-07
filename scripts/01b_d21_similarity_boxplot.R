# D21 similarity boxplot
# Input: BrayCurtisSimilarities_copepods.csv in project root
# Output: d21_similarity_boxplot.png in results directory

source("scripts/01_similarity_boxplot_helper.R", local = TRUE)
make_similarity_boxplot(
  "D21",
  "d21_similarity_boxplot.png",
  c("#8CCEF3", "#F6CB74", "#D9B3FF")
)
