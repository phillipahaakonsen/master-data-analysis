# D24 similarity boxplot
# Input: BrayCurtisSimilarities_copepods.csv in project root
# Output: d24_similarity_boxplot.png in results directory

source("scripts/01_similarity_boxplot_helper.R", local = TRUE)
make_similarity_boxplot(
  "D24",
  "d24_similarity_boxplot.png",
  c("#5D8FE8", "#EC9A43", "#D9B3FF")
)
