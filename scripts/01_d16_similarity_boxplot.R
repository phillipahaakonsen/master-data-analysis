# D16 similarity boxplot
# Input: BrayCurtisSimilarities_copepods.csv in project root
# Output: d16_similarity_boxplot.png in results directory

source("scripts/01_similarity_boxplot_helper.R", local = TRUE)
make_similarity_boxplot(
  "D16",
  "d16_similarity_boxplot.png",
  c("#D6F0FC", "#FDECBF", "#D9B3FF")
)
