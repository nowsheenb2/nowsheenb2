# Load libraries
library(VennDiagram)
library(readr)
library(dplyr)

# Read your CSV files
dataset1 <- read_csv("upregulated_ARRAYMRNA.csv")
dataset2 <- read_csv("upregulated_ARRAYMRNA_2.csv")
dataset3 <- read_csv("upregulated_RNASeqMRNA.csv")
dataset4 <- read_csv("upregulated_RNASeqMRNA_2.csv")


# Clean gene lists: remove NA and empty strings
genes1 <- dataset1$Gene %>% na.omit() %>% unique() %>% .[. != ""]
genes2 <- dataset2$Gene %>% na.omit() %>% unique() %>% .[. != ""]
genes3 <- dataset3$Gene %>% na.omit() %>% unique() %>% .[. != ""]
genes4 <- dataset4$Gene %>% na.omit() %>% unique() %>% .[. != ""]

# Find common genes across all 4 datasets
common_genes <- Reduce(intersect, list(genes1, genes2, genes3, genes4))

# Print them
cat("Number of common genes:", length(common_genes), "\n")
print(common_genes)

write.csv(data.frame(Gene = common_genes), "CommonGenes_up.csv", row.names = FALSE)

# Now safe to plot
venn.plot <- venn.diagram(
  x = list(
    GSE74530_UP = genes1,
    GSE30784_UP = genes2,
    GSE150469_UP = genes3,
    GSE246050_UP = genes4
  ),
  filename = "VennDiagram_2.png",
  imagetype = "png",
  height = 4500,
  width = 6000,
  resolution = 300,
  col = "transparent",
  fill = c("red", "blue", "green", "orange"),
  alpha = 0.5,
  cex = 3,
  cat.cex = 2.5,
  cat.col = c("red", "blue", "green", "orange")
)
