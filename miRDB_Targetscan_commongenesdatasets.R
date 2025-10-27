# Load libraries
library(VennDiagram)
library(readr)
library(dplyr)

# Read your CSV files
dataset1 <- read_csv("miRDB.csv")
dataset2 <- read_csv("Targetscan.csv")
dataset3 <- read_csv("DEmRNA_common.csv")


# Clean gene lists: remove NA and empty strings
genes1 <- dataset1$Gene %>% na.omit() %>% unique() %>% .[. != ""]
genes2 <- dataset2$Gene %>% na.omit() %>% unique() %>% .[. != ""]
genes3 <- dataset3$Gene %>% na.omit() %>% unique() %>% .[. != ""]


# Find common genes across all 4 datasets
common_genes <- Reduce(intersect, list(genes1, genes2, genes3))

# Print them
cat("Number of common genes:", length(common_genes), "\n")
print(common_genes)

write.csv(data.frame(Gene = common_genes), "common_genes_databases.csv", row.names = FALSE)

# Now safe to plot
venn.plot <- venn.diagram(
  x = list(miRDB
     = genes1,
    TargetScan = genes2,
    CommonGenesDatasets= genes3
  ),
  filename = "VennDiagram_3.png",
  imagetype = "png",
  height = 4000,
  width = 4000,
  resolution = 300,
  col = "transparent",
  fill = c("blue", "green", "orange"),
  alpha = 0.5,
  cex =3,
  cat.cex = 2.5,
  cat.col = c( "blue", "green", "orange")
)
