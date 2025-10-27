# Load libraries
library(VennDiagram)
library(readr)
library(dplyr)

# Read your CSV files
dataset1 <- read_csv("DElncRNA.csv")
dataset2 <- read_csv("Starbase.csv")



# Clean gene lists: remove NA and empty strings
genes1 <- dataset1$Symbol %>% na.omit() %>% unique() %>% .[. != ""]
genes2 <- dataset2$Symbol%>% na.omit() %>% unique() %>% .[. != ""]


# Find common genes across all 4 datasets
common_genes <- Reduce(intersect, list(genes1, genes2))

# Print them
cat("Number of common genes:", length(common_genes), "\n")
print(common_genes)

write.csv(data.frame(Gene = common_genes), "common_lncRNA.csv", row.names = FALSE)

# Now safe to plot
venn.plot <- venn.diagram(
  x = list(DElncRNA
           = genes1,
           starBase = genes2
           ),
  filename = "VennDiagram_4.png",
  imagetype = "png",
  height = 6000,
  width = 6000,
  resolution = 300,
  col = "transparent",
  fill = c("blue", "orange"),
  alpha = 0.5,
  cex =3,
  cat.cex = 2.5,
  cat.col = c( "blue", "orange")
)
