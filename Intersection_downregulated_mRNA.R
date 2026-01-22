##############################################
## 4-set Venn: Downregulated mRNAs
## GSE74530_Down, GSE30784_Down, GSE150469_Down, GSE246050_Down
##############################################

library(VennDiagram)
library(readr)
library(dplyr)

## 1. Read and clean gene lists -------------------------------
get_genes <- function(file, gene_col = "Gene") {
  df <- read_csv(file, show_col_types = FALSE)
  df[[gene_col]] %>%
    na.omit() %>%
    unique() %>%
    .[. != ""]
}

genes1 <- get_genes("downregulated_ARRAYMRNA.csv")     # GSE74530
genes2 <- get_genes("downregulated_ARRAYMRNA_2.csv")   # GSE30784
genes3 <- get_genes("downregulated_RNASeqMRNA.csv")    # GSE150469
genes4 <- get_genes("downregulated_RNASeqMRNA_2.csv")  # GSE246050

## 2. Common genes across all 4 datasets ----------------------
common_genes <- Reduce(intersect, list(genes1, genes2, genes3, genes4))

cat("Number of common downregulated genes:", length(common_genes), "\n")
print(common_genes)

write.csv(
  data.frame(Gene = common_genes),
  "CommonGenes_down.csv",
  row.names = FALSE
)

## 3. Colours --------------------------------------------------
fill_cols <- c("#F8766D", "#A3A500", "#00BFC4", "#C77CFF")

## 4. Venn -> TIFF 1200 dpi -----------------------------------
venn_down_tiff <- venn.diagram(
  x = list(
    GSE74530_Down  = genes1,
    GSE30784_Down  = genes2,
    GSE150469_Down = genes3,
    GSE246050_Down = genes4
  ),
  filename    = "Venn_down_4sets.tiff",
  imagetype   = "tiff",
  height      = 4800,      # 4 in @ 1200 dpi
  width       = 4800,
  resolution  = 1200,
  compression = "lzw",
  
  col   = "#555555",
  lwd   = 0.5,
  lty   = "solid",
  fill  = fill_cols,
  alpha = 0.45,
  
  # numbers inside regions
  cex        = 0.8,
  fontface   = "plain",
  fontfamily = "sans",
  
  # dataset names (a little smaller)
  cat.cex        = 0.7,   # << smaller than before
  cat.fontface   = "bold",
  cat.fontfamily = "sans",
  cat.col        = "black",
  
  margin = 0.07
)

## 5. Venn -> PNG 1200 dpi ------------------------------------
venn_down_png <- venn.diagram(
  x = list(
    GSE74530_Down  = genes1,
    GSE30784_Down  = genes2,
    GSE150469_Down = genes3,
    GSE246050_Down = genes4
  ),
  filename   = "Venn_down_4sets.png",
  imagetype  = "png",
  height     = 4800,
  width      = 4800,
  resolution = 1200,
  
  col        = "#555555",
  lwd        = 0.5,
  lty        = "solid",
  fill       = fill_cols,
  alpha      = 0.45,
  
  cex        = 0.8,
  fontface   = "plain",
  fontfamily = "sans",
  
  cat.cex        = 0.7,   # << smaller labels
  cat.fontface   = "bold",
  cat.fontfamily = "sans",
  cat.col        = "black",
  
  margin     = 0.07
)
