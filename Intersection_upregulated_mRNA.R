##############################################
## 4-set Venn: Upregulated mRNAs
## GSE74530_Up, GSE30784_Up, GSE150469_Up, GSE246050_Up
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

genes1 <- get_genes("upregulated_ARRAYMRNA.csv")     # GSE74530
genes2 <- get_genes("upregulated_ARRAYMRNA_2.csv")   # GSE30784
genes3 <- get_genes("upregulated_RNASeqMRNA.csv")    # GSE150469
genes4 <- get_genes("upregulated_RNASeqMRNA_2.csv")  # GSE246050

## 2. Common genes across all 4 datasets ----------------------
common_genes_up <- Reduce(intersect, list(genes1, genes2, genes3, genes4))

cat("Number of common upregulated genes:", length(common_genes_up), "\n")
print(common_genes_up)

write.csv(
  data.frame(Gene = common_genes_up),
  "CommonGenes_up.csv",
  row.names = FALSE
)

## 3. Colours (same palette as down) --------------------------
fill_cols <- c("#F8766D", "#A3A500", "#00BFC4", "#C77CFF")

## 4. Venn -> TIFF 1200 dpi -----------------------------------
venn_up_tiff <- venn.diagram(
  x = list(
    GSE74530_Up  = genes1,
    GSE30784_Up  = genes2,
    GSE150469_Up = genes3,
    GSE246050_Up = genes4
  ),
  filename    = "Venn_up_4sets.tiff",
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
  
  # dataset names (same smaller size as down)
  cat.cex        = 0.7,
  cat.fontface   = "bold",
  cat.fontfamily = "sans",
  cat.col        = "black",
  
  margin = 0.07
)

## 5. Venn -> PNG 1200 dpi ------------------------------------
venn_up_png <- venn.diagram(
  x = list(
    GSE74530_Up  = genes1,
    GSE30784_Up  = genes2,
    GSE150469_Up = genes3,
    GSE246050_Up = genes4
  ),
  filename   = "Venn_up_4sets.png",
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
  
  cat.cex        = 0.7,
  cat.fontface   = "bold",
  cat.fontfamily = "sans",
  cat.col        = "black",
  
  margin     = 0.07
)
