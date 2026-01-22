library(pheatmap)
library(grid)
library(ggplotify)
library(cowplot)
library(ggplot2)  # for ggsave()

gene_order <- c(
  "ADAM12", "ADORA2B", "BCAT1", "CXCL1", "EDIL3", "GLIPR1", "LOXL2", 
  "LY6K", "MMP13", "NRG1", "SLC25A32", "SLCO1B3", "TFRC", "BEX4", 
  "FAM221A", "GPX3", "MLPH", "RRAGD", "SVIP", "ZNF662", "ZNF677"
)

plot_and_save_heatmap <- function(expr_file, meta_file, dataset_name) {
  cat("Processing:", dataset_name, "\n")
  
  expr <- read.csv(expr_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- read.csv(meta_file, stringsAsFactors = FALSE)
  
  # Keep only two groups
  allowed_groups <- c("Normal", "Cancer")  # change to c("Normal", "Tumor") if needed
  meta <- meta[meta$Group %in% allowed_groups, , drop = FALSE]
  
  # Make sure gene names are consistent (upper case)
  rownames(expr) <- toupper(trimws(rownames(expr)))
  gene_order_upper <- toupper(gene_order)
  
  genes_present <- intersect(gene_order_upper, rownames(expr))
  if (length(genes_present) == 0) {
    stop(paste0("No matching genes in ", dataset_name))
  }
  
  # Subset & order genes
  expr_sub <- expr[match(genes_present, rownames(expr)), , drop = FALSE]
  expr_ordered <- expr_sub[order(match(rownames(expr_sub), gene_order_upper)), , drop = FALSE]
  
  # Z-score per gene
  expr_ordered <- t(scale(t(expr_ordered)))
  
  # Map back to original case for gene names
  gene_name_map <- setNames(gene_order, gene_order_upper)
  rownames(expr_ordered) <- gene_name_map[rownames(expr_ordered)]
  
  # Keep only common samples
  common_samples <- intersect(colnames(expr_ordered), meta$Sample)
  if (length(common_samples) == 0) {
    stop(paste0("No common samples in ", dataset_name, " after filtering groups"))
  }
  
  expr_ordered <- expr_ordered[, common_samples, drop = FALSE]
  meta_sub <- meta[match(common_samples, meta$Sample), , drop = FALSE]
  
  # Column annotation (Normal vs Cancer)
  annotation_col <- data.frame(Group = factor(meta_sub$Group))
  rownames(annotation_col) <- meta_sub$Sample
  
  # Row annotation: highlight ZNF662
  row_annotation <- data.frame(
    ZNF662_highlight = ifelse(rownames(expr_ordered) == "ZNF662", "ZNF662", "Other")
  )
  rownames(row_annotation) <- rownames(expr_ordered)
  
  # Annotation colors
  ann_colors <- list(
    Group = c("Normal" = "#4daf4a", "Cancer" = "#e41a1c"),
    ZNF662_highlight = c("ZNF662" = "#00BFFF", "Other" = "white")
  )
  
  # Create heatmap as pheatmap object
  p <- pheatmap(
    expr_ordered,
    show_colnames = TRUE,
    show_rownames = TRUE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_col = annotation_col,
    annotation_row = row_annotation,
    annotation_colors = ann_colors,
    main = dataset_name,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    border_color = NA,
    silent = TRUE
  )
  
  # Save individual PNG
  png(paste0(dataset_name, "_heatmap_ZNF662_highlight.png"),
      width = 8, height = 8, units = "in", res = 300)
  grid::grid.draw(p$gtable)
  dev.off()
  
  # Save individual PDF
  pdf(paste0(dataset_name, "_heatmap_ZNF662_highlight.pdf"),
      width = 8, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  dev.off()
  
  cat("✅ Saved heatmap for:", dataset_name, "\n")
  
  # Return as ggplot object for cowplot
  as.ggplot(p$gtable)
}

# Dataset list -------------------------------------------------
datasets <- list(
  GSE74530  = list(expr = "GSE74530_expression.csv",  meta = "GSE74530_meta.csv"),
  GSE150469 = list(expr = "GSE150469_expression.csv", meta = "GSE150469_meta.csv"),
  GSE30784  = list(expr = "GSE30784_expression.csv",  meta = "GSE30784_meta.csv"),
  GSE246050 = list(expr = "GSE246050_expression.csv", meta = "GSE246050_meta.csv")
)

plot_list  <- list()
label_list <- c("a", "b", "c", "d")  # lowercase panel labels

i <- 1
for (ds in names(datasets)) {
  tryCatch({
    gg <- plot_and_save_heatmap(datasets[[ds]]$expr, datasets[[ds]]$meta, ds)
    plot_list[[i]] <- gg
    i <- i + 1
  }, error = function(e) {
    message(paste0("❌ Error in ", ds, ": ", e$message))
  })
}

# Combined 2×2 panel with labels a, b, c, d
combined_plot <- cowplot::plot_grid(
  plotlist   = plot_list,
  labels     = label_list,
  label_size = 16,
  ncol       = 2
)

# Save combined as high-res PNG
ggsave(
  "Combined_Heatmaps_4in1.png",
  plot   = combined_plot,
  width  = 14,
  height = 14,
  dpi    = 600
)

# Save combined as high-res TIFF using base device (avoids TIFFOpen error)
tiff(
  "Combined_Heatmaps_4in1.tiff",
  width  = 14,
  height = 14,
  units  = "in",
  res    = 600,
  compression = "lzw"
)
print(combined_plot)
dev.off()

# Optional: PDF too
ggsave(
  "Combined_Heatmaps_4in1.pdf",
  plot   = combined_plot,
  width  = 14,
  height = 14
)

cat("✅ Saved combined 2x2 heatmap (PNG + TIFF + PDF)\n")
