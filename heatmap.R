library(pheatmap)
library(grid)
library(ggplotify)
library(cowplot)
library(ggplot2)  # for ggsave()

gene_order <- c("ADAM12", "ADORA2B", "BCAT1", "CXCL1", "EDIL3", "GLIPR1", "LOXL2", 
                "LY6K", "MMP13", "NRG1", "SLC25A32", "SLCO1B3", "TFRC", "BEX4", 
                "FAM221A", "GPX3", "MLPH", "RRAGD", "SVIP", "ZNF662", "ZNF677")

plot_and_save_heatmap <- function(expr_file, meta_file, dataset_name) {
  cat("Processing:", dataset_name, "\n")
  expr <- read.csv(expr_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- read.csv(meta_file, stringsAsFactors = FALSE)
  
  # Filter metadata for only two groups - modify as needed
  allowed_groups <- c("Normal", "Cancer")  # change to c("Normal", "Tumor") if desired
  meta <- meta[meta$Group %in% allowed_groups, , drop = FALSE]
  
  rownames(expr) <- toupper(trimws(rownames(expr)))
  gene_order_upper <- toupper(gene_order)
  genes_present <- intersect(gene_order_upper, rownames(expr))
  if (length(genes_present) == 0) stop(paste0("No matching genes in ", dataset_name))
  
  expr_sub <- expr[match(genes_present, rownames(expr)), , drop = FALSE]
  expr_ordered <- expr_sub[order(match(rownames(expr_sub), gene_order_upper)), , drop = FALSE]
  expr_ordered <- t(scale(t(expr_ordered)))  # z-score normalization
  
  gene_name_map <- setNames(gene_order, gene_order_upper)
  rownames(expr_ordered) <- gene_name_map[rownames(expr_ordered)]
  
  # Keep only samples present in filtered metadata
  common_samples <- intersect(colnames(expr_ordered), meta$Sample)
  if (length(common_samples) == 0) stop(paste0("No common samples in ", dataset_name, " after filtering groups"))
  
  expr_ordered <- expr_ordered[, common_samples, drop = FALSE]
  meta_sub <- meta[match(common_samples, meta$Sample), , drop = FALSE]
  
  annotation_col <- data.frame(Group = factor(meta_sub$Group))
  rownames(annotation_col) <- meta_sub$Sample
  ann_colors <- list(Group = c("Normal" = "#4daf4a", "Cancer" = "#e41a1c"))  # Update colors to match two groups
  

  row_annotation <- data.frame(ZNF662_highlight = ifelse(rownames(expr_ordered) == "ZNF662", "ZNF662", "Other"))
  rownames(row_annotation) <- rownames(expr_ordered)
  row_ann_colors <- list(ZNF662_highlight = c(ZNF662 = "#00BFFF", Other = "white"))
  
  
  p <- pheatmap(expr_ordered,
                show_colnames = TRUE,
                show_rownames = TRUE,
                cluster_cols = FALSE,
                cluster_rows = FALSE,
                annotation_col = annotation_col,
                annotation_colors = c(ann_colors, row_ann_colors),
                annotation_row = row_annotation,
                main = dataset_name,
                fontsize = 10,
                fontsize_row = 8,
                fontsize_col = 8,
                border_color = NA,
                silent = TRUE)
  
 
  png(paste0(dataset_name, "_heatmap_ZNF662_highlight.png"), width = 1200, height = 1200, res = 150)
  grid::grid.draw(p$gtable)
  dev.off()
  
  pdf(paste0(dataset_name, "_heatmap_ZNF662_highlight.pdf"), width = 10, height = 10)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  dev.off()
  
  
  cat("✅ Saved heatmap for:", dataset_name, "\n")
  return(as.ggplot(p$gtable))
}

datasets <- list(
  GSE74530 = list(expr = "GSE74530_expression.csv", meta = "GSE74530_meta.csv"),
  GSE150469 = list(expr = "GSE150469_expression.csv", meta = "GSE150469_meta.csv"),
  GSE30784 = list(expr = "GSE30784_expression.csv", meta = "GSE30784_meta.csv"),
  GSE246050 = list(expr = "GSE246050_expression.csv", meta = "GSE246050_meta.csv")
)

plot_list <- list()
label_list <- c("A", "B", "C", "D")

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

combined_plot <- cowplot::plot_grid(
  plotlist = plot_list,
  labels = paste0(label_list, ")"),
  label_size = 16,
  ncol = 2
)

ggsave("Combined_Heatmaps_4in1.png", plot = combined_plot, width = 14, height = 14, dpi = 300)
ggsave("Combined_Heatmaps_4in1.pdf", plot = combined_plot, width = 14, height = 14)

cat("✅ Saved combined 2x2 heatmap with only two groups (PNG + PDF)\n")


