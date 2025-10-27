# --- ROC utilities for gene/lncRNA/miRNA markers ------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(ggplot2)
  library(purrr)
  library(stringr)
})

# Core: compute oriented ROC (AUC >= 0.5), return roc, auc, ci
.compute_oriented_roc <- function(response_factor, predictor_numeric, ci_method = c("delong","bootstrap"),
                                  boot_n = 2000, boot_stratified = TRUE, quiet = TRUE) {
  ci_method <- match.arg(ci_method)
  # Original
  roc_orig <- roc(response_factor, predictor_numeric, levels = levels(response_factor), quiet = quiet)
  auc_orig <- as.numeric(auc(roc_orig))
  # Flipped
  roc_flip <- roc(response_factor, -predictor_numeric, levels = levels(response_factor), quiet = quiet)
  auc_flip <- as.numeric(auc(roc_flip))
  # Choose orientation with higher AUC (=> >= 0.5)
  if (auc_orig >= auc_flip) {
    roc_obj <- roc_orig; auc_val <- auc_orig
  } else {
    roc_obj <- roc_flip; auc_val <- auc_flip
  }
  # CI
  if (ci_method == "delong") {
    ci_val <- ci.auc(roc_obj, method = "delong")
  } else {
    # bootstrap CI as an option (slower, set seed upstream if needed)
    ci_val <- ci.auc(roc_obj, method = "bootstrap", boot.n = boot_n, boot.stratified = boot_stratified)
  }
  list(roc = roc_obj,
       auc = auc_val,
       ci_low = as.numeric(ci_val[1]),
       ci_high = as.numeric(ci_val[3]))
}

# Plot helper (adds diagonal and nice theming)
.plot_roc <- function(roc_obj, title, subtitle_color = "#2c3e50") {
  ggroc(roc_obj, size = 1.2, color = subtitle_color) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", linewidth = 0.6) +
    labs(title = title,
         x = "1 - Specificity",
         y = "Sensitivity") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

# Save with multiple formats
.save_plot_multi <- function(plot, filepath_no_ext, formats = c("tiff"), width = 7, height = 6, dpi = 300) {
  formats <- tolower(formats)
  for (fmt in formats) {
    device <- switch(fmt,
                     tiff = "tiff",
                     tif  = "tiff",
                     pdf  = "pdf",
                     svg  = "svg",
                     png  = "png",
                     stop("Unsupported filetype: ", fmt, ". Use tiff/pdf/svg/png.")
    )
    ggsave(paste0(filepath_no_ext, ".", fmt), plot = plot,
           width = width, height = height, units = "in", dpi = dpi, device = device)
  }
}

# General driver: one function for any marker list
# expr_file: CSV with columns: Gene, <Sample1>, <Sample2>, ...
# meta_file: CSV with columns: Sample, Group (values exactly 'Normal' or 'Cancer')
# markers: character vector of gene/lncRNA/miRNA names to evaluate
# dataset_name: used in plot titles and filenames
evaluate_markers_roc <- function(expr_file, meta_file, markers,
                                 dataset_name = NULL,
                                 save_dir = ".",
                                 formats = c("tiff"),
                                 ci_method = "delong",
                                 boot_n = 2000) {
  
  expr <- read_csv(expr_file, show_col_types = FALSE)
  meta <- read_csv(meta_file, show_col_types = FALSE)
  
  # Basic checks
  stopifnot(all(c("Gene") %in% names(expr)))
  stopifnot(all(c("Sample","Group") %in% names(meta)))
  
  # Factor order: controls (Normal) first, cases (Cancer) second
  meta <- meta %>%
    mutate(Group = factor(Group, levels = c("Normal","Cancer"))) %>%
    distinct(Sample, .keep_all = TRUE) # de-dup samples if any
  
  # Long form expression
  expr_long <- expr %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")
  
  # Optional: drop NAs
  expr_long <- expr_long %>% filter(!is.na(Expression))
  
  # Prepare output dir
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # Iterate markers
  results <- map_dfr(markers, function(mk) {
    x <- expr_long %>% filter(Gene == mk)
    if (nrow(x) == 0) {
      message("Marker not found: ", mk, " in ", expr_file, " — skipping")
      return(tibble(Marker = mk, AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                    N = NA_integer_, Dataset = dataset_name %||% basename(expr_file),
                    Status = "not_found"))
    }
    merged <- x %>% inner_join(meta, by = "Sample")
    
    # Check group coverage and counts
    if (!all(levels(meta$Group) %in% merged$Group)) {
      message("Marker ", mk, ": missing one of the groups after merge — skipping")
      return(tibble(Marker = mk, AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                    N = nrow(merged), Dataset = dataset_name %||% basename(expr_file),
                    Status = "missing_group"))
    }
    # Compute ROC
    roc_out <- .compute_oriented_roc(merged$Group, merged$Expression,
                                     ci_method = ci_method, boot_n = boot_n, quiet = TRUE)
    
    # Plot & save
    auc_str <- sprintf("AUC = %.3f (95%% CI: %.3f–%.3f)", roc_out$auc, roc_out$ci_low, roc_out$ci_high)
    ttl <- if (is.null(dataset_name)) {
      paste0("ROC: ", mk)
    } else {
      paste0("ROC: ", mk, " (", dataset_name, ")")
    }
    p <- .plot_roc(roc_out$roc, title = ttl) +
      labs(subtitle = auc_str)
    
    safe_marker <- gsub("[^A-Za-z0-9]+", "_", mk)
    safe_ds <- gsub("[^A-Za-z0-9]+", "_", dataset_name %||% "dataset")
    filepath_no_ext <- file.path(save_dir, paste0("ROC_", safe_marker, "_", safe_ds))
    .save_plot_multi(p, filepath_no_ext, formats = formats)
    
    tibble(Marker = mk,
           AUC = roc_out$auc,
           CI_low = roc_out$ci_low,
           CI_high = roc_out$ci_high,
           N = nrow(merged),
           Dataset = dataset_name %||% basename(expr_file),
           Status = "ok",
           PlotBase = filepath_no_ext)
  })
  
  results
}

# --- Example usage mirroring your current runs -------------------------------

# mRNA (ZNF662 across multiple datasets)
evaluate_markers_roc("GSE74530_expression.csv", "GSE74530_meta.csv",
                     markers = "ZNF662", dataset_name = "GSE74530",
                     save_dir = "plots/ZNF662", formats = c("tiff","pdf"))
evaluate_markers_roc("GSE150469_expression.csv", "GSE150469_meta.csv",
                     markers = "ZNF662", dataset_name = "GSE150469",
                     save_dir = "plots/ZNF662", formats = c("tiff","pdf"))
evaluate_markers_roc("GSE30784_expression.csv", "GSE30784_meta.csv",
                     markers = "ZNF662", dataset_name = "GSE30784",
                     save_dir = "plots/ZNF662", formats = c("tiff","pdf"))
evaluate_markers_roc("GSE246050_expression.csv", "GSE246050_meta.csv",
                     markers = "ZNF662", dataset_name = "GSE246050",
                     save_dir = "plots/ZNF662", formats = c("tiff","pdf"))

# lncRNA (XIST; uses same function)
evaluate_markers_roc("GSE125866_expression.csv", "GSE125866_meta.csv",
                     markers = "XIST", dataset_name = "GSE125866",
                     save_dir = "plots/XIST", formats = c("tiff","pdf"))

# miRNAs (vector)
mirnas <- c("hsa-miR-625-3p","hsa-miR-625-5p","hsa-miR-28-5p",
            "hsa-miR-146-5p","hsa-miR-424-5p","hsa-miR-374a-3p",
            "hsa-miR-16-5p","hsa-miR-15b-5p","hsa-miR-15a-5p")

evaluate_markers_roc("GSE28100_expression.csv","GSE28100_meta.csv",
                     markers = mirnas, dataset_name = "GSE28100",
                     save_dir = "plots/miRNA", formats = c("tiff","pdf"))
