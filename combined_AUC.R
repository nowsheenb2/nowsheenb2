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

# Portable `%||%`-like helper (use right-hand side if left is NULL/NA)
`%OR%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

# Core: compute oriented ROC (AUC >= 0.5), return roc, auc, ci, and best-threshold metrics
.compute_oriented_roc <- function(response_factor,
                                  predictor_numeric,
                                  ci_method = c("delong","bootstrap"),
                                  boot_n = 2000,
                                  boot_stratified = TRUE,
                                  bootstrap_seed = NULL,
                                  quiet = TRUE) {
  ci_method <- match.arg(ci_method)

  # Ensure expected level order: controls first, cases second
  lvls <- c("Normal","Cancer")
  response_factor <- factor(response_factor, levels = lvls)

  # Basic guardrails
  if (any(is.na(response_factor)) ||
      length(unique(na.omit(predictor_numeric))) < 2 ||
      length(unique(na.omit(response_factor))) < 2) {
    return(list(
      roc = NULL, auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      thr = NA_real_, sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_
    ))
  }

  # Compute both orientations
  roc_orig <- tryCatch(roc(response_factor, predictor_numeric, levels = lvls, quiet = quiet), error = function(e) NULL)
  roc_flip <- tryCatch(roc(response_factor, -predictor_numeric, levels = lvls, quiet = quiet), error = function(e) NULL)

  if (is.null(roc_orig) && is.null(roc_flip)) {
    return(list(
      roc = NULL, auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      thr = NA_real_, sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_
    ))
  }

  auc_orig <- if (!is.null(roc_orig)) as.numeric(auc(roc_orig)) else -Inf
  auc_flip <- if (!is.null(roc_flip)) as.numeric(auc(roc_flip)) else -Inf

  if (auc_orig >= auc_flip) { roc_obj <- roc_orig; auc_val <- auc_orig } else { roc_obj <- roc_flip; auc_val <- auc_flip }

  # Confidence interval
  if (ci_method == "delong") {
    ci_val <- tryCatch(ci.auc(roc_obj, method = "delong"), error = function(e) c(NA, NA, NA))
  } else {
    if (!is.null(bootstrap_seed)) set.seed(bootstrap_seed)
    ci_val <- tryCatch(ci.auc(roc_obj, method = "bootstrap",
                              boot.n = boot_n, boot.stratified = boot_stratified),
                       error = function(e) c(NA, NA, NA))
  }
  ci_low <- as.numeric(ci_val[1]); ci_high <- as.numeric(ci_val[3])

  # Best Youden threshold & metrics
  coords_best <- tryCatch(
    coords(roc_obj, "best",
           ret = c("threshold","sensitivity","specificity","ppv","npv"),
           best.method = "youden",
           transpose = TRUE),
    error = function(e) c(threshold = NA, sensitivity = NA, specificity = NA, ppv = NA, npv = NA)
  )

  list(
    roc = roc_obj,
    auc = auc_val,
    ci_low = ci_low,
    ci_high = ci_high,
    thr = unname(coords_best["threshold"]),
    sens = unname(coords_best["sensitivity"]),
    spec = unname(coords_best["specificity"]),
    ppv = unname(coords_best["ppv"]),
    npv = unname(coords_best["npv"])
  )
}

# Plot helper (adds correct diagonal and nice theming)
.plot_roc <- function(roc_obj, title, subtitle = NULL, line_color = "#2c3e50") {
  p <- ggroc(roc_obj, size = 1.2, color = line_color) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  if (!is.null(subtitle)) p <- p + labs(subtitle = subtitle)
  p
}

# Save with multiple formats; TIFF uses LZW compression
.save_plot_multi <- function(plot, filepath_no_ext,
                             formats = c("tiff"),
                             width = 7, height = 6, dpi = 300,
                             tiff_compression = "lzw") {
  formats <- tolower(formats)
  for (fmt in formats) {
    out <- paste0(filepath_no_ext, ".", fmt)
    if (fmt %in% c("tiff","tif")) {
      grDevices::tiff(filename = out, width = width, height = height,
                      units = "in", res = dpi, compression = tiff_compression)
      print(plot)
      grDevices::dev.off()
    } else if (fmt %in% c("pdf","svg","png")) {
      ggsave(out, plot = plot, width = width, height = height, units = "in", dpi = dpi, device = fmt)
    } else {
      stop("Unsupported filetype: ", fmt, ". Use tiff/pdf/svg/png.")
    }
  }
}

# General driver: evaluate one or more markers
# expr_file: CSV with columns: Gene, <Sample1>, <Sample2>, ...
# meta_file: CSV with columns: Sample, Group (values exactly 'Normal' or 'Cancer')
# markers: character vector of gene/lncRNA/miRNA names to evaluate
# dataset_name: used in titles and filenames
evaluate_markers_roc <- function(expr_file, meta_file, markers,
                                 dataset_name = NULL,
                                 save_dir = ".",
                                 formats = c("tiff"),
                                 ci_method = "delong",
                                 boot_n = 2000,
                                 bootstrap_seed = NULL,
                                 write_results_csv = TRUE) {

  expr <- read_csv(expr_file, show_col_types = FALSE)
  meta <- read_csv(meta_file, show_col_types = FALSE)

  # Basic checks
  stopifnot(all(c("Gene") %in% names(expr)))
  stopifnot(all(c("Sample","Group") %in% names(meta)))

  # Factor order: controls (Normal) first, cases (Cancer) second
  meta <- meta %>%
    mutate(Group = factor(Group, levels = c("Normal","Cancer"))) %>%
    distinct(Sample, .keep_all = TRUE)

  # Long-form expression and numeric coercion
  expr_long <- expr %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
    mutate(Expression = suppressWarnings(as.numeric(Expression))) %>%
    filter(!is.na(Expression))

  # Prepare output dir
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # Iterate markers
  results <- map_dfr(markers, function(mk) {
    x <- expr_long %>% filter(Gene == mk)
    if (nrow(x) == 0) {
      message("Marker not found: ", mk, " in ", expr_file, " — skipping")
      return(tibble(
        Marker = mk, AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        Threshold = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
        PPV = NA_real_, NPV = NA_real_,
        N = NA_integer_,
        Dataset = (dataset_name %OR% basename(expr_file)),
        Status = "not_found",
        PlotBase = NA_character_
      ))
    }

    merged <- x %>% inner_join(meta, by = "Sample")

    # Check both groups present
    if (!all(levels(meta$Group) %in% merged$Group)) {
      message("Marker ", mk, ": missing one of the groups after merge — skipping")
      return(tibble(
        Marker = mk, AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        Threshold = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
        PPV = NA_real_, NPV = NA_real_,
        N = nrow(merged),
        Dataset = (dataset_name %OR% basename(expr_file)),
        Status = "missing_group",
        PlotBase = NA_character_
      ))
    }

    # Compute ROC & metrics
    roc_out <- .compute_oriented_roc(
      response_factor = merged$Group,
      predictor_numeric = merged$Expression,
      ci_method = ci_method,
      boot_n = boot_n,
      bootstrap_seed = bootstrap_seed,
      quiet = TRUE
    )

    # If degenerate, record and move on
    if (is.null(roc_out$roc) || is.na(roc_out$auc)) {
      return(tibble(
        Marker = mk, AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
        Threshold = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
        PPV = NA_real_, NPV = NA_real_,
        N = nrow(merged),
        Dataset = (dataset_name %OR% basename(expr_file)),
        Status = "degenerate",
        PlotBase = NA_character_
      ))
    }

    # Plot & save
    auc_str <- sprintf("AUC = %.3f (95%% CI: %.3f–%.3f)", roc_out$auc, roc_out$ci_low, roc_out$ci_high)
    ttl <- if (is.null(dataset_name)) paste0("ROC: ", mk) else paste0("ROC: ", mk, " (", dataset_name, ")")
    p <- .plot_roc(roc_out$roc, title = ttl, subtitle = auc_str)

    safe_marker <- gsub("[^A-Za-z0-9]+", "_", mk)
    safe_ds <- gsub("[^A-Za-z0-9]+", "_", (dataset_name %OR% "dataset"))
    filepath_no_ext <- file.path(save_dir, paste0("ROC_", safe_marker, "_", safe_ds))
    .save_plot_multi(p, filepath_no_ext, formats = formats)

    tibble(
      Marker = mk,
      AUC = roc_out$auc,
      CI_low = roc_out$ci_low,
      CI_high = roc_out$ci_high,
      Threshold = roc_out$thr,
      Sensitivity = roc_out$sens,
      Specificity = roc_out$spec,
      PPV = roc_out$ppv,
      NPV = roc_out$npv,
      N = nrow(merged),
      Dataset = (dataset_name %OR% basename(expr_file)),
      Status = "ok",
      PlotBase = filepath_no_ext
    )
  })

  # Optional CSV summary
  if (isTRUE(write_results_csv)) {
    csv_name <- paste0("ROC_results_", gsub("[^A-Za-z0-9]+","_", dataset_name %OR% basename(expr_file)), ".csv")
    write_csv(results, file.path(save_dir, csv_name))
  }

  results
}

# mRNA (ZNF662 across multiple datasets)
# DeLong CI (default), export TIFF+PDF, write a CSV per call
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

# lncRNA (XIST)
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



