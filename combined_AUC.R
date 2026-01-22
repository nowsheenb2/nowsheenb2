###############################################
## Combined ROC figure:
## a) ZNF662 (4 datasets, 2×2)
## b) XIST (GSE125866)
## c) 9 miRNAs (GSE28100, 3×3)
###############################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(ggplot2)
  library(cowplot)
})

## ---------- Helper functions ---------------------------------------------

# Auto-tune bootstrap for small class sizes
.autotune_bootstrap <- function(y, boot_n, boot_stratified) {
  tab  <- table(y)
  minc <- if (length(tab)) min(tab) else 0
  
  if (minc <= 5) {                    # very small class
    boot_stratified <- FALSE
    boot_n <- min(boot_n, max(100L, 50L * minc))
  } else if (minc <= 10) {            # small-ish
    boot_n <- min(boot_n, max(200L, 80L * minc))
  }
  
  list(
    boot_n          = as.integer(boot_n),
    boot_stratified = isTRUE(boot_stratified)
  )
}

# CI with safe fallback to DeLong
.ci_auc_safe <- function(roc_obj,
                         method = c("delong", "bootstrap"),
                         boot_n = 2000,
                         boot_stratified = TRUE,
                         seed = NULL) {
  method <- match.arg(method)
  if (is.null(roc_obj)) return(c(NA_real_, NA_real_, NA_real_))
  if (!is.null(seed)) set.seed(seed)
  
  ci_try <- tryCatch({
    if (identical(method, "delong")) {
      pROC::ci.auc(roc_obj, method = "delong")
    } else {
      pROC::ci.auc(
        roc_obj,
        method           = "bootstrap",
        boot.n           = boot_n,
        boot.stratified  = boot_stratified
      )
    }
  }, error = function(e) e)
  
  if (!inherits(ci_try, "error")) return(as.numeric(ci_try))
  
  message("[warn] ci.auc(", method, ") failed: ", conditionMessage(ci_try),
          " — falling back to DeLong.")
  
  ci_fb <- tryCatch(
    pROC::ci.auc(roc_obj, method = "delong"),
    error = function(e) c(NA_real_, NA_real_, NA_real_)
  )
  as.numeric(ci_fb)
}

# Plot ROC with nice theme
.plot_roc <- function(roc_obj, title, subtitle = NULL, line_color = "#2c3e50") {
  p <- pROC::ggroc(roc_obj, size = 1.2, color = line_color) +
    ggplot2::geom_abline(
      slope     = 1,
      intercept = 0,
      linetype  = "dashed",
      size      = 0.6
    ) +
    ggplot2::labs(
      title = title,
      x     = "1 - Specificity",
      y     = "Sensitivity"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
  
  if (!is.null(subtitle)) {
    p <- p + ggplot2::labs(subtitle = subtitle)
  }
  
  p
}

# Save plot in multiple formats
.save_plot_multi <- function(plot,
                             filepath_no_ext,
                             formats = c("png", "tiff", "pdf"),
                             width   = 10,
                             height  = 15,
                             dpi     = 600,
                             tiff_compression = "lzw") {
  formats <- tolower(formats)
  
  for (fmt in formats) {
    out <- paste0(filepath_no_ext, ".", fmt)
    if (fmt %in% c("tif", "tiff")) {
      grDevices::tiff(
        filename    = out,
        width       = width,
        height      = height,
        units       = "in",
        res         = dpi,
        compression = tiff_compression
      )
      print(plot)
      grDevices::dev.off()
    } else {
      ggplot2::ggsave(
        filename = out,
        plot     = plot,
        width    = width,
        height   = height,
        units    = "in",
        dpi      = dpi,
        device   = fmt
      )
    }
  }
}

# Core oriented ROC (single run, orientation fixed so AUC ≥ 0.5)
.compute_oriented_roc <- function(response_factor,
                                  predictor_numeric,
                                  ci_method       = c("delong", "bootstrap"),
                                  boot_n          = 2000,
                                  boot_stratified = TRUE,
                                  bootstrap_seed  = NULL,
                                  quiet           = TRUE) {
  ci_method <- match.arg(ci_method)
  lvls <- c("Normal", "Cancer")
  response_factor <- factor(response_factor, levels = lvls)
  
  # basic sanity checks
  if (any(is.na(response_factor)) ||
      length(unique(na.omit(predictor_numeric))) < 2 ||
      length(unique(na.omit(response_factor)))   < 2) {
    return(list(
      roc          = NULL,
      auc          = NA_real_,
      ci_low       = NA_real_,
      ci_high      = NA_real_,
      thr          = NA_real_,
      thr_original = NA_real_,
      sens         = NA_real_,
      spec         = NA_real_,
      ppv          = NA_real_,
      npv          = NA_real_
    ))
  }
  
  roc_orig <- tryCatch(
    pROC::roc(response_factor, predictor_numeric,
              levels = lvls, quiet = quiet),
    error = function(e) NULL
  )
  roc_flip <- tryCatch(
    pROC::roc(response_factor, -predictor_numeric,
              levels = lvls, quiet = quiet),
    error = function(e) NULL
  )
  
  if (is.null(roc_orig) && is.null(roc_flip)) {
    return(list(
      roc          = NULL,
      auc          = NA_real_,
      ci_low       = NA_real_,
      ci_high      = NA_real_,
      thr          = NA_real_,
      thr_original = NA_real_,
      sens         = NA_real_,
      spec         = NA_real_,
      ppv          = NA_real_,
      npv          = NA_real_
    ))
  }
  
  auc_orig <- if (!is.null(roc_orig)) as.numeric(pROC::auc(roc_orig)) else -Inf
  auc_flip <- if (!is.null(roc_flip)) as.numeric(pROC::auc(roc_flip)) else -Inf
  
  used_flip <- isTRUE(auc_flip > auc_orig)
  roc_obj   <- if (used_flip) roc_flip else roc_orig
  auc_val   <- if (used_flip) auc_flip else auc_orig
  
  # CI (auto-tuned)
  tuned  <- .autotune_bootstrap(response_factor, boot_n, boot_stratified)
  ci_val <- .ci_auc_safe(
    roc_obj,
    method          = ci_method,
    boot_n          = tuned$boot_n,
    boot_stratified = tuned$boot_stratified,
    seed            = bootstrap_seed
  )
  ci_low  <- as.numeric(ci_val[1])
  ci_high <- as.numeric(ci_val[3])
  
  coords_best <- tryCatch(
    pROC::coords(
      roc_obj,
      "best",
      ret          = c("threshold","sensitivity","specificity","ppv","npv"),
      best.method  = "youden",
      transpose    = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(coords_best)) {
    thr_used <- sens <- spec <- ppv <- npv <- NA_real_
  } else {
    thr_used <- as.numeric(coords_best[1, "threshold"])
    if (!is.finite(thr_used)) thr_used <- NA_real_
    sens <- as.numeric(coords_best[1, "sensitivity"])
    spec <- as.numeric(coords_best[1, "specificity"])
    ppv  <- as.numeric(coords_best[1, "ppv"])
    npv  <- as.numeric(coords_best[1, "npv"])
  }
  
  thr_original <- if (is.na(thr_used)) {
    NA_real_
  } else if (used_flip) {
    -thr_used
  } else {
    thr_used
  }
  
  list(
    roc          = roc_obj,
    auc          = auc_val,
    ci_low       = ci_low,
    ci_high      = ci_high,
    thr          = thr_used,
    thr_original = thr_original,
    sens         = sens,
    spec         = spec,
    ppv          = ppv,
    npv          = npv
  )
}

# Build one ROC ggplot for a marker in one dataset
# NOTE: title now uses two lines so "(GSE28100)" is fully visible
make_marker_roc_plot <- function(expr_file,
                                 meta_file,
                                 marker,
                                 dataset_name  = NULL,
                                 ci_method     = "delong",
                                 boot_n        = 2000,
                                 boot_stratified = TRUE) {
  expr <- readr::read_csv(expr_file, show_col_types = FALSE)
  meta <- readr::read_csv(meta_file, show_col_types = FALSE)
  
  stopifnot("Gene" %in% colnames(expr))
  stopifnot(all(c("Sample","Group") %in% colnames(meta)))
  
  meta$Group <- factor(meta$Group, levels = c("Normal","Cancer"))
  meta <- dplyr::distinct(meta, Sample, .keep_all = TRUE)
  
  expr_long <- tidyr::pivot_longer(
    expr,
    cols      = -Gene,
    names_to  = "Sample",
    values_to = "Expression"
  )
  expr_long$Expression <- suppressWarnings(as.numeric(expr_long$Expression))
  expr_long <- expr_long[!is.na(expr_long$Expression), , drop = FALSE]
  
  x <- expr_long[expr_long$Gene == marker, , drop = FALSE]
  if (nrow(x) == 0L) {
    stop("Marker ", marker, " not found in ", expr_file)
  }
  
  merged <- dplyr::inner_join(x, meta, by = "Sample")
  merged <- merged[!is.na(merged$Group) & is.finite(merged$Expression), , drop = FALSE]
  
  if (dplyr::n_distinct(merged$Group) < 2) {
    stop("Only one class present for ", marker,
         " in ", ifelse(is.null(dataset_name), basename(expr_file), dataset_name))
  }
  
  roc_out <- .compute_oriented_roc(
    response_factor   = merged$Group,
    predictor_numeric = merged$Expression,
    ci_method         = ci_method,
    boot_n            = boot_n,
    boot_stratified   = boot_stratified,
    bootstrap_seed    = NULL,
    quiet             = TRUE
  )
  
  if (is.null(roc_out$roc) || is.na(roc_out$auc)) {
    stop("Degenerate ROC for ", marker,
         " in ", ifelse(is.null(dataset_name), basename(expr_file), dataset_name))
  }
  
  auc_str <- sprintf(
    "AUC = %.3f (95%% CI: %.3f–%.3f)",
    roc_out$auc, roc_out$ci_low, roc_out$ci_high
  )
  
  # DATASET ON NEW LINE TO AVOID TRUNCATION (e.g. GSE28100)
  if (is.null(dataset_name)) {
    ttl <- paste0("ROC: ", marker)
  } else {
    ttl <- paste0("ROC: ", marker, "\n(", dataset_name, ")")
  }
  
  .plot_roc(roc_out$roc, title = ttl, subtitle = auc_str)
}

## ---------- Build combined (a,b,c) figure ------------------------------

run_combined_roc_abc <- function(out_dir = "plots_cv/Combined_ROC_abc") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## a) ZNF662 – 4 datasets, 2×2 grid
  znf_cfg <- list(
    list(expr = "GSE74530_expression.csv",   meta = "GSE74530_meta.csv",   ds = "GSE74530"),
    list(expr = "GSE150469_expression.csv",  meta = "GSE150469_meta.csv",  ds = "GSE150469"),
    list(expr = "GSE30784_expression.csv",   meta = "GSE30784_meta.csv",   ds = "GSE30784"),
    list(expr = "GSE246050_expression.csv",  meta = "GSE246050_meta.csv",  ds = "GSE246050")
  )
  
  znf_plots <- lapply(znf_cfg, function(cfg) {
    make_marker_roc_plot(
      expr_file       = cfg$expr,
      meta_file       = cfg$meta,
      marker          = "ZNF662",
      dataset_name    = cfg$ds,
      ci_method       = "delong",
      boot_n          = 2000,
      boot_stratified = TRUE
    )
  })
  
  panel_a <- cowplot::plot_grid(
    plotlist = znf_plots,
    ncol     = 2
  )
  
  ## b) XIST – single ROC
  panel_b <- make_marker_roc_plot(
    expr_file       = "GSE125866_expression.csv",
    meta_file       = "GSE125866_meta.csv",
    marker          = "XIST",
    dataset_name    = "GSE125866",
    ci_method       = "delong",
    boot_n          = 2000,
    boot_stratified = TRUE
  )
  
  ## c) miRNAs – 9 ROC curves, 3×3 grid
  mirnas <- c(
    "hsa-miR-625-3p", "hsa-miR-625-5p", "hsa-miR-28-5p",
    "hsa-miR-146b-5p","hsa-miR-424-5p", "hsa-miR-374a-3p",
    "hsa-miR-16-5p",  "hsa-miR-15b-5p", "hsa-miR-15a-5p"
  )
  
  mirna_plots <- lapply(mirnas, function(mn) {
    make_marker_roc_plot(
      expr_file       = "GSE28100_expression.csv",
      meta_file       = "GSE28100_meta.csv",
      marker          = mn,
      dataset_name    = "GSE28100",
      ci_method       = "delong",
      boot_n          = 2000,
      boot_stratified = TRUE
    )
  })
  
  panel_c <- cowplot::plot_grid(
    plotlist = mirna_plots,
    ncol     = 3
  )
  
  ## Combine a, b, c with labels
  combined_roc_abc <- cowplot::plot_grid(
    panel_a,
    panel_b,
    panel_c,
    labels      = c("a", "b", "c"),
    label_size  = 18,
    ncol        = 1,
    rel_heights = c(1.3, 0.7, 1.6)   # tweak if you want more/less space
  )
  
  .save_plot_multi(
    plot            = combined_roc_abc,
    filepath_no_ext = file.path(out_dir, "ROC_ZNF662_XIST_miRNAs_abc"),
    formats         = c("png", "tiff", "pdf"),
    width           = 10,
    height          = 15,
    dpi             = 600
  )
  
  message("✅ Combined ROC figure saved in: ",
          normalizePath(out_dir, mustWork = FALSE))
}

## ---------- Run to generate the figure -----------------------

run_combined_roc_abc()
