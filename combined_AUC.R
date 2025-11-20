# ====================== ROC + K-fold CV utilities (auto-tuned bootstrap) ======================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(pROC)
  library(ggplot2); library(purrr); library(stringr)
})

`%OR%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x
.safe_div <- function(x, y) ifelse(y > 0, x / y, NA_real_)
.mode_chr <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_character_); ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

# --- NEW: auto-tune bootstrap for tiny classes/folds ---
.autotune_bootstrap <- function(y, boot_n, boot_stratified) {
  tab <- table(y); minc <- if (length(tab)) min(tab) else 0
  if (minc <= 5) {                        # very small class ⇒ lighten bootstrap
    boot_stratified <- FALSE
    boot_n <- min(boot_n, max(100L, 50L * minc))
  } else if (minc <= 10) {                # smallish ⇒ moderate
    boot_n <- min(boot_n, max(200L, 80L * minc))
  }
  list(boot_n = as.integer(boot_n), boot_stratified = isTRUE(boot_stratified))
}

# Safe CI with fallback to DeLong (never aborts)
.ci_auc_safe <- function(roc_obj, method = c("delong","bootstrap"),
                         boot_n = 2000, boot_stratified = TRUE, seed = NULL) {
  method <- match.arg(method)
  if (is.null(roc_obj)) return(c(NA_real_, NA_real_, NA_real_))
  if (!is.null(seed)) set.seed(seed)
  ci_try <- tryCatch({
    if (identical(method, "delong")) {
      pROC::ci.auc(roc_obj, method = "delong")
    } else {
      pROC::ci.auc(roc_obj, method = "bootstrap",
                   boot.n = boot_n, boot.stratified = boot_stratified)
    }
  }, error = function(e) e)
  if (!inherits(ci_try, "error")) return(as.numeric(ci_try))
  message(sprintf("[warn] ci.auc(%s) failed: %s — falling back to DeLong.",
                  method, conditionMessage(ci_try)))
  ci_fb <- tryCatch(pROC::ci.auc(roc_obj, method = "delong"),
                    error = function(e) c(NA_real_, NA_real_, NA_real_))
  as.numeric(ci_fb)
}

.plot_roc <- function(roc_obj, title, subtitle = NULL, line_color = "#2c3e50") {
  p <- pROC::ggroc(roc_obj, size = 1.2, color = line_color) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.6) +
    ggplot2::labs(title = title, x = "1 - Specificity", y = "Sensitivity") +
    ggplot2::theme_minimal(base_size = 14) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  if (!is.null(subtitle)) p <- p + ggplot2::labs(subtitle = subtitle)
  p
}

.save_plot_multi <- function(plot, filepath_no_ext,
                             formats = c("png","tiff","pdf"),
                             width = 7, height = 6, dpi = 300,
                             tiff_compression = "lzw") {
  formats <- tolower(formats)
  for (fmt in formats) {
    out <- paste0(filepath_no_ext, ".", fmt)
    tryCatch({
      if (fmt %in% c("tiff","tif")) {
        grDevices::tiff(filename = out, width = width, height = height,
                        units = "in", res = dpi, compression = tiff_compression)
        print(plot); grDevices::dev.off()
      } else if (fmt %in% c("pdf","svg","png")) {
        if (fmt == "svg" && !requireNamespace("svglite", quietly = TRUE)) {
          warning("Format 'svg' requires 'svglite'; skipping.")
        } else {
          ggplot2::ggsave(out, plot = plot, width = width, height = height,
                          units = "in", dpi = dpi, device = fmt)
        }
      } else stop("Unsupported filetype: ", fmt)
    }, error = function(e) {
      warning("Failed to save ", out, ": ", conditionMessage(e))
      try({ if (grDevices::dev.cur() > 1) grDevices::dev.off() }, silent = TRUE)
    })
  }
}

# ---------- Core oriented ROC (single run) ----------
.compute_oriented_roc <- function(response_factor, predictor_numeric,
                                  ci_method = c("delong","bootstrap"),
                                  boot_n = 2000, boot_stratified = TRUE,
                                  bootstrap_seed = NULL, quiet = TRUE) {
  ci_method <- match.arg(ci_method)
  lvls <- c("Normal","Cancer")
  response_factor <- factor(response_factor, levels = lvls)
  if (any(is.na(response_factor)) ||
      length(unique(na.omit(predictor_numeric))) < 2 ||
      length(unique(na.omit(response_factor))) < 2) {
    return(list(roc=NULL, auc=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
                thr=NA_real_, thr_original=NA_real_, sens=NA_real_, spec=NA_real_,
                ppv=NA_real_, npv=NA_real_, used_flip=NA, direction=NA_character_,
                rule=NA_character_, n_ctrl=sum(response_factor=="Normal",na.rm=TRUE),
                n_case=sum(response_factor=="Cancer",na.rm=TRUE), prevalence=NA_real_))
  }
  roc_orig <- tryCatch(pROC::roc(response_factor, predictor_numeric, levels = lvls, quiet = quiet), error=function(e) NULL)
  roc_flip <- tryCatch(pROC::roc(response_factor, -predictor_numeric, levels = lvls, quiet = quiet), error=function(e) NULL)
  if (is.null(roc_orig) && is.null(roc_flip)) {
    return(list(roc=NULL, auc=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
                thr=NA_real_, thr_original=NA_real_, sens=NA_real_, spec=NA_real_,
                ppv=NA_real_, npv=NA_real_, used_flip=NA, direction=NA_character_,
                rule=NA_character_, n_ctrl=sum(response_factor=="Normal",na.rm=TRUE),
                n_case=sum(response_factor=="Cancer",na.rm=TRUE), prevalence=NA_real_))
  }
  auc_orig <- if (!is.null(roc_orig)) as.numeric(pROC::auc(roc_orig)) else -Inf
  auc_flip <- if (!is.null(roc_flip)) as.numeric(pROC::auc(roc_flip)) else -Inf
  used_flip <- isTRUE(auc_flip > auc_orig)
  roc_obj <- if (used_flip) roc_flip else roc_orig
  auc_val <- if (used_flip) auc_flip else auc_orig
  
  # --- auto-tune bootstrap just for CI (uses true labels) ---
  tuned <- .autotune_bootstrap(response_factor, boot_n, boot_stratified)
  ci_val <- .ci_auc_safe(roc_obj, method = ci_method,
                         boot_n = tuned$boot_n, boot_stratified = tuned$boot_stratified,
                         seed = bootstrap_seed)
  ci_low <- as.numeric(ci_val[1]); ci_high <- as.numeric(ci_val[3])
  
  coords_best <- tryCatch(
    pROC::coords(roc_obj, "best",
                 ret = c("threshold","sensitivity","specificity","ppv","npv"),
                 best.method = "youden", transpose = FALSE),
    error = function(e) NULL
  )
  if (is.null(coords_best)) {
    thr_used <- NA_real_; sens <- spec <- ppv <- npv <- NA_real_
  } else {
    thr_used <- as.numeric(coords_best[1,"threshold"]); if (!is.finite(thr_used)) thr_used <- NA_real_
    sens <- as.numeric(coords_best[1,"sensitivity"]); spec <- as.numeric(coords_best[1,"specificity"])
    ppv  <- as.numeric(coords_best[1,"ppv"]);        npv  <- as.numeric(coords_best[1,"npv"])
  }
  thr_original <- if (is.na(thr_used)) NA_real_ else if (used_flip) -thr_used else thr_used
  dir_used <- if (is.null(roc_obj)) NA_character_ else as.character(roc_obj$direction)
  direction_original <- if (is.na(dir_used)) NA_character_ else { if (used_flip) if (dir_used==">") "<" else ">" else dir_used }
  rule <- if (is.na(thr_original) || is.na(direction_original)) NA_character_
  else paste0("Predict Cancer if Expression ", direction_original, " ", signif(thr_original, 5))
  n_ctrl <- sum(response_factor=="Normal", na.rm=TRUE); n_case <- sum(response_factor=="Cancer", na.rm=TRUE)
  prevalence <- if ((n_ctrl + n_case) > 0) n_case / (n_ctrl + n_case) else NA_real_
  
  list(roc=roc_obj, auc=auc_val, ci_low=ci_low, ci_high=ci_high,
       thr=thr_used, thr_original=thr_original, sens=sens, spec=spec, ppv=ppv, npv=npv,
       used_flip=used_flip, direction=direction_original, rule=rule,
       n_ctrl=n_ctrl, n_case=n_case, prevalence=prevalence)
}

# ---------- Public: single-run ROC for markers ----------
evaluate_markers_roc <- function(expr_file, meta_file, markers,
                                 dataset_name = NULL, save_dir = ".",
                                 formats = c("png","tiff","pdf"),
                                 ci_method = "delong", boot_n = 2000,
                                 bootstrap_seed = NULL, boot_stratified = TRUE,
                                 write_results_csv = TRUE) {
  expr <- readr::read_csv(expr_file, show_col_types = FALSE)
  meta <- readr::read_csv(meta_file, show_col_types = FALSE)
  stopifnot(all(c("Gene") %in% names(expr)))
  stopifnot(all(c("Sample","Group") %in% names(meta)))
  meta$Group <- factor(meta$Group, levels = c("Normal","Cancer"))
  meta <- dplyr::distinct(meta, Sample, .keep_all = TRUE)
  
  expr_long <- tidyr::pivot_longer(expr, -Gene, names_to = "Sample", values_to = "Expression")
  expr_long$Expression <- suppressWarnings(as.numeric(expr_long$Expression))
  expr_long <- expr_long[!is.na(expr_long$Expression), , drop = FALSE]
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  results_list <- lapply(markers, function(mk) {
    x <- expr_long[expr_long$Gene == mk, , drop = FALSE]
    if (!nrow(x)) return(tibble::tibble(
      Marker=mk, AUC=NA_real_, CI_low=NA_real_, CI_high=NA_real_,
      Threshold=NA_real_, Sensitivity=NA_real_, Specificity=NA_real_, PPV=NA_real_, NPV=NA_real_,
      N=NA_integer_, Dataset=(dataset_name %OR% basename(expr_file)),
      Status="not_found", PlotBase=NA_character_, DecisionRule=NA_character_, Prevalence=NA_real_
    ))
    
    merged <- dplyr::inner_join(x, meta, by="Sample")
    merged <- merged[!is.na(merged$Group) & is.finite(merged$Expression), , drop = FALSE]
    if (dplyr::n_distinct(merged$Group) < 2) return(tibble::tibble(
      Marker=mk, AUC=NA_real_, CI_low=NA_real_, CI_high=NA_real_,
      Threshold=NA_real_, Sensitivity=NA_real_, Specificity=NA_real_, PPV=NA_real_, NPV=NA_real_,
      N=nrow(merged), Dataset=(dataset_name %OR% basename(expr_file)),
      Status="missing_group", PlotBase=NA_character_, DecisionRule=NA_character_, Prevalence=NA_real_
    ))
    
    # compute ROC (with auto-tuned bootstrap)
    tuned <- .autotune_bootstrap(merged$Group, boot_n, boot_stratified)
    roc_out <- .compute_oriented_roc(
      response_factor = merged$Group, predictor_numeric = merged$Expression,
      ci_method = ci_method, boot_n = tuned$boot_n,
      boot_stratified = tuned$boot_stratified, bootstrap_seed = bootstrap_seed, quiet = TRUE
    )
    if (is.null(roc_out$roc) || is.na(roc_out$auc)) return(tibble::tibble(
      Marker=mk, AUC=NA_real_, CI_low=NA_real_, CI_high=NA_real_,
      Threshold=NA_real_, Sensitivity=NA_real_, Specificity=NA_real_, PPV=NA_real_, NPV=NA_real_,
      N=nrow(merged), Dataset=(dataset_name %OR% basename(expr_file)),
      Status="degenerate", PlotBase=NA_character_, DecisionRule=NA_character_, Prevalence=NA_real_
    ))
    
    auc_str <- sprintf("AUC = %.3f (95%% CI: %.3f–%.3f)", roc_out$auc, roc_out$ci_low, roc_out$ci_high)
    ttl <- if (is.null(dataset_name)) paste0("ROC: ", mk) else paste0("ROC: ", mk, " (", dataset_name, ")")
    p <- .plot_roc(roc_out$roc, title = ttl, subtitle = auc_str)
    
    safe_marker <- gsub("[^A-Za-z0-9]+", "_", mk)
    safe_ds <- gsub("[^A-Za-z0-9]+", "_", (dataset_name %OR% "dataset"))
    filepath_no_ext <- file.path(save_dir, paste0("ROC_", safe_marker, "_", safe_ds))
    .save_plot_multi(p, filepath_no_ext, formats = formats)
    
    tibble::tibble(Marker=mk, AUC=roc_out$auc, CI_low=roc_out$ci_low, CI_high=roc_out$ci_high,
                   Threshold=roc_out$thr_original, Sensitivity=roc_out$sens, Specificity=roc_out$spec,
                   PPV=roc_out$ppv, NPV=roc_out$npv, N=nrow(merged),
                   Dataset=(dataset_name %OR% basename(expr_file)), Status="ok",
                   PlotBase=filepath_no_ext, DecisionRule=roc_out$rule, Prevalence=roc_out$prevalence)
  })
  
  results <- dplyr::bind_rows(results_list)
  if (isTRUE(write_results_csv)) {
    csv_name <- paste0("ROC_results_", gsub("[^A-Za-z0-9]+","_", dataset_name %OR% basename(expr_file)), ".csv")
    readr::write_csv(results, file.path(save_dir, csv_name))
  }
  results
}

# ---------- K-fold CV helpers ----------
.make_stratified_folds <- function(y, k = 5, seed = NULL) {
  y <- factor(y, levels = c("Normal","Cancer"))
  if (!is.null(seed)) set.seed(seed)
  n <- length(y); folds <- integer(n); tab <- table(y)
  k_eff <- max(2, min(k, max(1, min(tab))))
  if (k_eff < k) message("Reducing k from ", k, " to ", k_eff, " due to small class counts."); k <- k_eff
  for (lev in levels(y)) {
    idx <- which(y == lev); if (!length(idx)) next
    idx <- sample(idx); bucket <- rep(1:k, length.out = length(idx)); folds[idx] <- bucket
  }
  folds
}

.train_oriented_rule <- function(group_train, x_train,
                                 ci_method = c("delong","bootstrap"),
                                 boot_n = 2000, boot_stratified = TRUE,
                                 seed = NULL, quiet = TRUE) {
  lvls <- c("Normal","Cancer")
  group_train <- factor(group_train, levels = lvls)
  roc_orig <- tryCatch(pROC::roc(group_train, x_train, levels = lvls, quiet = quiet), error=function(e) NULL)
  roc_flip <- tryCatch(pROC::roc(group_train, -x_train, levels = lvls, quiet = quiet), error=function(e) NULL)
  if (is.null(roc_orig) && is.null(roc_flip)) return(NULL)
  auc_orig <- if (!is.null(roc_orig)) as.numeric(pROC::auc(roc_orig)) else -Inf
  auc_flip <- if (!is.null(roc_flip)) as.numeric(pROC::auc(roc_flip)) else -Inf
  use_flip <- auc_flip > auc_orig
  roc_train <- if (use_flip) roc_flip else roc_orig
  
  coords_best <- tryCatch(pROC::coords(roc_train, "best",
                                       ret=c("threshold","sensitivity","specificity"),
                                       best.method="youden", transpose=FALSE), error=function(e) NULL)
  thr_used <- if (is.null(coords_best)) NA_real_ else as.numeric(coords_best[1,"threshold"])
  if (!is.finite(thr_used)) thr_used <- NA_real_
  thr_original <- if (isTRUE(use_flip) && !is.na(thr_used)) -thr_used else thr_used
  dir_used <- as.character(roc_train$direction)
  direction_original <- if (isTRUE(use_flip)) if (dir_used==">") "<" else ">" else dir_used
  
  apply_rule <- function(x) {
    if (is.na(thr_original) || is.na(direction_original)) return(rep(NA_integer_, length(x)))
    if (direction_original == ">") as.integer(x >= thr_original) else as.integer(x <= thr_original)
  }
  orient_score <- function(x) if (isTRUE(use_flip)) -x else x
  
  list(used_flip=isTRUE(use_flip), thr_original=thr_original, direction_original=direction_original,
       roc_train=roc_train, apply_rule=apply_rule, orient_score=orient_score)
}

.fold_confusion <- function(y_true, y_pred01) {
  y_true <- factor(y_true, levels = c("Normal","Cancer"))
  tp <- sum(y_pred01 == 1 & y_true == "Cancer", na.rm = TRUE)
  tn <- sum(y_pred01 == 0 & y_true == "Normal", na.rm = TRUE)
  fp <- sum(y_pred01 == 1 & y_true == "Normal", na.rm = TRUE)
  fn <- sum(y_pred01 == 0 & y_true == "Cancer", na.rm = TRUE)
  list(TP=tp, FP=fp, TN=tn, FN=fn)
}

# ---------- K-fold CV core ----------
.kfold_cv_run <- function(merged, marker_name, k = 5, repeats = 1, seed = NULL,
                          ci_method = c("delong","bootstrap"),
                          boot_n = 2000, boot_stratified = TRUE,
                          quiet = TRUE, use_global_threshold = c("per_fold","median_train"),
                          target_prevalence = NULL) {
  ci_method <- match.arg(ci_method)
  use_global_threshold <- match.arg(use_global_threshold)
  lvls <- c("Normal","Cancer")
  merged$Group <- factor(merged$Group, levels = lvls)
  merged <- merged[!is.na(merged$Group) & is.finite(merged$Expression), , drop = FALSE]
  class_tab <- table(merged$Group)
  if (length(class_tab) < 2 || min(class_tab) < 2) {
    return(list(status="insufficient", marker=marker_name, cv_auc=NA_real_, cv_ci_low=NA_real_, cv_ci_high=NA_real_,
                sens=NA_real_, spec=NA_real_, ppv=NA_real_, npv=NA_real_, k=k, repeats=repeats,
                fold_results=tibble::tibble(), cv_roc=NULL, global_thr=NA_real_, global_dir=NA_character_,
                ppv_adj=NA_real_, npv_adj=NA_real_, target_prev=target_prevalence))
  }
  
  oof_scores <- numeric(0); oof_truth <- factor(character(0), levels = lvls); oof_expr <- numeric(0)
  agg <- list(TP=0L, FP=0L, TN=0L, FN=0L); train_thresholds <- c(); train_dirs <- character()
  fold_rows <- list(); fold_id <- 1L
  
  for (rep_i in seq_len(repeats)) {
    seed_rep <- if (is.null(seed)) NULL else (seed + rep_i - 1L)
    folds <- .make_stratified_folds(merged$Group, k = k, seed = seed_rep)
    for (fold in sort(unique(folds))) {
      idx_test  <- which(folds == fold); idx_train <- which(folds != fold)
      train <- merged[idx_train, , drop = FALSE]; test <- merged[idx_test, , drop = FALSE]
      tr <- .train_oriented_rule(train$Group, train$Expression,
                                 ci_method=ci_method, boot_n=boot_n,
                                 boot_stratified=boot_stratified, seed=seed_rep, quiet=quiet)
      if (is.null(tr)) {
        fold_rows[[fold_id]] <- tibble::tibble(Repeat=rep_i, Fold=fold, UsedFlip=NA,
                                               Thr=NA_real_, Dir=NA_character_, N_train=nrow(train), N_test=nrow(test),
                                               Prev_test=mean(test$Group=="Cancer"), AUC_test=NA_real_, Sens_test=NA_real_,
                                               Spec_test=NA_real_, PPV_test=NA_real_, NPV_test=NA_real_)
        fold_id <- fold_id + 1L; next
      }
      train_thresholds <- c(train_thresholds, tr$thr_original); train_dirs <- c(train_dirs, tr$direction_original)
      pred01 <- tr$apply_rule(test$Expression); cf <- .fold_confusion(test$Group, pred01)
      sens <- .safe_div(cf$TP, cf$TP + cf$FN); spec <- .safe_div(cf$TN, cf$TN + cf$FP)
      ppv  <- .safe_div(cf$TP, cf$TP + cf$FP); npv <- .safe_div(cf$TN, cf$TN + cf$FN)
      test_scores <- tr$orient_score(test$Expression)
      auc_test <- tryCatch({
        if (dplyr::n_distinct(test$Group) < 2) NA_real_
        else as.numeric(pROC::auc(pROC::roc(test$Group, test_scores, levels = lvls, quiet = TRUE)))
      }, error=function(e) NA_real_)
      oof_scores <- c(oof_scores, test_scores)
      oof_truth  <- factor(c(as.character(oof_truth), as.character(test$Group)), levels = lvls)
      oof_expr   <- c(oof_expr, test$Expression)
      agg$TP <- agg$TP + cf$TP; agg$FP <- agg$FP + cf$FP; agg$TN <- agg$TN + cf$TN; agg$FN <- agg$FN + cf$FN
      fold_rows[[fold_id]] <- tibble::tibble(Repeat=rep_i, Fold=fold, UsedFlip=tr$used_flip, Thr=tr$thr_original,
                                             Dir=tr$direction_original, N_train=nrow(train), N_test=nrow(test),
                                             Prev_test=mean(test$Group=="Cancer"), AUC_test=auc_test,
                                             Sens_test=sens, Spec_test=spec, PPV_test=ppv, NPV_test=npv)
      fold_id <- fold_id + 1L
    }
  }
  
  fold_results <- if (length(fold_rows)) dplyr::bind_rows(fold_rows) else tibble::tibble()
  cv_roc <- tryCatch(pROC::roc(oof_truth, oof_scores, levels = lvls, quiet = TRUE), error=function(e) NULL)
  if (is.null(cv_roc)) {
    return(list(status="degenerate", marker=marker_name, cv_auc=NA_real_, cv_ci_low=NA_real_, cv_ci_high=NA_real_,
                sens=NA_real_, spec=NA_real_, ppv=NA_real_, npv=NA_real_, k=k, repeats=repeats,
                fold_results=fold_results, cv_roc=NULL, global_thr=NA_real_, global_dir=NA_character_,
                ppv_adj=NA_real_, npv_adj=NA_real_, target_prev=target_prevalence))
  }
  
  # --- auto-tune bootstrap for CV using all OOF truth ---
  tuned_cv <- .autotune_bootstrap(oof_truth, boot_n, boot_stratified)
  ci_val <- .ci_auc_safe(cv_roc, method = ci_method,
                         boot_n = tuned_cv$boot_n, boot_stratified = tuned_cv$boot_stratified,
                         seed = seed)
  cv_auc <- as.numeric(pROC::auc(cv_roc)); cv_ci_low <- as.numeric(ci_val[1]); cv_ci_high <- as.numeric(ci_val[3])
  
  sens_cv <- .safe_div(agg$TP, agg$TP + agg$FN); spec_cv <- .safe_div(agg$TN, agg$TN + agg$FP)
  ppv_cv  <- .safe_div(agg$TP, agg$TP + agg$FP); npv_cv <- .safe_div(agg$TN, agg$TN + agg$FN)
  
  global_thr <- NA_real_; global_dir <- NA_character_
  if (identical(use_global_threshold, "median_train") && length(train_thresholds)) {
    global_thr <- stats::median(train_thresholds, na.rm = TRUE)
    global_dir <- .mode_chr(train_dirs); if (is.na(global_dir)) global_dir <- ">"
    dir_mode <- .mode_chr(train_dirs); dir_agreement <- mean(train_dirs == dir_mode, na.rm = TRUE)
    if (is.nan(dir_agreement)) dir_agreement <- NA_real_
    if (!is.na(dir_agreement) && dir_agreement < 0.6) {
      message(sprintf("Warning: low agreement on direction across folds (%.0f%%).", 100 * dir_agreement))
    }
    pred01_g <- if (global_dir == ">") as.integer(oof_expr >= global_thr) else as.integer(oof_expr <= global_thr)
    cf_g <- .fold_confusion(oof_truth, pred01_g)
    sens_cv <- .safe_div(cf_g$TP, cf_g$TP + cf_g$FN); spec_cv <- .safe_div(cf_g$TN, cf_g$TN + cf_g$FP)
    ppv_cv  <- .safe_div(cf_g$TP, cf_g$TP + cf_g$FP); npv_cv <- .safe_div(cf_g$TN, cf_g$TN + cf_g$FN)
  }
  
  ppv_adj <- npv_adj <- NA_real_
  if (!is.null(target_prevalence) && is.finite(sens_cv) && is.finite(spec_cv)) {
    p <- target_prevalence
    ppv_adj <- (sens_cv * p) / (sens_cv * p + (1 - spec_cv) * (1 - p))
    npv_adj <- (spec_cv * (1 - p)) / ((1 - sens_cv) * p + spec_cv * (1 - p))
  }
  
  list(status="ok", marker=marker_name, cv_auc=cv_auc, cv_ci_low=cv_ci_low, cv_ci_high=cv_ci_high,
       sens=sens_cv, spec=spec_cv, ppv=ppv_cv, npv=npv_cv, k=k, repeats=repeats,
       fold_results=fold_results, cv_roc=cv_roc, global_thr=global_thr, global_dir=global_dir,
       ppv_adj=ppv_adj, npv_adj=npv_adj, target_prev=target_prevalence)
}

# ---------- Public: CV driver ----------
evaluate_markers_roc_cv <- function(expr_file, meta_file, markers,
                                    dataset_name = NULL, k = 5, repeats = 1, seed = 123,
                                    save_dir = ".", formats = c("png","tiff","pdf"),
                                    ci_method = "delong", boot_n = 2000,
                                    write_results_csv = TRUE, boot_stratified = TRUE,
                                    use_global_threshold = c("per_fold","median_train"),
                                    target_prevalence = NULL) {
  use_global_threshold <- match.arg(use_global_threshold)
  expr <- readr::read_csv(expr_file, show_col_types = FALSE)
  meta <- readr::read_csv(meta_file,  show_col_types = FALSE)
  stopifnot(all(c("Gene") %in% names(expr)))
  stopifnot(all(c("Sample","Group") %in% names(meta)))
  meta$Group <- factor(meta$Group, levels = c("Normal","Cancer"))
  meta <- dplyr::distinct(meta, Sample, .keep_all = TRUE)
  
  expr_long <- tidyr::pivot_longer(expr, -Gene, names_to = "Sample", values_to = "Expression")
  expr_long$Expression <- suppressWarnings(as.numeric(expr_long$Expression))
  expr_long <- expr_long[!is.na(expr_long$Expression), , drop = FALSE]
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  out_rows <- purrr::map_dfr(markers, function(mk) {
    x <- expr_long[expr_long$Gene == mk, , drop = FALSE]
    if (!nrow(x)) return(tibble::tibble(
      Marker=mk, CV_AUC=NA_real_, CV_CI_low=NA_real_, CV_CI_high=NA_real_,
      CV_Sensitivity=NA_real_, CV_Specificity=NA_real_, CV_PPV=NA_real_, CV_NPV=NA_real_,
      K=k, Repeats=repeats, N=NA_integer_, Dataset=(dataset_name %OR% basename(expr_file)),
      Status="not_found", CVPlotBase=NA_character_, ThresholdMode=use_global_threshold,
      GlobalThr=NA_real_, GlobalDir=NA_character_, CV_PPV_adj=NA_real_, CV_NPV_adj=NA_real_, TargetPrev=target_prevalence
    ))
    
    merged <- dplyr::inner_join(x, meta, by = "Sample")
    merged <- merged[!is.na(merged$Group) & is.finite(merged$Expression), , drop = FALSE]
    class_tab <- table(merged$Group)
    if (length(class_tab) < 2 || min(class_tab) < 2) return(tibble::tibble(
      Marker=mk, CV_AUC=NA_real_, CV_CI_low=NA_real_, CV_CI_high=NA_real_,
      CV_Sensitivity=NA_real_, CV_Specificity=NA_real_, CV_PPV=NA_real_, CV_NPV=NA_real_,
      K=k, Repeats=repeats, N=nrow(merged), Dataset=(dataset_name %OR% basename(expr_file)),
      Status="insufficient", CVPlotBase=NA_character_, ThresholdMode=use_global_threshold,
      GlobalThr=NA_real_, GlobalDir=NA_character_, CV_PPV_adj=NA_real_, CV_NPV_adj=NA_real_, TargetPrev=target_prevalence
    ))
    
    cv <- .kfold_cv_run(merged, mk, k=k, repeats=repeats, seed=seed,
                        ci_method=ci_method, boot_n=boot_n, boot_stratified=boot_stratified,
                        quiet=TRUE, use_global_threshold=use_global_threshold,
                        target_prevalence=target_prevalence)
    
    if (!identical(cv$status, "ok")) return(tibble::tibble(
      Marker=mk, CV_AUC=NA_real_, CV_CI_low=NA_real_, CV_CI_high=NA_real_,
      CV_Sensitivity=NA_real_, CV_Specificity=NA_real_, CV_PPV=NA_real_, CV_NPV=NA_real_,
      K=k, Repeats=repeats, N=nrow(merged), Dataset=(dataset_name %OR% basename(expr_file)),
      Status=cv$status %OR% "degenerate", CVPlotBase=NA_character_, ThresholdMode=use_global_threshold,
      GlobalThr=cv$global_thr, GlobalDir=cv$global_dir, CV_PPV_adj=cv$ppv_adj, CV_NPV_adj=cv$npv_adj, TargetPrev=target_prevalence
    ))
    
    auc_str <- sprintf("CV AUC = %.3f (95%% CI: %.3f–%.3f)", cv$cv_auc, cv$cv_ci_low, cv$cv_ci_high)
    ttl <- if (is.null(dataset_name)) paste0("CV ROC: ", mk) else paste0("CV ROC: ", mk, " (", dataset_name, ")")
    p_cv <- .plot_roc(cv$cv_roc, title = ttl, subtitle = auc_str)
    
    safe_marker <- gsub("[^A-Za-z0-9]+", "_", mk)
    safe_ds <- gsub("[^A-Za-z0-9]+", "_", (dataset_name %OR% "dataset"))
    filepath_no_ext <- file.path(save_dir, paste0("CVROC_", safe_marker, "_", safe_ds))
    .save_plot_multi(p_cv, filepath_no_ext, formats = formats)
    
    fold_csv <- file.path(save_dir, paste0("CV_folds_", safe_marker, "_", safe_ds, ".csv"))
    readr::write_csv(cv$fold_results, fold_csv)
    
    tibble::tibble(Marker=mk, CV_AUC=cv$cv_auc, CV_CI_low=cv$cv_ci_low, CV_CI_high=cv$cv_ci_high,
                   CV_Sensitivity=cv$sens, CV_Specificity=cv$spec, CV_PPV=cv$ppv, CV_NPV=cv$npv,
                   K=k, Repeats=repeats, N=nrow(merged), Dataset=(dataset_name %OR% basename(expr_file)),
                   Status="ok", CVPlotBase=filepath_no_ext, ThresholdMode=use_global_threshold,
                   GlobalThr=cv$global_thr, GlobalDir=cv$global_dir,
                   CV_PPV_adj=cv$ppv_adj, CV_NPV_adj=cv$npv_adj, TargetPrev=target_prevalence)
  })
  
  if (isTRUE(write_results_csv)) {
    csv_name <- paste0("CV_ROC_results_", gsub("[^A-Za-z0-9]+","_", dataset_name %OR% basename(expr_file)), "_", use_global_threshold, ".csv")
    readr::write_csv(out_rows, file.path(save_dir, csv_name))
  }
  out_rows
}

# ======= Run all 6 datasets (ZNF662 + XIST + miRNAs) =======

run_all_six_roc_cv <- function(save_root = "plots_cv",
                               ci_method = "bootstrap",
                               boot_n = 400,
                               boot_stratified = TRUE,
                               formats = c("png","pdf")) {
  dir.create(save_root, showWarnings = FALSE, recursive = TRUE)
  
  # 1) ZNF662 across four datasets
  cfg_znf662 <- list(
    list(expr="GSE74530_expression.csv",  meta="GSE74530_meta.csv",  ds="GSE74530",  k=5, repeats=1, seed=42),
    list(expr="GSE150469_expression.csv", meta="GSE150469_meta.csv", ds="GSE150469", k=2, repeats=10, seed=123),
    list(expr="GSE30784_expression.csv",  meta="GSE30784_meta.csv",  ds="GSE30784",  k=5, repeats=10, seed=2025),
    list(expr="GSE246050_expression.csv", meta="GSE246050_meta.csv", ds="GSE246050", k=2, repeats=10, seed=42)
  )
  invisible(lapply(cfg_znf662, function(cfg) {
    evaluate_markers_roc_cv(
      expr_file = cfg$expr, meta_file = cfg$meta,
      markers = "ZNF662", dataset_name = cfg$ds,
      k = cfg$k, repeats = cfg$repeats, seed = cfg$seed,
      save_dir = file.path(save_root, "ZNF662"),
      formats = formats,
      ci_method = ci_method, boot_n = boot_n, boot_stratified = boot_stratified,
      write_results_csv = TRUE, use_global_threshold = "per_fold"
    )
  }))
  
  # 2) XIST (lncRNA)
  evaluate_markers_roc_cv(
    "GSE125866_expression.csv", "GSE125866_meta.csv",
    markers = "XIST", dataset_name = "GSE125866",
    k = 5, repeats = 5, seed = 7,
    save_dir = file.path(save_root, "XIST"),
    formats = formats,
    ci_method = ci_method, boot_n = boot_n, boot_stratified = boot_stratified,
    write_results_csv = TRUE, use_global_threshold = "per_fold"
  )
  
  # 3) miRNAs (batch)
  mirnas <- c("hsa-miR-625-3p","hsa-miR-625-5p","hsa-miR-28-5p",
              "hsa-miR-146b-5p","hsa-miR-424-5p","hsa-miR-374a-3p",
              "hsa-miR-16-5p","hsa-miR-15b-5p","hsa-miR-15a-5p")
  evaluate_markers_roc_cv(
    "GSE28100_expression.csv","GSE28100_meta.csv",
    markers = mirnas, dataset_name = "GSE28100",
    k = 5, repeats = 5, seed = 99,
    save_dir = file.path(save_root, "miRNA"),
    formats = formats,
    ci_method = ci_method, boot_n = boot_n, boot_stratified = boot_stratified,
    write_results_csv = TRUE, use_global_threshold = "per_fold"
  )
  
  message("All 6 dataset runs completed. Outputs in: ", normalizePath(save_root, mustWork = FALSE))
}

# Optional: run automatically when sourcing, if you want:
if (identical(Sys.getenv("RUN_EXAMPLES"), "1")) {
  run_all_six_roc_cv(save_root = "plots_cv",
                     ci_method = "bootstrap",
                     boot_n = 400,
                     boot_stratified = TRUE,
                     formats = c("png","pdf"))
}

run_all_six_roc_cv()
