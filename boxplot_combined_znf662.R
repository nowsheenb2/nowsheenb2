# =========================
# Unified panels: mRNA (ZNF662), lncRNA (XIST), miRNAs (6 listed)
# Features: sample counts on x labels, uniform y-scale per section, optional p-values
# =========================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(cowplot)
  library(ggpubr)  # only used if show_pvals = TRUE
})

# ---- Global config ----
show_pvals <- FALSE  # <- set TRUE to show Welch t-test p-values on each panel
allowed_groups <- c("Normal", "Cancer")

# ------------------------------------------------------------
# Helper: build boxplot given a plot_df (Sample, Expression, Group), title, y_limits
# ------------------------------------------------------------
make_boxplot <- function(plot_df, panel_title, y_limits, show_pvals = FALSE) {
  # sample counts for labels
  nN <- sum(plot_df$Group == "Normal")
  nC <- sum(plot_df$Group == "Cancer")
  x_labels <- c("Normal" = paste0("Normal (n=", nN, ")"),
                "Cancer" = paste0("Cancer (n=", nC, ")"))
  
  p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.85) +
    scale_fill_manual(values = c("Normal" = "#4daf4a", "Cancer" = "#e41a1c")) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    labs(title = panel_title, x = "Group", y = "Expression") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0))
  
  if (isTRUE(show_pvals)) {
    y_lab <- y_limits[2] - 0.02 * diff(y_limits)  # consistent label position
    p <- p + stat_compare_means(method = "t.test", label = "p.format", label.y = y_lab)
  }
  p
}

# =========================
# SECTION 1 — mRNA: ZNF662 across four datasets
# =========================
mRNA_gene <- "ZNF662"
mRNA_datasets <- list(
  GSE74530  = list(expr = "GSE74530_expression.csv",  meta = "GSE74530_meta.csv"),
  GSE150469 = list(expr = "GSE150469_expression.csv", meta = "GSE150469_meta.csv"),
  GSE30784  = list(expr = "GSE30784_expression.csv",  meta = "GSE30784_meta.csv"),
  GSE246050 = list(expr = "GSE246050_expression.csv", meta = "GSE246050_meta.csv")
)

# loader that returns plot_df or NULL
load_mRNA_df <- function(ds_name, files, gene = "ZNF662") {
  expr <- read.csv(files$expr, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- read.csv(files$meta, stringsAsFactors = FALSE)
  rownames(expr) <- toupper(rownames(expr))
  g_upper <- toupper(gene)
  if (!(g_upper %in% rownames(expr))) return(NULL)
  common_samples <- intersect(colnames(expr), meta$Sample)
  if (length(common_samples) < 3) return(NULL)
  
  expr <- expr[, common_samples, drop = FALSE]
  meta <- meta %>%
    filter(Sample %in% common_samples, Group %in% allowed_groups) %>%
    mutate(Group = factor(Group, levels = c("Normal","Cancer")))
  expr <- expr[, meta$Sample, drop = FALSE]
  
  data.frame(
    Dataset = ds_name,
    Sample = meta$Sample,
    Expression = as.numeric(expr[g_upper, meta$Sample]),
    Group = meta$Group
  )
}

# gather all mRNA data and global y
mRNA_list <- lapply(names(mRNA_datasets), function(nm) load_mRNA_df(nm, mRNA_datasets[[nm]], mRNA_gene))
names(mRNA_list) <- names(mRNA_datasets)
mRNA_list <- mRNA_list[!vapply(mRNA_list, is.null, logical(1))]
stopifnot(length(mRNA_list) > 0)

mRNA_all   <- bind_rows(mRNA_list)
mRNA_range <- range(mRNA_all$Expression, na.rm = TRUE)
mRNA_pad   <- diff(mRNA_range); if (is.na(mRNA_pad) || mRNA_pad == 0) mRNA_pad <- 1
mRNA_ylim  <- c(mRNA_range[1] - 0.05*mRNA_pad, mRNA_range[2] + 0.10*mRNA_pad)

# build/save per-dataset plots
mRNA_plots <- list()
for (ds in names(mRNA_list)) {
  df <- mRNA_list[[ds]]
  title <- paste0(mRNA_gene, " Expression (", ds, ")")
  p <- make_boxplot(df, title, mRNA_ylim, show_pvals)
  mRNA_plots[[ds]] <- p
  ggsave(
    filename = paste0(ds, "_", mRNA_gene, if (show_pvals) "_boxplot_ttest" else "_boxplot_NOPVAL", ".png"),
    plot     = p, width = 7, height = 7, dpi = 300
  )
  ggsave(
    filename = paste0(ds, "_", mRNA_gene, if (show_pvals) "_boxplot_ttest" else "_boxplot_NOPVAL", ".pdf"),
    plot     = p, width = 7, height = 7
  )
}

# combined mRNA 2×2 panel
mRNA_combined <- cowplot::plot_grid(plotlist = mRNA_plots, labels = "AUTO", ncol = 2)
mRNA_base <- paste0(mRNA_gene, "_boxplot_all_mRNA_", if (show_pvals) "ttest" else "NOPVAL")
ggsave(filename = paste0(mRNA_base, ".png"), plot = mRNA_combined, width = 14, height = 14, dpi = 300)
ggsave(filename = paste0(mRNA_base, ".pdf"), plot = mRNA_combined, width = 14, height = 14)

# =========================
# SECTION 2 — miRNAs: 6 targets (GSE28100)
# =========================
miRNA_ds <- list(expr = "GSE28100_expression.csv", meta = "GSE28100_meta.csv")
mirnas <- c("hsa-miR-625-3p", "hsa-miR-424-5p", "hsa-miR-374a-3p",
            "hsa-miR-16-5p",  "hsa-miR-15b-5p", "hsa-miR-15a-5p",
            "hsa-miR-625-5p", ("hsa-miR-28-5p", "hsa-miR-146-5p")

normalize_mir <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("^hsa-", "", x)
  x <- gsub("^mirna-", "mir-", x)
  x <- gsub("_", "-", x)
  x
}

expr_mi <- read.csv(miRNA_ds$expr, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
meta_mi <- read.csv(miRNA_ds$meta, stringsAsFactors = FALSE)
common_mi <- intersect(colnames(expr_mi), meta_mi$Sample)
stopifnot(length(common_mi) > 0)
expr_mi <- expr_mi[, common_mi, drop = FALSE]
meta_mi <- meta_mi %>%
  filter(Sample %in% common_mi, Group %in% allowed_groups) %>%
  mutate(Group = factor(Group, levels = c("Normal","Cancer")))
expr_mi <- expr_mi[, meta_mi$Sample, drop = FALSE]

# Prepare miRNA plots & uniform y-limits across all selected miRNAs
mi_plots <- list()
mi_dfs   <- list()

# Find rows for each miRNA robustly
rn_norm <- normalize_mir(rownames(expr_mi))

for (mir in mirnas) {
  hits <- which(rn_norm == normalize_mir(mir))
  if (length(hits) < 1) { warning("miRNA not found: ", mir); next }
  hit <- hits[1]
  df <- data.frame(
    Dataset = "GSE28100",
    Sample = meta_mi$Sample,
    Expression = as.numeric(expr_mi[hit, meta_mi$Sample]),
    Group = meta_mi$Group,
    mir = mir
  )
  mi_dfs[[mir]] <- df
}

stopifnot(length(mi_dfs) > 0)
mi_all   <- bind_rows(mi_dfs)
mi_range <- range(mi_all$Expression, na.rm = TRUE)
mi_pad   <- diff(mi_range); if (is.na(mi_pad) || mi_pad == 0) mi_pad <- 1
mi_ylim  <- c(mi_range[1] - 0.05*mi_pad, mi_range[2] + 0.10*mi_pad)

# Build/save per-miRNA plots
for (mir in names(mi_dfs)) {
  df <- mi_dfs[[mir]]
  title <- paste0(mir, " Expression (GSE28100)")
  p <- make_boxplot(df, title, mi_ylim, show_pvals)
  mi_plots[[mir]] <- p
  safe <- gsub("[^A-Za-z0-9_\\-]", "_", mir)
  ggsave(filename = paste0("GSE28100_", safe, if (show_pvals) "_ttest" else "_NOPVAL", ".png"),
         plot = p, width = 7, height = 7, dpi = 300)
  ggsave(filename = paste0("GSE28100_", safe, if (show_pvals) "_ttest" else "_NOPVAL", ".pdf"),
         plot = p, width = 7, height = 7)
}

# combined miRNA grid
mi_combined <- cowplot::plot_grid(plotlist = mi_plots, labels = "AUTO", ncol = 3)
mi_base <- paste0("miRNA_boxplots_GSE28100_", if (show_pvals) "ttest" else "NOPVAL")
ggsave(filename = paste0(mi_base, ".png"), plot = mi_combined, width = 14, height = 12, dpi = 300)
ggsave(filename = paste0(mi_base, ".pdf"), plot = mi_combined, width = 14, height = 12)

# =========================
# SECTION 3 — lncRNA: XIST (GSE125866)
# =========================
lnc_ds <- list(expr = "GSE125866_expression.csv", meta = "GSE125866_meta.csv")
lnc_target <- "XIST"

expr_lnc <- read.csv(lnc_ds$expr, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
meta_lnc <- read.csv(lnc_ds$meta, stringsAsFactors = FALSE)

common_lnc <- intersect(colnames(expr_lnc), meta_lnc$Sample)
stopifnot(length(common_lnc) > 0)
expr_lnc <- expr_lnc[, common_lnc, drop = FALSE]
meta_lnc <- meta_lnc %>%
  filter(Sample %in% common_lnc, Group %in% allowed_groups) %>%
  mutate(Group = factor(Group, levels = c("Normal","Cancer")))
expr_lnc <- expr_lnc[, meta_lnc$Sample, drop = FALSE]

rownames(expr_lnc) <- toupper(trimws(rownames(expr_lnc)))
if (!toupper(lnc_target) %in% rownames(expr_lnc)) stop("XIST not found in lncRNA file.")

lnc_df <- data.frame(
  Dataset = "GSE125866",
  Sample = meta_lnc$Sample,
  Expression = as.numeric(expr_lnc[toupper(lnc_target), meta_lnc$Sample]),
  Group = meta_lnc$Group
)

lnc_range <- range(lnc_df$Expression, na.rm = TRUE)
lnc_pad   <- diff(lnc_range); if (is.na(lnc_pad) || lnc_pad == 0) lnc_pad <- 1
lnc_ylim  <- c(lnc_range[1] - 0.05*lnc_pad, lnc_range[2] + 0.10*lnc_pad)

lnc_title <- paste0(lnc_target, " Expression (GSE125866)")
lnc_plot  <- make_boxplot(lnc_df, lnc_title, lnc_ylim, show_pvals)
ggsave(filename = paste0("GSE125866_", lnc_target, if (show_pvals) "_ttest" else "_NOPVAL", ".png"),
       plot = lnc_plot, width = 7, height = 7, dpi = 300)
ggsave(filename = paste0("GSE125866_", lnc_target, if (show_pvals) "_ttest" else "_NOPVAL", ".pdf"),
       plot = lnc_plot, width = 7, height = 7)

# =========================
# GRAND COMBINED FIGURE — mRNA (4) + miRNA (6) + lncRNA (1) = 11 panels
# =========================
grand_plots <- c(mRNA_plots, mi_plots, list(lnc_plot))
grand_combined <- cowplot::plot_grid(plotlist = grand_plots, labels = "AUTO", ncol = 3)
grand_base <- paste0("ALL_panels_mRNA_miRNA_lncRNA_", if (show_pvals) "ttest" else "NOPVAL")
ggsave(filename = paste0(grand_base, ".png"), plot = grand_combined, width = 21, height = 21, dpi = 300)
ggsave(filename = paste0(grand_base, ".pdf"), plot = grand_combined, width = 21, height = 21)

cat("✅ Done. Saved:\n",
    "- Per-dataset ZNF662 panels + combined mRNA panel\n",
    "- 6 miRNA panels + combined miRNA panel\n",
    "- XIST (GSE125866)\n",
    "- GRAND combined figure with all panels\n")
