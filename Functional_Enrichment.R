#!/usr/bin/env Rscript
# ============================================================
# GO / KEGG / Reactome ORA (clusterProfiler) — hardened
# Main ORA: minGSSize = 3 (statistical)
# ZNF662 panel: minGSSize = 1 (exploratory, keep even if non-sig)
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(enrichplot)
  library(DOSE)
  library(ReactomePA)
  library(msigdbr)
  library(KEGGREST)
})

# -------- 0) Paths / logging --------------------------------
wd <- getwd()
outdir <- file.path(wd, "ora_outputs")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
logf <- file.path(outdir, "ORA_log.txt")
logmsg <- function(...) {
  cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", sprintf(...), "\n"),
      file = logf, append = TRUE)
  message(sprintf(...))
}
set.seed(12345)
sink(file(logf, open = "at"))
cat("----- sessionInfo -----\n"); print(sessionInfo()); cat("\n"); sink()
logmsg("Working directory: %s", wd)

# -------- 1) Inputs ------------------------------------------
genes_symbols <- c(
  "BCAT1","SLC25A32","EDIL3","TFRC","ADAM12","RRAGD","LOXL2","LY6K","GPX3","NRG1",
  "ADORA2B","GLIPR1","ZNF662","SVIP","SLCO1B3","CXCL1","MMP13","MLPH","ZNF677","FAM221A","BEX4"
)

background_file <- "Union_Genes_All4.csv"
if (!file.exists(background_file)) stop(sprintf("Background file not found: %s", background_file))
logmsg("Background file: %s", normalizePath(background_file))

background_data <- tryCatch(
  read.csv(background_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) stop(sprintf("Failed to read background CSV: %s", e$message))
)

# Prefer a standard symbol column if present
prefer <- c("Union_Genes", "Gene", "SYMBOL", "gene", "symbol")
pick <- prefer[prefer %in% names(background_data)]
if (length(pick)) {
  symbol_col <- pick[1]
  bg_symbols <- unique(as.character(background_data[[symbol_col]]))
  logmsg("Detected background column: '%s' (%d unique entries).", symbol_col, length(bg_symbols))
} else {
  bg_symbols <- unique(as.character(unlist(background_data)))
  logmsg("No obvious symbol column; flattened all cells (%d unique entries).", length(bg_symbols))
}
bg_symbols <- bg_symbols[!is.na(bg_symbols) & nzchar(bg_symbols)]

# -------- 2) ID mapping --------------------------------------
logmsg("Converting SYMBOL -> ENTREZID (foreground + background)...")
conv_genes <- suppressWarnings(bitr(genes_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
conv_bg    <- suppressWarnings(bitr(bg_symbols,    fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))

genes_entrez <- unique(conv_genes$ENTREZID)
bg_entrez    <- unique(conv_bg$ENTREZID)

mapped_fg   <- data.frame(SYMBOL = conv_genes$SYMBOL, ENTREZID = conv_genes$ENTREZID)
unmapped_fg <- setdiff(genes_symbols, conv_genes$SYMBOL)
mapped_bg   <- data.frame(SYMBOL = conv_bg$SYMBOL, ENTREZID = conv_bg$ENTREZID)
unmapped_bg <- setdiff(bg_symbols, conv_bg$SYMBOL)

write.csv(mapped_fg, file.path(outdir, "foreground_mapped.csv"), row.names = FALSE)
write.csv(data.frame(SYMBOL = unmapped_fg), file.path(outdir, "foreground_unmapped.csv"), row.names = FALSE)
write.csv(mapped_bg, file.path(outdir, "background_mapped.csv"), row.names = FALSE)
write.csv(data.frame(SYMBOL = unmapped_bg), file.path(outdir, "background_unmapped.csv"), row.names = FALSE)

logmsg("Foreground: %d SYMBOL -> %d ENTREZ", length(genes_symbols), length(genes_entrez))
logmsg("Background: %d SYMBOL -> %d ENTREZ", length(bg_symbols), length(bg_entrez))
in_bg <- intersect(genes_entrez, bg_entrez)
logmsg("Overlap (foreground in background): %d/%d", length(in_bg), length(genes_entrez))

# -------- 3) Parameters --------------------------------------
pCutoff   <- 0.05
qCutoff   <- 0.20
pAdj      <- "BH"
minGS_main <- 3     # <<< main analysis threshold
minGS_znf  <- 1     # <<< ZNF662 panel threshold
maxGSSize <- 5000

# -------- 4) Enrichment functions ----------------------------
run_enrichGO <- function(ont, pCut=pCutoff, qCut=qCutoff, minGS=minGS_main) {
  tryCatch({
    enrichGO(
      gene          = genes_entrez,
      universe      = bg_entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pvalueCutoff  = pCut,
      qvalueCutoff  = qCut,
      pAdjustMethod = pAdj,
      minGSSize     = minGS,
      maxGSSize     = maxGSSize,
      readable      = TRUE
    )
  }, error = function(e) { logmsg("enrichGO(%s) error: %s", ont, e$message); NULL })
}

simplify_safe <- function(x) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(x)
  set.seed(12345)
  tryCatch(simplify(x, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) x)
}

# ---- MAIN ORA (minGS = 3) ----
ego_bp <- run_enrichGO("BP", minGS = minGS_main) %>% simplify_safe()
ego_mf <- run_enrichGO("MF", minGS = minGS_main) %>% simplify_safe()
ego_cc <- run_enrichGO("CC", minGS = minGS_main) %>% simplify_safe()

# KEGG (clusterProfiler flavor; no custom universe)
kegg <- tryCatch({
  enrichKEGG(
    gene          = genes_entrez,
    organism      = "hsa",
    pvalueCutoff  = pCutoff,
    pAdjustMethod = pAdj,
    minGSSize     = minGS_main,
    maxGSSize     = maxGSSize
  )
}, error = function(e) { logmsg("enrichKEGG error: %s", e$message); NULL })
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# KEGG (universe-respecting enricher flavor)
kegg_enricher <- NULL
try({
  kids <- names(keggList("pathway", "hsa"))
  t2g <- do.call(rbind, lapply(kids, function(pid) {
    gs <- keggGet(pid)[[1]]$GENE
    if (is.null(gs)) return(NULL)
    entrez <- gs[seq(1, length(gs), by = 2)]
    if (length(entrez)) data.frame(term = pid, entrez = entrez, stringsAsFactors = FALSE)
  }))
  if (!is.null(t2g) && nrow(t2g)) {
    kegg_enricher <- enricher(
      gene          = genes_entrez,
      TERM2GENE     = t2g[, c("term","entrez")],
      pAdjustMethod = pAdj,
      pvalueCutoff  = pCutoff,
      minGSSize     = minGS_main,
      universe      = bg_entrez
    )
    if (!is.null(kegg_enricher) && nrow(as.data.frame(kegg_enricher)) > 0) {
      kegg_enricher <- setReadable(kegg_enricher, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }
  }
}, silent = TRUE)

# Reactome
react <- tryCatch({
  enrichPathway(gene = genes_entrez, organism = "human",
                pvalueCutoff = pCutoff, pAdjustMethod = pAdj, minGSSize = minGS_main)
}, error = function(e) { logmsg("Reactome enrichPathway error: %s", e$message); NULL })
if (!is.null(react) && nrow(as.data.frame(react)) > 0) {
  react <- setReadable(react, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

# MSigDB (optional)
get_msig <- function(category, subcategory = NULL) msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory)
make_list <- function(msig_df) if (is.null(msig_df)) NULL else split(msig_df$entrez_gene, msig_df$gs_name)

msig_h <- tryCatch(get_msig("H"), error = function(e) NULL)
msig_c2_react <- tryCatch(get_msig("C2", "REACTOME"), error = function(e) NULL)

msig_h_list <- make_list(msig_h)
msig_react_list <- make_list(msig_c2_react)

enr_h <- tryCatch(enricher(genes_entrez, TERM2GENE = stack(msig_h_list)[, 2:1],
                           pAdjustMethod = pAdj, pvalueCutoff = pCutoff, minGSSize = minGS_main),
                  error = function(e) NULL)
enr_reactome_msig <- tryCatch(enricher(genes_entrez, TERM2GENE = stack(msig_react_list)[, 2:1],
                                       pAdjustMethod = pAdj, pvalueCutoff = pCutoff, minGSSize = minGS_main),
                              error = function(e) NULL)

# -------- 5) Save tables & plots -----------------------------
save_table <- function(res, stem) {
  outfile <- file.path(outdir, paste0(stem, ".csv"))
  if (is.null(res)) { write.csv(data.frame(), outfile, row.names = FALSE); logmsg("%s: NULL -> empty CSV.", stem); return() }
  df <- as.data.frame(res)
  if (nrow(df) > 0 && "geneID" %in% names(df)) {
    tmp <- df[, c("ID","Description","pvalue","p.adjust","qvalue","Count","GeneRatio","BgRatio","geneID")]
    tmp$geneID <- as.character(tmp$geneID)
    lst <- strsplit(tmp$geneID, "/", fixed = TRUE)
    long <- data.frame(
      ID          = rep(tmp$ID, lengths(lst)),
      Description = rep(tmp$Description, lengths(lst)),
      pvalue      = rep(tmp$pvalue, lengths(lst)),
      p_adjust    = rep(tmp$p.adjust, lengths(lst)),
      qvalue      = rep(tmp$qvalue, lengths(lst)),
      gene_symbol = unlist(lst),
      stringsAsFactors = FALSE
    )
    write.csv(long, file.path(outdir, paste0(stem, "_term_gene_table.csv")), row.names = FALSE)
    logmsg("%s: wrote expanded term–gene table (%d rows).", stem, nrow(long))
  }
  write.csv(df, outfile, row.names = FALSE)
  logmsg("%s: wrote %d rows.", stem, nrow(df))
}

pub_theme <- theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_line(linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        legend.position = "right")

save_plot <- function(plot_obj, name, width=7, height=6, dpi=300) {
  tryCatch({
    ggsave(file.path(outdir, paste0(name, ".png")), plot=plot_obj, width=width, height=height, dpi=dpi)
    ggsave(file.path(outdir, paste0(name, ".svg")), plot=plot_obj, width=width, height=height)
    logmsg("Saved %s.[png|svg]", name)
  }, error = function(e) logmsg("Failed to save %s: %s", name, e$message))
}

plot_wrap <- function(res, title_prefix) {
  if (is.null(res)) { logmsg("%s: NULL; skipping plots.", title_prefix); return() }
  df <- as.data.frame(res)
  if (nrow(df) == 0) { logmsg("%s: 0 rows; skipping plots.", title_prefix); return() }
  topN <- min(20, nrow(df))
  p_bar <- barplot(res, showCategory = topN, title = paste0(title_prefix, " — Barplot")) + pub_theme
  p_dot <- dotplot(res, showCategory = topN, title = paste0(title_prefix, " — Dotplot")) + pub_theme
  save_plot(p_bar, paste0(title_prefix, "_barplot"))
  save_plot(p_dot, paste0(title_prefix, "_dotplot"))
}

# Main saves & plots (minGS = 3)
save_table(ego_bp, "GO_BP_enrichment");  plot_wrap(ego_bp, "GO_BP")
save_table(ego_mf, "GO_MF_enrichment");  plot_wrap(ego_mf, "GO_MF")
save_table(ego_cc, "GO_CC_enrichment");  plot_wrap(ego_cc, "GO_CC")
save_table(kegg,   "KEGG_clusterProfiler");  plot_wrap(kegg, "KEGG_clusterProfiler")
if (!is.null(kegg_enricher)) { save_table(kegg_enricher, "KEGG_enricher_with_universe"); plot_wrap(kegg_enricher, "KEGG_enricher_with_universe") }
save_table(react,  "Reactome_enrichment"); plot_wrap(react, "Reactome")
if (!is.null(enr_h))             { save_table(enr_h, "MSigDB_Hallmark_enrichment");         plot_wrap(enr_h, "MSigDB_Hallmark") }
if (!is.null(enr_reactome_msig)) { save_table(enr_reactome_msig, "MSigDB_Reactome_enrichment"); plot_wrap(enr_reactome_msig, "MSigDB_Reactome") }

# -------- 6) ZNF662 panel (minGS = 1; keep even if non-sig) --
znf_filter <- function(res) {
  if (is.null(res)) return(res)
  df <- as.data.frame(res)
  if (!nrow(df)) return(df)
  if (!"geneID" %in% names(df)) return(df[0, ])
  gid  <- toupper(df$geneID)
  desc <- toupper(df$Description)
  hit <- grepl("\\bZNF662\\b", gid) | grepl("\\bZNF662\\b", desc)
  df[hit, , drop = FALSE]
}

ensure_znf <- function(res, fallback_fun) {
  out <- znf_filter(res)
  if (!is.null(out) && nrow(out) > 0) return(out)
  logmsg("ZNF662: no hits at default thresholds -> rerun permissive (p=1,q=1,minGS=1).")
  r2 <- fallback_fun()
  znf_filter(r2)
}

# Run ZNF662 panels
ego_bp_znf <- ensure_znf(ego_bp, function() run_enrichGO("BP", pCut=1, qCut=1, minGS=minGS_znf))
ego_mf_znf <- ensure_znf(ego_mf, function() run_enrichGO("MF", pCut=1, qCut=1, minGS=minGS_znf))
ego_cc_znf <- ensure_znf(ego_cc, function() run_enrichGO("CC", pCut=1, qCut=1, minGS=minGS_znf))
kegg_znf <- ensure_znf(kegg, function() {
  tmp <- tryCatch(enrichKEGG(gene = genes_entrez, organism="hsa", pvalueCutoff=1,
                             pAdjustMethod=pAdj, minGSSize=minGS_znf), error=function(e) NULL)
  if (!is.null(tmp) && nrow(as.data.frame(tmp)) > 0)
    setReadable(tmp, OrgDb=org.Hs.eg.db, keyType="ENTREZID") else tmp
})
react_znf <- ensure_znf(react, function() {
  tryCatch({ tmp <- enrichPathway(genes_entrez, "human", pvalueCutoff=1,
                                  pAdjustMethod=pAdj, minGSSize=minGS_znf)
  if (!is.null(tmp) && nrow(as.data.frame(tmp)) > 0)
    setReadable(tmp, OrgDb=org.Hs.eg.db, keyType="ENTREZID") else tmp },
  error=function(e) NULL)
})

# ---- Save + plot ZNF panels (top 10 only) ----
plot_wrap_top10 <- function(res, title_prefix) {
  if (is.null(res)) { logmsg("%s: NULL; skipping plots.", title_prefix); return() }
  df <- as.data.frame(res)
  if (nrow(df) == 0) { logmsg("%s: 0 rows; skipping plots.", title_prefix); return() }
  topN <- min(10, nrow(df))  # <<< limit to top 10
  p_bar <- barplot(res, showCategory = topN,
                   title = paste0(title_prefix, " — Top ", topN, " Barplot")) + pub_theme
  p_dot <- dotplot(res, showCategory = topN,
                   title = paste0(title_prefix, " — Top ", topN, " Dotplot")) + pub_theme
  save_plot(p_bar, paste0(title_prefix, "_barplot_top10"))
  save_plot(p_dot, paste0(title_prefix, "_dotplot_top10"))
}

# Save + plot (only top 10 terms)
save_table(ego_bp_znf, "ZNF662_GO_BP");      plot_wrap_top10(ego_bp_znf, "ZNF662_GO_BP")
save_table(ego_mf_znf, "ZNF662_GO_MF");      plot_wrap_top10(ego_mf_znf, "ZNF662_GO_MF")
save_table(ego_cc_znf, "ZNF662_GO_CC");      plot_wrap_top10(ego_cc_znf, "ZNF662_GO_CC")
save_table(kegg_znf,   "ZNF662_KEGG");       plot_wrap_top10(kegg_znf,   "ZNF662_KEGG")
save_table(react_znf,  "ZNF662_Reactome");   plot_wrap_top10(react_znf,  "ZNF662_Reactome")
