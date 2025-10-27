suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(forestmodel)
  library(broom)
  library(cowplot)
})

# =========================
# I/O (update paths)
# =========================
expr_file <- "expression_filtered.csv"
clin_file <- "clinical_filtered.csv"

dir.create("KM_Plots_MED",    showWarnings = FALSE)
dir.create("Cox_Results_MED", showWarnings = FALSE)
dir.create("Cox_Plots_MED",   showWarnings = FALSE)

# =========================
# Config (MEDIAN SPLIT)
# =========================
min_group_n   <- 10       # require >= this many in each arm or skip gene
km_show_pval  <- TRUE     # show p-value on KM panels
sig_p_cutoff  <- 0.05     # highlight if log-rank p < this
title_suffix  <- "Survival (Median Split: High vs Low)"

# =========================
# Load & harmonize
# =========================
expr_df <- read_csv(expr_file, show_col_types = FALSE)
colnames(expr_df)[1] <- "sampleID"
expr_df <- expr_df %>% mutate(sampleID = str_to_upper(str_sub(sampleID, 1, 15)))

clin_df_raw <- read_csv(clin_file, show_col_types = FALSE)
clin_df <- clin_df_raw %>%
  mutate(
    sampleID       = str_to_upper(str_sub(sampleID, 1, 15)),
    days_to_death  = suppressWarnings(as.numeric(days_to_death)),
    days_to_fu     = suppressWarnings(as.numeric(days_to_last_followup)),
    vital_status   = ifelse(str_to_upper(vital_status) == "DECEASED", 1, 0),
    survival_time  = ifelse(!is.na(days_to_death), days_to_death, days_to_fu)
  ) %>%
  select(sampleID, survival_time, vital_status) %>%
  filter(!is.na(sampleID), !is.na(survival_time), !is.na(vital_status),
         survival_time > 0)

# de-dup
expr_df <- expr_df %>% distinct(sampleID, .keep_all = TRUE)
clin_df <- clin_df %>% distinct(sampleID, .keep_all = TRUE)

merged_df_all <- inner_join(expr_df, clin_df, by = "sampleID")
cat("Merged dataset dimensions: ", paste(dim(merged_df_all), collapse = " x "), "\n")

genes <- setdiff(colnames(expr_df), "sampleID")
stopifnot(length(genes) > 0)

# =========================
# Helpers
# =========================
median_split <- function(vec) {
  med <- median(vec, na.rm = TRUE)
  grp <- ifelse(vec >= med, "High", "Low")
  tabs <- table(grp)
  if (length(tabs) < 2 || any(tabs == 0)) {
    ord <- order(vec, na.last = NA)
    n <- length(ord); cut <- floor(n/2)
    grp <- rep(NA_character_, n)
    grp[ord[seq_len(cut)]]     <- "Low"
    grp[ord[seq(cut+1, n)]]    <- "High"
  }
  factor(grp, levels = c("Low","High"))
}

# Consistent KM plotter (returns list(plot=ggplot, p_logrank=numeric))
make_km_panel <- function(df, gene, show_pval = TRUE, sig_p = 0.05) {
  # log-rank p (for highlight)
  sd_obj <- survdiff(Surv(survival_time, vital_status) ~ group, data = df)
  p_lr   <- 1 - pchisq(sd_obj$chisq, df = length(sd_obj$n) - 1)
  
  # p-value placement: near top-left of each panel
  px <- 0.05 * max(df$survival_time, na.rm = TRUE)
  py <- 0.95
  
  fit_km <- survfit(Surv(survival_time, vital_status) ~ group, data = df)
  p_km   <- ggsurvplot(
    fit_km, data = df,
    pval = show_pval,
    pval.coord = c(px, py),
    risk.table = TRUE,
    conf.int = FALSE,
    title = paste0(gene, " - ", title_suffix),
    legend.title = "Expression",
    palette = c("dodgerblue3", "firebrick3"),
    risk.table.y.text.col = TRUE,
    risk.table.height = 0.22
  )
  
  # highlight significant panels: red border + star
  if (!is.na(p_lr) && p_lr < sig_p) {
    panel <- ggdraw(p_km$plot) +
      draw_plot_label("★", x = 0.02, y = 0.98, size = 18) +
      draw_rect(color = "firebrick3", size = 2, fill = NA)
  } else {
    panel <- p_km$plot
  }
  
  list(plot = panel, p_logrank = p_lr)
}

cox_per_sd <- function(x, surv_time, status) {
  z <- scale(x)
  m <- coxph(Surv(surv_time, status) ~ z)
  broom::tidy(m) %>%
    mutate(HR = exp(estimate),
           HR_low = exp(estimate - 1.96*std.error),
           HR_up  = exp(estimate + 1.96*std.error),
           model = "Per 1 SD")
}

cox_group <- function(group, surv_time, status) {
  m <- coxph(Surv(survival_time, vital_status) ~ group,
             data = data.frame(survival_time = surv_time, vital_status = status, group = group))
  broom::tidy(m) %>% filter(term == "groupHigh") %>%
    mutate(HR = exp(estimate),
           HR_low = exp(estimate - 1.96*std.error),
           HR_up  = exp(estimate + 1.96*std.error),
           model = "Median split")
}

# =========================
# Main loop
# =========================
km_panels <- list()
cox_rows  <- list()

for (gene in genes) {
  cat("Analyzing gene:", gene, "\n")
  df <- merged_df_all %>% filter(!is.na(.data[[gene]]))
  if (nrow(df) < (2*min_group_n)) {
    cat("  Skipped: too few with expression\n"); next
  }
  
  df <- df %>% mutate(group = median_split(.data[[gene]]))
  tabs <- table(df$group)
  if (any(tabs < min_group_n)) {
    cat("  Skipped: insufficient group sizes after median split\n\n"); next
  }
  if (all(df$vital_status == 0)) {
    cat("  Skipped: no events\n\n"); next
  }
  
  # KM panel (with auto-highlight)
  km <- make_km_panel(df, gene, show_pval = km_show_pval, sig_p = sig_p_cutoff)
  
  # Save individual KM
  ggsave(file.path("KM_Plots_MED", paste0(gene, "_KM_MED.png")),
         plot = km$plot, width = 7, height = 7, dpi = 300)
  ggsave(file.path("KM_Plots_MED", paste0(gene, "_KM_MED.pdf")),
         plot = km$plot, width = 7, height = 7, device = cairo_pdf)
  
  km_panels[[gene]] <- km$plot
  
  # Cox models
  row_group <- cox_group(df$group, df$survival_time, df$vital_status) %>% mutate(gene = gene)
  row_cont  <- cox_per_sd(df[[gene]], df$survival_time, df$vital_status) %>% mutate(gene = gene)
  cox_rows[[gene]] <- bind_rows(row_group, row_cont)
  
  # Per-gene forest (optional)
  g_forest <- forest_model(coxph(Surv(survival_time, vital_status) ~ group, data = df))
  ggsave(file.path("Cox_Plots_MED", paste0(gene, "_Cox_forest_MED.png")),
         plot = g_forest, width = 6, height = 4, dpi = 300)
  
  # Save raw Cox summary text
  sink(file.path("Cox_Results_MED", paste0(gene, "_Cox_median.txt"))); 
  print(summary(coxph(Surv(survival_time, vital_status) ~ group, data = df)))
  sink()
  
  cat("  Done.\n\n")
}

# =========================
# Combined forest plot
# =========================
stopifnot(length(cox_rows) > 0)

forest_df <- bind_rows(cox_rows) %>%
  mutate(
    p_adj_BH  = p.adjust(p.value, method = "BH"),
    dir       = ifelse(HR >= 1, "Higher risk (HR>1)", "Lower risk (HR<1)"),
    label_hr  = sprintf("%.2f (%.2f–%.2f)", HR, HR_low, HR_up),
    label_p   = ifelse(p.value < 0.001, "p<0.001", paste0("p=", sprintf("%.3f", p.value))),
    label_all = paste0("HR ", label_hr, "; ", label_p),
    row_lab   = paste(gene, model)
  ) %>%
  arrange(model, HR) %>%
  mutate(row_lab = factor(row_lab, levels = unique(row_lab)))

pad  <- 0.06 * diff(range(c(forest_df$HR_low, forest_df$HR_up), na.rm = TRUE))
laby <- forest_df$HR_up + pad

combined_forest_plot <- ggplot(forest_df, aes(x = row_lab, y = HR, color = dir)) +
  geom_errorbar(aes(ymin = HR_low, ymax = HR_up), width = 0.15, linewidth = 0.55) +
  geom_point(shape = 16, size = 2.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(aes(y = laby, label = label_all), size = 3, hjust = 0, color = "black") +
  coord_flip(clip = "off") +
  scale_color_manual(values = c("Lower risk (HR<1)" = "steelblue4",
                                "Higher risk (HR>1)" = "firebrick3"),
                     name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Combined Cox Proportional Hazards Forest Plot",
    x = "Gene (Model)",
    y = "Hazard Ratio"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text  = element_text(size = 9),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 28, 10, 10)
  )

ggsave(file.path("Cox_Plots_MED", "Combined_Cox_Forest_median.png"),
       plot = combined_forest_plot, width = 8.5, height = 5.2, dpi = 300)
ggsave(file.path("Cox_Plots_MED", "Combined_Cox_Forest_median.pdf"),
       plot = combined_forest_plot, width = 8.5, height = 5.2, device = cairo_pdf)

readr::write_csv(forest_df, file.path("Cox_Results_MED", "Combined_Cox_Tidy_Summary_median.csv"))

# =========================
# ONE-PAGE COMPOSITE: all KM (with highlights) + combined forest
# =========================
stopifnot(length(km_panels) > 0)
n_panels  <- length(km_panels)
ncol_grid <- if (n_panels <= 2) 2 else if (n_panels <= 4) 2 else 3

km_grid <- cowplot::plot_grid(plotlist = km_panels, labels = "AUTO", ncol = ncol_grid)

composite <- cowplot::plot_grid(
  km_grid,
  combined_forest_plot,
  ncol = 1,
  rel_heights = c(2.0, 1.0)
)

ggsave("ALL_Survival_and_Cox_MEDIAN_COMPOSITE.png", composite, width = 12, height = 16, dpi = 300)
ggsave("ALL_Survival_and_Cox_MEDIAN_COMPOSITE.pdf", composite, width = 12, height = 16)

cat("\n✅ Done.\nHighlighted KM panels where log-rank p <", sig_p_cutoff, ".\n",
    "Saved per-gene KMs, combined forest, and one-page composite.\n")
