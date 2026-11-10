# ============================================================
# 01d_ssGSEA_plots.R
#
# READ-AND-PLOT ONLY: loads the scored-data bundle written by
# 01c_ssGSEA_compute.R and produces:
#   - overview + publication pathway-activity panels (score_sd_all)
#   - limma statistics on RAW ssGSEA scores (17vs0 / 68vs0 / 68vs17)
#   - ctrl-referenced panel with significance annotations (score_sd_ctrl)
#   - diagnostic: ssGSEA score vs mean pathway expression
#
# Input : processed/bulkRNAseq/ssGSEA/ssGSEA_scored_data.rds
# Output: output/bulkRNAseq/tables/ssGSEA/...
#         output/bulkRNAseq/figures/ssGSEA/...
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(grid)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
processed_dir <- file.path("processed", "bulkRNAseq", "ssGSEA")
tables_dir    <- file.path("output", "bulkRNAseq", "tables", "ssGSEA")
figures_dir   <- file.path("output", "bulkRNAseq", "figures", "ssGSEA")
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

bundle_path <- file.path(processed_dir, "ssGSEA_scored_data.rds")
stopifnot(file.exists(bundle_path))
bundle <- readRDS(bundle_path)

df_long       <- bundle$df_long
meta          <- bundle$meta
cond_cols     <- bundle$cond_cols
gene_sets_use <- bundle$gene_sets_use
expr_for_gsva <- bundle$expr_for_gsva
gsva_mat      <- bundle$gsva_mat

# ------------------------------------------------------------
# write the raw + standardized pathway score table
# (computed in 01c; 01c itself never writes to output/, only processed/)
# ------------------------------------------------------------
readr::write_csv(
  df_long |>
    dplyr::select(sample, condition_nm, pathway, pathway_label,
                  ssgsea_score, score_sd_all, ctrl_mean, ctrl_sd, score_sd_ctrl),
  file.path(tables_dir, "PATHWAY_SCORES_GSVA_ssGSEA_raw_allSD_and_ctrlSD.csv")
)

# ------------------------------------------------------------
# small helpers (kept local to this script; pure + cheap)
# ------------------------------------------------------------
std_limma_cols <- function(df) {
  df <- as.data.frame(df)
  n0 <- colnames(df)
  n1 <- tolower(n0)
  n1 <- gsub("[^a-z0-9]+", "_", n1)
  n1 <- gsub("^_+|_+$", "", n1)
  n1 <- sub("^adj_p_val$", "fdr", n1)
  n1 <- sub("^logfc$", "logFC", n1)
  colnames(df) <- n1
  df
}

p_to_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

make_pathway_label <- function(x, wrap_width = 32) {
  x2 <- x |>
    stringr::str_replace("\\s*\\(n=\\d+\\)\\s*$", "") |>
    stringr::str_replace("^HALLMARK__", "HM: ") |>
    stringr::str_replace("^GOBP__", "BP: ") |>
    stringr::str_replace("^GOMF__", "MF: ") |>
    stringr::str_replace("^KEGG__", "KEGG: ") |>
    stringr::str_replace_all("__", ": ")
  stringr::str_wrap(x2, width = wrap_width)
}

# ------------------------------------------------------------
# 10) overview panel using score_sd_all
#
# IMPORTANT: score_sd_all is a pathway-wise z-score across ALL
# samples, i.e. NOT ctrl-referenced.
# ------------------------------------------------------------
p_overview <- ggplot(df_long, aes(x = condition_nm, y = score_sd_all)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = condition_nm),
             position = position_jitter(width = 0.12),
             size = 2) +
  scale_color_manual(values = cond_cols) +
  facet_wrap(~ pathway_label, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 9)
  ) +
  labs(y = "Pathway activity (ssGSEA; z-score across all samples)")

ggsave(file.path(figures_dir, "PATHWAY_SCORE_PANEL_GSVA_allSampleSD.pdf"),
       p_overview, width = 11, height = 8)

# ------------------------------------------------------------
# 11) publication panel using score_sd_all (selected 8-9 pathways)
# ------------------------------------------------------------
main_order <- c(
  "HM: OXPHOS",
  "HM: E2F_TARGETS",
  "HM: CHOLESTEROL_HOMEOSTASIS",
  "HM: MYC",
  "HM: DNA_REPAIR",
  "HM: G2M",
  "HM: UNFOLDED"
)

ferro_label <- "KEGG: hsa04216: Ferroptosis"
iron_label  <- "MF: GO:0005506: iron.ion.binding"
order8 <- c(main_order, ferro_label, iron_label)

avail <- sort(unique(df_long$pathway_label))

cat("\n--- AVAILABLE pathway labels ---\n")
print(avail)

cat("\n--- EXPECTED selected labels ---\n")
print(order8)

cat("\n--- MISSING selected labels ---\n")
print(setdiff(order8, avail))

cat("\n--- EXTRA available labels ---\n")
print(setdiff(avail, order8))

order8 <- order8[order8 %in% unique(df_long$pathway_label)]
stopifnot(length(order8) == 9)

df_plot <- df_long |>
  dplyr::filter(pathway_label %in% order8) |>
  dplyr::mutate(
    pathway_label = factor(pathway_label, levels = order8),
    condition_nm  = factor(condition_nm, levels = c("0nM", "17nM", "68nM"))
  )

use_fixed_y <- TRUE
use_ylim <- TRUE
ylim_vals <- c(-3, 3)

p_all <- ggplot(df_plot, aes(x = condition_nm, y = score_sd_all)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  geom_point(aes(color = condition_nm),
             position = position_jitter(width = 0.12),
             size = 2.0, alpha = 1) +
  scale_color_manual(values = cond_cols) +
  facet_wrap(
    ~ pathway_label,
    ncol = 2,
    scales = if (use_fixed_y) "fixed" else "free_y"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  ) +
  labs(y = "Pathway activity (ssGSEA; z-score across all samples)")

if (use_ylim) {
  p_all <- p_all + coord_cartesian(ylim = ylim_vals)
}

ggsave(file.path(figures_dir, "PATHWAY_PANEL_ALL8_2col_allSampleSD.pdf"),
       p_all, width = 8, height = 17, units = "in", limitsize = FALSE)

base_size_pts <- 8
p_all_pub <- p_all +
  theme_bw(base_size = base_size_pts) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = base_size_pts + 8),
    axis.text.x  = element_text(size = base_size_pts + 4),
    axis.text.y  = element_text(size = base_size_pts + 4),
    strip.text   = element_text(size = base_size_pts + 4),
    panel.spacing = unit(0.6, "lines")
  )

fig_width <- 6
fig_height <- fig_width * 2

ggsave(file.path(figures_dir, "PATHWAY_PANEL_ALL8_publication_allSampleSD.pdf"),
       p_all_pub, width = fig_width, height = fig_height, units = "in",
       device = grDevices::cairo_pdf)

ggsave(file.path(figures_dir, "PATHWAY_PANEL_ALL8_publication_allSampleSD.svg"),
       p_all_pub, width = fig_width, height = fig_height, units = "in",
       device = svglite::svglite)

# ------------------------------------------------------------
# 12) test pathway score differences with limma on RAW ssGSEA scores
#
# IMPORTANT: limma is fit on RAW ssGSEA scores (ssgsea_score),
# not on score_sd_all and not on score_sd_ctrl.
#
# F_anyDose/P_anyDose/FDR_anyDose = omnibus test: does at least one of
# the three dose contrasts (17vs0, 68vs0, 68vs17) differ from zero?
# This replaces an earlier version of this script that ran the F-test
# on the raw group-mean fit instead of the contrasts fit -- that tested
# "are the group means all zero?" (nearly always true for real scores)
# rather than "do the groups differ?", which made it misleadingly
# significant regardless of any real treatment effect.
# ------------------------------------------------------------
Y <- gsva_mat
meta_fit <- meta[colnames(Y), , drop = FALSE]
stopifnot(all(colnames(Y) == rownames(meta_fit)))

meta_fit$condition_nm <- factor(meta_fit$condition_nm, levels = c("0nM", "17nM", "68nM"))
meta_fit$condition_safe <- factor(
  meta_fit$condition_nm,
  levels = c("0nM", "17nM", "68nM"),
  labels = c("nM0", "nM17", "nM68")
)

design <- model.matrix(~0 + condition_safe, data = meta_fit)
colnames(design) <- levels(meta_fit$condition_safe)

fit <- limma::lmFit(Y, design)
cont <- limma::makeContrasts(
  `17vs0`  = nM17 - nM0,
  `68vs0`  = nM68 - nM0,
  `68vs17` = nM68 - nM17,
  levels = design
)

fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))

# NOTE: the omnibus/overall test below must be run on the CONTRASTS fit
# (fit2), not on the raw group-mean fit. `design` has no intercept
# (~0 + condition_safe), so each raw coefficient is a group MEAN, not a
# difference -- an F-test on those (coef = NULL on the raw fit) asks
# "are the group means all zero?", which is nearly always true for any
# real score, not "do the groups differ?" and would be misleadingly
# significant. Testing coef = 1:3 on fit2 correctly asks "does at least
# one of the three condition contrasts differ from zero?"
tt_F <- limma::topTable(fit2, coef = 1:3, number = Inf, sort.by = "F") |>
  tibble::rownames_to_column("pathway") |>
  std_limma_cols()

tt_17vs0 <- limma::topTable(fit2, coef = "17vs0", number = Inf, sort.by = "P") |>
  tibble::rownames_to_column("pathway") |>
  std_limma_cols()

tt_68vs0 <- limma::topTable(fit2, coef = "68vs0", number = Inf, sort.by = "P") |>
  tibble::rownames_to_column("pathway") |>
  std_limma_cols()

tt_68vs17 <- limma::topTable(fit2, coef = "68vs17", number = Inf, sort.by = "P") |>
  tibble::rownames_to_column("pathway") |>
  std_limma_cols()

stats_tbl <- tt_F |>
  dplyr::select(pathway, f, p_value, fdr) |>
  dplyr::rename(F_anyDose = f, P_anyDose = p_value, FDR_anyDose = fdr) |>
  dplyr::left_join(
    tt_17vs0 |>
      dplyr::select(pathway, logFC, p_value, fdr) |>
      dplyr::rename(logFC_17vs0 = logFC, P_17vs0 = p_value, FDR_17vs0 = fdr),
    by = "pathway"
  ) |>
  dplyr::left_join(
    tt_68vs0 |>
      dplyr::select(pathway, logFC, p_value, fdr) |>
      dplyr::rename(logFC_68vs0 = logFC, P_68vs0 = p_value, FDR_68vs0 = fdr),
    by = "pathway"
  ) |>
  dplyr::left_join(
    tt_68vs17 |>
      dplyr::select(pathway, logFC, p_value, fdr) |>
      dplyr::rename(logFC_68vs17 = logFC, P_68vs17 = p_value, FDR_68vs17 = fdr),
    by = "pathway"
  )

readr::write_csv(stats_tbl, file.path(tables_dir, "PATHWAY_STATS_limma_rawssGSEA.csv"))

# ------------------------------------------------------------
# 13) ctrl-referenced panel with significance annotations
#
# IMPORTANT: this plot uses score_sd_ctrl:
# (sample score - ctrl mean) / ctrl SD, so all y-values are
# explicitly relative to the ctrl group.
# ------------------------------------------------------------
ann <- stats_tbl |>
  dplyr::mutate(
    star_17vs0 = p_to_stars(FDR_17vs0),
    star_68vs0 = p_to_stars(FDR_68vs0)
  ) |>
  dplyr::select(pathway, star_17vs0, star_68vs0)

ypos <- df_long |>
  dplyr::group_by(pathway) |>
  dplyr::summarise(
    y_max = max(score_sd_ctrl, na.rm = TRUE),
    y_min = min(score_sd_ctrl, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    span = pmax(1, y_max - y_min),
    y1 = y_max + 0.10 * span,
    y2 = y_max + 0.22 * span
  )

ann_plot <- ypos |>
  dplyr::left_join(ann, by = "pathway") |>
  dplyr::left_join(
    df_long |>
      dplyr::distinct(pathway, pathway_label),
    by = "pathway"
  ) |>
  dplyr::mutate(
    x1 = "0nM",
    x2 = "17nM",
    x3 = "68nM"
  )

p_ctrl <- ggplot(df_long, aes(x = condition_nm, y = score_sd_ctrl)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = condition_nm),
             position = position_jitter(width = 0.12),
             size = 2) +
  scale_color_manual(values = cond_cols) +
  facet_wrap(~ pathway_label, scales = "fixed", ncol = 2) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 9)
  ) +
  labs(y = "Pathway activity (ssGSEA; SD relative to ctrl)")

p_ctrl <- p_ctrl +
  geom_segment(
    data = ann_plot,
    aes(x = x1, xend = x2, y = y1, yend = y1),
    inherit.aes = FALSE, linewidth = 0.3
  ) +
  geom_segment(
    data = ann_plot,
    aes(x = x1, xend = x1, y = y1, yend = y1 - 0.03 * span),
    inherit.aes = FALSE, linewidth = 0.3
  ) +
  geom_segment(
    data = ann_plot,
    aes(x = x2, xend = x2, y = y1, yend = y1 - 0.03 * span),
    inherit.aes = FALSE, linewidth = 0.3
  ) +
  geom_text(
    data = ann_plot,
    aes(x = 1.5, y = y1 + 0.03 * span, label = star_17vs0),
    inherit.aes = FALSE, size = 3
  ) +
  geom_segment(
    data = ann_plot,
    aes(x = x1, xend = x3, y = y2, yend = y2),
    inherit.aes = FALSE, linewidth = 0.3
  ) +
  geom_segment(
    data = ann_plot,
    aes(x = x1, xend = x1, y = y2, yend = y2 - 0.03 * span),
    inherit.aes = FALSE, linewidth = 0.3
  ) +
  geom_segment(
    data = ann_plot,
    aes(x = x3, xend = x3, y = y2, yend = y2 - 0.03 * span),
    inherit.aes = FALSE, linewidth = 0.3
  ) +
  geom_text(
    data = ann_plot,
    aes(x = 2.0, y = y2 + 0.03 * span, label = star_68vs0),
    inherit.aes = FALSE, size = 3
  )

ggsave(file.path(figures_dir, "PATHWAY_SCORE_PANEL_GSVA_ctrlSD_withStats_limmaFDR.pdf"),
       p_ctrl, width = 11, height = 8)

# ------------------------------------------------------------
# 14) diagnostic: ssGSEA score vs mean pathway expression
#
# Compares per-sample mean pathway expression (same genes that
# entered ssGSEA) to raw ssGSEA scores.
#   ssgsea_z:    z-score of raw ssGSEA across all samples within pathway
#   mean_expr_z: z-score of mean pathway expression across all samples
# ------------------------------------------------------------
scores_raw_diag <- as.data.frame(t(gsva_mat))
scores_raw_diag$sample <- rownames(scores_raw_diag)
scores_raw_diag$condition_nm <- meta[scores_raw_diag$sample, "condition_nm"]

ssgsea_long_diag <- scores_raw_diag |>
  tidyr::pivot_longer(
    cols = -c(sample, condition_nm),
    names_to = "pathway",
    values_to = "ssgsea_score"
  )

mean_expr_long <- lapply(names(gene_sets_use), function(pw) {
  genes <- gene_sets_use[[pw]]
  submat <- expr_for_gsva[genes, , drop = FALSE]
  
  data.frame(
    sample = colnames(submat),
    pathway = pw,
    pathway_mean_expr = colMeans(submat, na.rm = TRUE),
    pathway_median_expr = apply(submat, 2, median, na.rm = TRUE),
    n_genes_used = nrow(submat),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

diag_long <- ssgsea_long_diag |>
  dplyr::left_join(mean_expr_long, by = c("sample", "pathway")) |>
  dplyr::mutate(
    pathway_label = make_pathway_label(pathway, wrap_width = 34),
    condition_nm = factor(condition_nm, levels = c("0nM", "17nM", "68nM"))
  ) |>
  dplyr::group_by(pathway) |>
  dplyr::mutate(
    ssgsea_z = as.numeric(scale(ssgsea_score)),
    mean_expr_z = as.numeric(scale(pathway_mean_expr))
  ) |>
  dplyr::ungroup()

readr::write_csv(
  diag_long |>
    dplyr::select(
      sample, condition_nm, pathway, pathway_label,
      ssgsea_score, ssgsea_z,
      pathway_mean_expr, pathway_median_expr, mean_expr_z,
      n_genes_used
    ),
  file.path(tables_dir, "DIAGNOSTIC_ssGSEA_vs_meanExpression_perSample.csv")
)

cor_tbl <- diag_long |>
  dplyr::group_by(pathway, pathway_label) |>
  dplyr::summarise(
    n_samples = sum(complete.cases(ssgsea_score, pathway_mean_expr)),
    pearson_r = suppressWarnings(cor(
      ssgsea_score, pathway_mean_expr,
      use = "pairwise.complete.obs",
      method = "pearson"
    )),
    spearman_rho = suppressWarnings(cor(
      ssgsea_score, pathway_mean_expr,
      use = "pairwise.complete.obs",
      method = "spearman"
    )),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(spearman_rho))

readr::write_csv(cor_tbl, file.path(tables_dir, "DIAGNOSTIC_ssGSEA_vs_meanExpression_correlations.csv"))

p_scatter <- ggplot(
  diag_long,
  aes(x = pathway_mean_expr, y = ssgsea_score, color = condition_nm)
) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
  scale_color_manual(values = cond_cols) +
  facet_wrap(~ pathway_label, scales = "free", ncol = 2) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 9),
    axis.text = element_text(size = 8)
  ) +
  labs(
    x = "Mean pathway gene expression",
    y = "Raw ssGSEA score",
    color = "Condition",
    title = "Per-sample ssGSEA score vs mean pathway expression"
  )

ggsave(file.path(figures_dir, "DIAGNOSTIC_ssGSEA_vs_meanExpression_scatter.pdf"),
       p_scatter, width = 11, height = 8.5)

p_mean_box <- ggplot(
  diag_long,
  aes(x = condition_nm, y = pathway_mean_expr, color = condition_nm)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), size = 2) +
  scale_color_manual(values = cond_cols) +
  facet_wrap(~ pathway_label, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 9)
  ) +
  labs(
    x = NULL,
    y = "Mean pathway gene expression",
    title = "Mean pathway expression across conditions"
  )

ggsave(file.path(figures_dir, "DIAGNOSTIC_meanExpression_byCondition.pdf"),
       p_mean_box, width = 11, height = 8.5)

diag_long_z <- diag_long |>
  dplyr::select(sample, condition_nm, pathway, pathway_label, ssgsea_z, mean_expr_z) |>
  tidyr::pivot_longer(
    cols = c(ssgsea_z, mean_expr_z),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      ssgsea_z = "ssGSEA",
      mean_expr_z = "Mean expression"
    )
  )

p_compare <- ggplot(
  diag_long_z,
  aes(x = condition_nm, y = value, color = condition_nm)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), size = 1.8) +
  scale_color_manual(values = cond_cols) +
  facet_grid(metric ~ pathway_label, scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 8)
  ) +
  labs(
    x = NULL,
    y = "Within-pathway z-score across all samples",
    title = "ssGSEA vs mean pathway expression"
  )

ggsave(file.path(figures_dir, "DIAGNOSTIC_ssGSEA_vs_meanExpression_zscoreCompare.pdf"),
       p_compare, width = 14, height = 8.5, limitsize = FALSE)

message("Done.")
message("Tables written under:  ", tables_dir)
message("Figures written under: ", figures_dir)
