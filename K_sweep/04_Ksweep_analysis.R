# ============================================================
# K-sensitivity of the sparse NMF decomposition
#
# Tests the choice of K = 30 (L1 = 0.5) across K = 15, 20, 25, 30, 35, 40
# in two ways: in-sample reconstruction error, and correspondence with
# RCTD-derived cell-type fractions. For each K the NMF factor most
# positively correlated with the RCTD-derived Schwann-cell fraction and
# the one most positively correlated with the malignant-cell fraction are
# identified, together with the summed activity of the top three
# positively correlated factors. Cosine similarity of gene loadings then
# asks whether the selected factor recovers a similar programme at
# different K.
#
# RCTD is independent of the NMF factorisation but is estimated from the
# same spatial expression data, so this is a reference-based consistency
# assessment rather than independent biological validation.
#
# Inputs (input/ksweep/): one NMF fit and one recon_error_k*.csv per K,
#   produced by run_hvg_matrix_container.sh and run_ksweep_container.sh
#   on the cluster; copy them here keeping the job's filenames.
#   input/RCTD_weights.Rds: RCTD cell-type weights per spot.
#
# Outputs (output/04_K_sensitivity/):
#   tables/NMF_K_reconstruction_error.csv  per-K errors combined
#   tables/K_sensitivity_summary.csv       main table: selected factors,
#                                          correlations, error per K
#   tables/top_factor_redundancy.csv       cosine similarity within K
#   tables/crossK_cosine_Nerve.csv / _Tumour.csv
#   figures/K_sensitivity_panel.pdf        supplementary figure
#   figures/heatmaps/crossK_cosine_Nerve.pdf / _Tumour.pdf
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Package 'pheatmap' is required: install.packages(\"pheatmap\")")
}
library(pheatmap)

# ---- 0. Directory structure -------------------------------------------
# Run from the project folder (the one containing input/), or from scripts/.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input")) stop("Run this script from the project folder (the one containing input/).")

INPUT_DIR     <- "input"

ANALYSIS_NAME <- "04_K_sensitivity"
TABLES_DIR    <- file.path("output", ANALYSIS_NAME, "tables")
FIGURES_DIR   <- file.path("output", ANALYSIS_NAME, "figures")
HEATMAPS_DIR  <- file.path(FIGURES_DIR, "heatmaps")
dir.create(TABLES_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(HEATMAPS_DIR, recursive = TRUE, showWarnings = FALSE)

# The cluster K-sweep outputs are produced elsewhere -- copy them into
# input/ksweep/ before running, keeping the cluster job's filenames
# (PNI_NMF_L105-k*.rds and recon_error_k*.csv).
KSWEEP_DIR <- file.path(INPUT_DIR, "ksweep")

Ks       <- c(15, 20, 25, 30, 35, 40)
CHOSEN_K <- 30   # the manuscript's K, marked in the figures

# Factors summed for the "top-N" composite. The compositional framework
# represents the malignant programme as a sum of several factors, so the
# single best-matching factor is the conservative summary and the summed
# composite the less conservative one. Applied to nerve and tumour alike.
N_TOP <- 3

nmf_fit_path <- function(K) {
  file.path(KSWEEP_DIR, sprintf("PNI_NMF_L105-k%02d.rds", K))
}

# Standardise factor names to NMF_<n>; create them if the fit has none
factor_names <- function(nm, n) {
  if (is.null(nm)) return(sprintf("NMF_%d", seq_len(n)))
  sub("^nmf(\\d+)$", "NMF_\\1", nm, perl = TRUE)
}

cosine_matrix <- function(w_mat) {
  norms <- sqrt(colSums(w_mat^2))
  if (any(!is.finite(norms)) || any(norms == 0)) {
    stop("Cannot calculate cosine similarity: a loading vector has zero or non-finite norm (",
         paste(colnames(w_mat)[!is.finite(norms) | norms == 0], collapse = ", "), ").")
  }
  crossprod(w_mat) / outer(norms, norms)
}

# ------------------------------------------------------------
# 1. Load RCTD-derived cell-type fractions
# ------------------------------------------------------------
rctd <- readRDS(file.path(INPUT_DIR, "RCTD_weights.Rds"))
stopifnot("barcode" %in% colnames(rctd))
stopifnot(all(c("Schwann", "Ductal_malignant") %in% colnames(rctd)))

rctd_sub <- rctd %>%
  dplyr::select(barcode, schwann_frac = Schwann, malignant_frac = Ductal_malignant)
if (anyDuplicated(rctd_sub$barcode)) stop("Duplicate barcodes in the RCTD table.")

cat("RCTD table:", nrow(rctd_sub), "spots\n")

# ------------------------------------------------------------
# 2. Combine the per-K reconstruction-error csvs into one table,
#    checking exactly one row per K that matches its filename
# ------------------------------------------------------------
recon_list <- list()
for (K in Ks) {
  f <- file.path(KSWEEP_DIR, sprintf("recon_error_k%02d.csv", K))
  # A partial sweep must not be interpreted as a complete sensitivity analysis
  if (!file.exists(f)) {
    stop("Missing reconstruction-error file: ", f)
  }
  df <- read.csv(f)
  if (!all(c("K", "rel_frob_error") %in% colnames(df))) {
    stop(f, " lacks the required columns K and rel_frob_error.")
  }
  if (nrow(df) != 1) {
    stop(f, " has ", nrow(df), " rows; expected exactly 1 (the saved fit's error).")
  }
  if (is.na(df$K) || df$K != K) {
    stop(f, " reports an incorrect or missing K value (", df$K, "); its filename says K = ", K, ".")
  }
  recon_list[[length(recon_list) + 1]] <- df
}
recon_error <- bind_rows(recon_list)
if (!setequal(recon_error$K, Ks) || anyDuplicated(recon_error$K)) {
  stop("Reconstruction-error table does not contain exactly the expected K values (",
       paste(Ks, collapse = ", "), ").")
}

write.csv(recon_error, file.path(TABLES_DIR, "NMF_K_reconstruction_error.csv"), row.names = FALSE)
cat("Combined reconstruction-error table (", nrow(recon_error), "rows) saved to:\n ",
    file.path(TABLES_DIR, "NMF_K_reconstruction_error.csv"), "\n\n")
print(recon_error)

# ------------------------------------------------------------
# 3. For each K: load the fit, compute scale-aware H, correlate every
#    factor with the RCTD-derived Schwann / malignant fractions.
#    Gene loadings (W) are cached here, after checking that every fit
#    uses the same genes in the same order.
# ------------------------------------------------------------

# Top-N summed composite: sum the exposures of the (up to) N_TOP factors
# with the highest POSITIVE correlation, then correlate the sum with the
# RCTD-derived fraction. Negatively correlated factors are never summed
# (they would dilute, not strengthen, the composite); undefined (NA)
# correlations are ignored.
top_n_summed_cor <- function(cor_vec, Hc, target) {
  valid_idx <- which(is.finite(cor_vec) & cor_vec > 0)
  if (length(valid_idx) == 0) {
    return(list(cor = NA_real_, factors = character(0)))
  }
  ord     <- valid_idx[order(cor_vec[valid_idx], decreasing = TRUE)]
  top_idx <- head(ord, N_TOP)
  summed  <- colSums(Hc[top_idx, , drop = FALSE])
  list(cor = suppressWarnings(cor(summed, target, use = "complete.obs")),
       factors = rownames(Hc)[top_idx])
}

# Correlation of every factor with a target; zero-variance factors give
# NA instead of a warning and are excluded from selection
factor_cor <- function(Hc, target) {
  r <- apply(Hc, 1, function(x) suppressWarnings(cor(x, target, use = "complete.obs")))
  r[!is.finite(r)] <- NA_real_
  r
}

w_by_K   <- list()   # gene loadings per K, names "K15", "K20", ...
gene_ids <- NULL     # reference gene order (from the first fit loaded)
k_results_list <- list()

for (K in Ks) {
  fit_path <- nmf_fit_path(K)
  if (!file.exists(fit_path)) {
    stop(sprintf("Missing NMF fit for K=%d: %s", K, fit_path))
  }
  fit <- readRDS(fit_path)
  
  # ---- Gene alignment check (W rows must be the same genes, same order)
  w <- fit@w
  colnames(w) <- factor_names(colnames(w), ncol(w))
  if (is.null(rownames(w))) {
    stop(sprintf("K=%d: gene loadings have no gene names, so gene alignment across K cannot be verified.", K))
  }
  if (is.null(gene_ids)) {
    gene_ids <- rownames(w)
  } else if (!identical(rownames(w), gene_ids)) {
    stop(sprintf("K=%d: gene identities or order differ from the first fit; align genes before comparing loadings.", K))
  }
  w_by_K[[sprintf("K%02d", K)]] <- w
  
  # ---- Scale-aware exposures, same convention as the main workflow:
  # H <- sweep(h, 1, d, `*`)  (factors x spots)
  H <- sweep(fit@h, 1, fit@d, `*`)
  rownames(H) <- factor_names(rownames(H), nrow(H))
  stopifnot(identical(rownames(H), colnames(w)))
  rm(fit)
  
  common <- intersect(colnames(H), rctd_sub$barcode)
  if (length(common) == 0) {
    warning(sprintf("K=%d: no overlapping barcodes between NMF fit and RCTD table -- check barcode formatting", K))
    next
  }
  cat(sprintf("K=%d: %d / %d NMF spots matched to RCTD\n", K, length(common), ncol(H)))
  
  Hc <- H[, common, drop = FALSE]
  rc <- rctd_sub[match(common, rctd_sub$barcode), ]
  
  schwann_cor <- factor_cor(Hc, rc$schwann_frac)
  malig_cor   <- factor_cor(Hc, rc$malignant_frac)
  
  if (all(is.na(schwann_cor)) || all(is.na(malig_cor))) {
    warning(sprintf("K=%d: no valid factor correlations for one or both cell types -- skipping", K))
    next
  }
  
  # Signed maximum, not maximum absolute value: a strongly negative
  # correlation does not identify a nerve (or tumour) factor.
  nerve_idx  <- which.max(replace(schwann_cor, is.na(schwann_cor), -Inf))
  tumour_idx <- which.max(replace(malig_cor,   is.na(malig_cor),   -Inf))
  
  nerve_top  <- top_n_summed_cor(schwann_cor, Hc, rc$schwann_frac)
  tumour_top <- top_n_summed_cor(malig_cor,   Hc, rc$malignant_frac)
  
  k_results_list[[length(k_results_list) + 1]] <- data.frame(
    K                  = K,
    n_spots_matched    = length(common),
    n_factors_na_cor   = sum(is.na(schwann_cor) | is.na(malig_cor)),
    nerve_factor       = names(schwann_cor)[nerve_idx],
    nerve_max_cor      = unname(schwann_cor[nerve_idx]),
    nerve_top_factors  = paste(nerve_top$factors, collapse = ";"),
    nerve_top_n_used   = length(nerve_top$factors),
    nerve_top_cor      = nerve_top$cor,
    tumour_factor      = names(malig_cor)[tumour_idx],
    tumour_max_cor     = unname(malig_cor[tumour_idx]),
    tumour_top_factors = paste(tumour_top$factors, collapse = ";"),
    tumour_top_n_used  = length(tumour_top$factors),
    tumour_top_cor     = tumour_top$cor
  )
}

if (length(k_results_list) == 0) stop("No K could be evaluated -- check the input files.")
k_results <- bind_rows(k_results_list)

# Every fit exists (checked above), but a K can still be skipped for a
# data problem (no matching barcodes, no valid correlations). Say so loudly.
skipped_K <- setdiff(Ks, k_results$K)
if (length(skipped_K) > 0) {
  warning("Correlation results are missing for K = ", paste(skipped_K, collapse = ", "),
          " (see warnings above). The sweep is incomplete for the correlation analysis.",
          call. = FALSE)
}
print(k_results)

# ------------------------------------------------------------
# 4. Combine with reconstruction error, save summary table
# ------------------------------------------------------------
summary_tbl <- left_join(k_results, recon_error, by = "K")
stopifnot(nrow(summary_tbl) == nrow(k_results))   # one recon row per K (checked in Section 2)
write.csv(summary_tbl, file.path(TABLES_DIR, "K_sensitivity_summary.csv"), row.names = FALSE)
cat("\n--- Final summary table ---\n")
print(summary_tbl)

# ------------------------------------------------------------
# 4b. Within-K redundancy of the summed top-N factors
#
# Cosine similarity of gene loadings (W) rather than spot exposures (H),
# as in the main workflow: a value near 1 means near-duplicate
# programmes, well below 1 means distinct programmes contributing to the
# same RCTD-derived fraction.
# ------------------------------------------------------------
redundancy_rows <- list()
for (i in seq_len(nrow(summary_tbl))) {
  K <- summary_tbl$K[i]
  w <- w_by_K[[sprintf("K%02d", K)]]
  for (ct in c("nerve", "tumour")) {
    facs <- strsplit(summary_tbl[[paste0(ct, "_top_factors")]][i], ";", fixed = TRUE)[[1]]
    if (length(facs) < 2) next
    cm <- cosine_matrix(w[, facs, drop = FALSE])
    pairs <- which(upper.tri(cm), arr.ind = TRUE)
    redundancy_rows[[length(redundancy_rows) + 1]] <- data.frame(
      K         = K,
      cell_type = ct,
      factor_a  = facs[pairs[, 1]],
      factor_b  = facs[pairs[, 2]],
      cosine_w  = round(cm[pairs], 3)
    )
  }
}
top_factor_redundancy <- bind_rows(redundancy_rows)
write.csv(top_factor_redundancy, file.path(TABLES_DIR, "top_factor_redundancy.csv"), row.names = FALSE)
cat("\n--- Within-K cosine similarity of the summed top factors (gene loadings) ---\n")
print(top_factor_redundancy)

# ------------------------------------------------------------
# 4c. Cross-K consistency of the selected factor's gene loadings
#
# Each K is a separate NMF run, and its factor is selected on spot-level
# exposures (H) while the comparison is made on gene loadings (W), which
# the selection does not directly constrain. Gene order was verified to
# be identical across fits in Section 3, so loadings are comparable.
# ------------------------------------------------------------
collect_best_w <- function(factor_col_name) {
  vecs <- list()
  for (i in seq_len(nrow(summary_tbl))) {
    K   <- summary_tbl$K[i]
    fac <- summary_tbl[[factor_col_name]][i]
    if (is.na(fac)) next
    vecs[[sprintf("K%02d_%s", K, fac)]] <- w_by_K[[sprintf("K%02d", K)]][, fac]
  }
  do.call(cbind, vecs)
}

cross_k_cosine <- function(w_mat, label) {
  m <- cosine_matrix(w_mat)
  diag(m) <- 1
  cat(sprintf("\n--- %s: cross-K cosine similarity (selected factor's gene loadings) ---\n", label))
  print(round(m, 3))
  write.csv(round(m, 4), file.path(TABLES_DIR, sprintf("crossK_cosine_%s.csv", label)))
  
  pdf(file.path(HEATMAPS_DIR, sprintf("crossK_cosine_%s.pdf", label)), width = 6, height = 5)
  pheatmap(m,
           cluster_rows = FALSE, cluster_cols = FALSE,   # keep K order
           display_numbers = TRUE, number_format = "%.2f",
           main = sprintf("%s: gene-loading similarity across K fits\n(rows/cols in K order)", label),
           fontsize = 9)
  dev.off()
  m
}

nerve_crossK_mat  <- cross_k_cosine(collect_best_w("nerve_factor"),  "Nerve")
tumour_crossK_mat <- cross_k_cosine(collect_best_w("tumour_factor"), "Tumour")

cat("\nCross-K heatmaps saved to", HEATMAPS_DIR, "\n")
cat("============================================================================\n")

# ------------------------------------------------------------
# 5. Two-panel supplementary figure
# ------------------------------------------------------------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.line  = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    plot.title = element_text(face = "bold", size = 12)
  )

p_recon <- ggplot(summary_tbl, aes(x = K, y = rel_frob_error)) +
  geom_line(color = "grey30") +
  geom_point(size = 2.5, color = "black") +
  geom_vline(xintercept = CHOSEN_K, linetype = "dashed", color = "firebrick", alpha = 0.6) +
  scale_x_continuous(breaks = Ks) +
  labs(title = "Reconstruction error vs K",
       x = "K (number of factors)",
       y = "Relative Frobenius\nreconstruction error (in-sample)") +
  theme_pub

# Long format: colour = cell type, line type / shape = summary method,
# so all four lines are distinguishable in the legend
TOP_LABEL    <- sprintf("Top %d summed", N_TOP)
SINGLE_LABEL <- "Single best factor"
cor_long <- bind_rows(
  data.frame(K = summary_tbl$K, cell_type = "Nerve (Schwann)",  summary = SINGLE_LABEL, cor = summary_tbl$nerve_max_cor),
  data.frame(K = summary_tbl$K, cell_type = "Nerve (Schwann)",  summary = TOP_LABEL,    cor = summary_tbl$nerve_top_cor),
  data.frame(K = summary_tbl$K, cell_type = "Tumour (malignant)", summary = SINGLE_LABEL, cor = summary_tbl$tumour_max_cor),
  data.frame(K = summary_tbl$K, cell_type = "Tumour (malignant)", summary = TOP_LABEL,    cor = summary_tbl$tumour_top_cor)
) %>%
  mutate(summary = factor(summary, levels = c(SINGLE_LABEL, TOP_LABEL)))

y_min <- min(0, cor_long$cor, na.rm = TRUE)

p_corr <- ggplot(cor_long, aes(x = K, y = cor, color = cell_type, linetype = summary, shape = summary,
                               group = interaction(cell_type, summary))) +
  geom_line() +
  geom_point(size = 2.3) +
  geom_vline(xintercept = CHOSEN_K, linetype = "dashed", color = "firebrick", alpha = 0.6) +
  scale_x_continuous(breaks = Ks) +
  coord_cartesian(ylim = c(y_min, 1)) +
  scale_color_manual(values = c("Nerve (Schwann)" = "#E08214", "Tumour (malignant)" = "gray35"),
                     name = NULL) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), c(SINGLE_LABEL, TOP_LABEL)), name = NULL) +
  scale_shape_manual(values = setNames(c(16, 1), c(SINGLE_LABEL, TOP_LABEL)), name = NULL) +
  labs(title = "Correspondence with RCTD-derived fractions vs K",
       x = "K (number of factors)",
       y = "Correlation with RCTD-derived\ncell-type fraction") +
  theme_pub +
  theme(legend.position = "bottom", legend.box = "vertical")

p_panel <- p_recon + p_corr + plot_annotation(
  title = "NMF K-sensitivity: reconstruction quality and correspondence with RCTD-derived cell-type fractions",
  theme = theme(plot.title = element_text(size = 12, face = "bold"))
)
print(p_panel)

ggsave(file.path(FIGURES_DIR, "K_sensitivity_panel.pdf"), p_panel, width = 13, height = 5.5)

cat("\nDone. Outputs:\n")
cat(" ", file.path(TABLES_DIR, "K_sensitivity_summary.csv"), "\n")
cat(" ", file.path(TABLES_DIR, "top_factor_redundancy.csv"), "\n")
cat(" ", file.path(FIGURES_DIR, "K_sensitivity_panel.pdf"), "\n")
cat("Red dashed line marks K =", CHOSEN_K, "(the manuscript's chosen K) in both panels.\n")
