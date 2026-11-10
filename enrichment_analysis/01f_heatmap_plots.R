# ============================================================
# 01f_heatmap_plots.R
#
# READ-AND-PLOT ONLY: loads the merged leading-edge bundle written by
# 01e_heatmap_compute.R and draws one heatmap per target term (batch-
# corrected expression, z-scored per row, columns split by condition).
#
# Input : processed/bulkRNAseq/heatmap/merged_leadingedge_bundle.rds
# Output: output/bulkRNAseq/figures/heatmap/*.pdf
# ============================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(matrixStats)
})

# ------------------------------------------------------------
# paths + bundle
# ------------------------------------------------------------
processed_dir <- file.path("processed", "bulkRNAseq", "heatmap")
figures_dir   <- file.path("output", "bulkRNAseq", "figures", "heatmap")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

bundle_path <- file.path(processed_dir, "merged_leadingedge_bundle.rds")
stopifnot(file.exists(bundle_path))
bundle <- readRDS(bundle_path)

expr_mat       <- bundle$expr_mat_bc
meta           <- bundle$meta
batch_var      <- bundle$batch_var
exp_var        <- bundle$exp_var
comparisons    <- bundle$comparisons
targets        <- bundle$targets
merged_results <- bundle$merged_results

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------
zscore_rows <- function(m) t(scale(t(as.matrix(m))))

# ------------------------------------------------------------
# color palettes
# ------------------------------------------------------------
cond_cols <- c("0nM" = "#1b9e77", "17nM" = "#d95f02", "68nM" = "#7570b3")

make_level_colors <- function(x) {
  lv <- levels(as.factor(x))
  cols <- circlize::rand_color(length(lv))
  names(cols) <- lv
  cols
}
batch_cols <- make_level_colors(meta[[batch_var]])
exp_cols   <- make_level_colors(meta[[exp_var]])

# ------------------------------------------------------------
# heatmap plotting (MERGED) + optional row membership annotation
# ------------------------------------------------------------
plot_leadingedge_heatmap_merged <- function(expr_mat, meta, ids, label_map, mem_df,
                                            title, out_pdf,
                                            top_n = 80,
                                            show_membership = TRUE) {
  if (length(ids) < 3) {
    message("  Skipping (too few mapped genes): ", title, " (n=", length(ids), ")")
    return(invisible(NULL))
  }

  # sample order: condition -> batch -> exp -> sample
  keep_samps <- rownames(meta)[meta$condition_nm %in% c("0nM", "17nM", "68nM")]
  ord_df <- meta[keep_samps, , drop = FALSE]
  ord_df$.sample <- rownames(ord_df)

  o <- order(ord_df$condition_nm, ord_df[[batch_var]], ord_df[[exp_var]], ord_df$.sample)
  samp_order <- ord_df$.sample[o]

  m <- expr_mat[ids, samp_order, drop = FALSE]

  # row labels
  rlab <- unname(label_map[rownames(m)])
  rlab[is.na(rlab) | rlab == ""] <- rownames(m)[is.na(rlab) | rlab == ""]
  rownames(m) <- make.unique(rlab)

  # if too many genes: keep top variable genes, but keep membership aligned
  if (is.finite(top_n) && nrow(m) > top_n) {
    vv <- matrixStats::rowVars(as.matrix(m), na.rm = TRUE)
    ord <- order(vv, decreasing = TRUE)
    keep <- rownames(m)[ord][seq_len(top_n)]
    m <- m[keep, , drop = FALSE]

    label_for_id <- unname(label_map[mem_df$expr_id])
    label_for_id[is.na(label_for_id) | label_for_id == ""] <- mem_df$expr_id[is.na(label_for_id) | label_for_id == ""]
    label_for_id <- make.unique(label_for_id)
    mem_df$._label_unique <- label_for_id
    mem_df <- mem_df[mem_df$._label_unique %in% rownames(m), , drop = FALSE]
  }

  mz <- zscore_rows(m)
  meta_aligned <- meta[colnames(mz), , drop = FALSE]

  ha_top <- HeatmapAnnotation(
    condition = meta_aligned$condition_nm,
    batch     = meta_aligned[[batch_var]],
    exp       = meta_aligned[[exp_var]],
    col = list(condition = cond_cols, batch = batch_cols, exp = exp_cols),
    annotation_legend_param = list(
      condition = list(title = "Condition"),
      batch     = list(title = batch_var),
      exp       = list(title = exp_var)
    )
  )

  ra <- NULL
  if (isTRUE(show_membership) && nrow(mem_df)) {
    mem_use <- mem_df
    if (!"._label_unique" %in% names(mem_use)) {
      label_for_id <- unname(label_map[mem_use$expr_id])
      label_for_id[is.na(label_for_id) | label_for_id == ""] <- mem_use$expr_id[is.na(label_for_id) | label_for_id == ""]
      mem_use$._label_unique <- make.unique(label_for_id)
    }
    mem_use <- mem_use[match(rownames(mz), mem_use$._label_unique), , drop = FALSE]

    for (cc in comparisons) mem_use[[cc]][is.na(mem_use[[cc]])] <- FALSE

    ra <- rowAnnotation(
      `LE in...` = ComplexHeatmap::anno_simple(
        as.matrix(mem_use[, comparisons, drop = FALSE]) * 1,
        col = c("0" = "white", "1" = "black"),
        gp = grid::gpar(col = NA)
      ),
      annotation_name_side = "top",
      annotation_name_gp = grid::gpar(fontsize = 8)
    )
  }

  ht <- Heatmap(
    mz, name = "z",
    top_annotation = ha_top,
    left_annotation = ra,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows    = FALSE,
    show_row_names  = TRUE,
    row_names_gp    = grid::gpar(fontsize = 7),
    column_split    = meta_aligned$condition_nm,
    column_title    = title
  )

  grDevices::pdf(out_pdf, width = 12, height = 6.5, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()
  message("  Wrote: ", out_pdf)
}

# ------------------------------------------------------------
# run: one heatmap per term (already merged across comparisons by 01e)
# ------------------------------------------------------------
for (tgt in targets) {
  message("\n===============================")
  message("Merged term: ", tgt$label, " | ", tgt$collection)
  message("===============================")

  merged <- merged_results[[tgt$label]]

  if (is.null(merged) || length(merged$ids) < 3) {
    message("  Skipping: insufficient genes after merge for ", tgt$label)
    next
  }

  title <- paste0(tgt$label, " | LE union across: ", paste(comparisons, collapse = " / "))
  out_pdf <- file.path(figures_dir, paste0(tgt$label, "__MERGED_LEADINGEDGE_ALLCONDS_BATCHCORR_heatmap.pdf"))

  plot_leadingedge_heatmap_merged(
    expr_mat  = expr_mat,
    meta      = meta,
    ids       = merged$ids,
    label_map = merged$label_map,
    mem_df    = merged$mem,
    title     = title,
    out_pdf   = out_pdf,
    top_n     = 80,          # set Inf to show ALL union genes
    show_membership = FALSE  # set TRUE to add row block showing which comps contributed
  )
}

message("\nDone. Merged leading-edge heatmaps in: ", figures_dir)

