# ============================================================
# 01b_GSEA_plots.R
#
# READ-AND-PLOT ONLY: loads the .rds GSEA result objects written
# by 01a_GSEA_compute.R and regenerates:
#   - dotplots split by NES positive / NES negative
#   - readable (gene-symbol) result tables
# One flat output folder per comparison.
#
# Input : processed/bulkRNAseq/GSEA/<run_label>/<comparison>/**/*.rds
# Output: output/bulkRNAseq/tables/GSEA/<run_label>/<comparison>/...
#         output/bulkRNAseq/figures/GSEA/<run_label>/<comparison>/...
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
step_root <- file.path("processed", "bulkRNAseq", "GSEA")

# Pick up the run written by 01a_GSEA_compute.R automatically.
# To regenerate plots for an older run instead, just hardcode run_label.
run_label <- readLines(file.path(step_root, "_LATEST_RUN_LABEL.txt"), n = 1)

in_dir       <- file.path(step_root, run_label)
tables_dir   <- file.path("output", "bulkRNAseq", "tables", "GSEA", run_label, "REGENERATED_BY_COMPARISON")
figures_dir  <- file.path("output", "bulkRNAseq", "figures", "GSEA", run_label, "REGENERATED_BY_COMPARISON")
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

comparisons <- c("DRAXIN_vs_ctrl", "diluted_vs_ctrl", "DRAXIN_vs_diluted")

regen_cfg <- list(
  showCategory = 20,
  width = 9,
  height = 6,
  dpi = 300,
  base_size = 8,
  title_size = 10,
  axis_text_size = 7,
  axis_title_size = 8,
  strip_text_size = 8,
  legend_text_size = 7,
  legend_title_size = 8
)

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------
safe_ggsave <- function(path, plot, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi),
           error = function(e) message("WARN: ggsave failed: ", path, " | ", e$message))
  invisible(TRUE)
}

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(readr::write_csv(x, path), error = function(e) message("WARN: write_csv failed: ", path, " | ", e$message))
  invisible(TRUE)
}

is_enrichment_rds <- function(x) inherits(x, "gseaResult") || inherits(x, "enrichResult")

safe_result_tbl <- function(res) {
  tryCatch(tibble::as_tibble(res@result), error = function(e) {
    message("WARN: could not extract @result | ", e$message)
    tibble::tibble()
  })
}

make_readable_if_possible <- function(res) {
  tryCatch(
    clusterProfiler::setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) {
      message("WARN: setReadable failed | ", e$message)
      res
    }
  )
}

subset_result_object <- function(res, keep_idx) {
  res2 <- res
  res_tbl <- safe_result_tbl(res2)

  if (nrow(res_tbl) == 0 || length(keep_idx) != nrow(res_tbl)) return(NULL)

  res2@result <- res_tbl[keep_idx, , drop = FALSE]
  if (nrow(res2@result) == 0) return(NULL)
  res2
}

make_flat_name <- function(file_path, comparison_name, out_dir_root) {
  rel <- sub(
    paste0("^", gsub("([\\^\\$\\.\\|\\(\\)\\[\\]\\*\\+\\?\\\\])", "\\\\\\1", out_dir_root), "/?"),
    "",
    file_path
  )
  rel <- sub("\\.rds$", "", rel)
  rel <- sub(paste0("^", comparison_name, "/?"), "", rel)
  rel <- gsub("/", "__", rel)
  rel <- gsub("\\s+", "_", rel)
  rel
}

build_dotplot <- function(res, title_text, cfg = regen_cfg) {
  res_tbl <- safe_result_tbl(res)
  if (nrow(res_tbl) == 0) return(NULL)

  ncat <- min(cfg$showCategory, nrow(res_tbl))

  tryCatch({
    enrichplot::dotplot(res, showCategory = ncat) +
      ggplot2::ggtitle(title_text) +
      ggplot2::theme_bw(base_size = cfg$base_size) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(size = cfg$title_size),
        axis.text    = ggplot2::element_text(size = cfg$axis_text_size),
        axis.title   = ggplot2::element_text(size = cfg$axis_title_size),
        strip.text   = ggplot2::element_text(size = cfg$strip_text_size),
        legend.text  = ggplot2::element_text(size = cfg$legend_text_size),
        legend.title = ggplot2::element_text(size = cfg$legend_title_size)
      )
  }, error = function(e) {
    message("WARN: dotplot failed for ", title_text, " | ", e$message)
    NULL
  })
}

save_plot_pair <- function(plot_obj, out_base, cfg = regen_cfg) {
  if (is.null(plot_obj)) return(invisible(FALSE))
  safe_ggsave(paste0(out_base, ".png"), plot_obj, width = cfg$width, height = cfg$height, dpi = cfg$dpi)
  safe_ggsave(paste0(out_base, ".pdf"), plot_obj, width = cfg$width, height = cfg$height, dpi = cfg$dpi)
  invisible(TRUE)
}

regen_one_rds <- function(f, comparison_name, comparison_tables_dir, comparison_figures_dir, out_dir_root, cfg = regen_cfg) {
  res <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(res) || !is_enrichment_rds(res)) return(invisible(FALSE))

  # readable CSV export -> tables
  res_readable <- make_readable_if_possible(res)
  readable_tbl <- safe_result_tbl(res_readable)

  flat_name <- make_flat_name(f, comparison_name = comparison_name, out_dir_root = out_dir_root)

  if (nrow(readable_tbl) > 0) {
    safe_write_csv(readable_tbl, file.path(comparison_tables_dir, paste0(flat_name, "__readable.csv")))
  }

  # plots -> figures
  res_tbl <- safe_result_tbl(res)
  if (nrow(res_tbl) == 0) return(invisible(TRUE))

  if ("NES" %in% colnames(res_tbl)) {
    pos_idx <- !is.na(res_tbl$NES) & res_tbl$NES > 0
    neg_idx <- !is.na(res_tbl$NES) & res_tbl$NES < 0

    res_pos <- subset_result_object(res, pos_idx)
    res_neg <- subset_result_object(res, neg_idx)

    if (!is.null(res_pos) && nrow(res_pos@result) > 0) {
      p_pos <- build_dotplot(res_pos, title_text = paste0(flat_name, " | NES positive"), cfg = cfg)
      save_plot_pair(p_pos, file.path(comparison_figures_dir, paste0(flat_name, "__NES_positive_dotplot_REGEN")), cfg = cfg)
    }
    if (!is.null(res_neg) && nrow(res_neg@result) > 0) {
      p_neg <- build_dotplot(res_neg, title_text = paste0(flat_name, " | NES negative"), cfg = cfg)
      save_plot_pair(p_neg, file.path(comparison_figures_dir, paste0(flat_name, "__NES_negative_dotplot_REGEN")), cfg = cfg)
    }
  } else {
    p_all <- build_dotplot(res, title_text = paste0(flat_name, " | all terms"), cfg = cfg)
    save_plot_pair(p_all, file.path(comparison_figures_dir, paste0(flat_name, "__dotplot_REGEN")), cfg = cfg)
  }

  invisible(TRUE)
}

# ------------------------------------------------------------
# main loop: one flat output folder per comparison
# ------------------------------------------------------------
regen_status <- list()

for (nm in comparisons) {
  message("=== Regenerating outputs for comparison: ", nm, " ===")

  comp_dir <- file.path(in_dir, nm)
  if (!dir.exists(comp_dir)) {
    warning("Comparison output folder not found: ", comp_dir)
    next
  }

  comparison_tables_dir  <- file.path(tables_dir, nm)
  comparison_figures_dir <- file.path(figures_dir, nm)
  dir.create(comparison_tables_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(comparison_figures_dir, recursive = TRUE, showWarnings = FALSE)

  rds_files <- list.files(comp_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  rds_files <- rds_files[!grepl("ranks_ENTREZ\\.rds$", rds_files)]

  # if both X.rds and X__readable.rds exist, keep only X.rds
  plain_rds    <- rds_files[!grepl("__readable\\.rds$", rds_files)]
  readable_rds <- rds_files[grepl("__readable\\.rds$", rds_files)]
  plain_key    <- sub("\\.rds$", "", plain_rds)
  readable_key <- sub("__readable\\.rds$", "", readable_rds)
  chosen_rds   <- c(plain_rds, readable_rds[!readable_key %in% plain_key])

  n_done <- 0L
  n_skip <- 0L

  for (f in chosen_rds) {
    ok <- tryCatch(
      regen_one_rds(
        f = f,
        comparison_name = nm,
        comparison_tables_dir = comparison_tables_dir,
        comparison_figures_dir = comparison_figures_dir,
        out_dir_root = in_dir
      ),
      error = function(e) {
        message("WARN: failed on ", f, " | ", e$message)
        FALSE
      }
    )
    if (isTRUE(ok)) n_done <- n_done + 1L else n_skip <- n_skip + 1L
  }

  regen_status[[nm]] <- tibble::tibble(
    comparison = nm,
    n_rds_processed = n_done,
    n_rds_skipped = n_skip,
    tables_folder = comparison_tables_dir,
    figures_folder = comparison_figures_dir
  )
}

regen_status_df <- dplyr::bind_rows(regen_status)
safe_write_csv(regen_status_df, file.path(tables_dir, "REGEN_STATUS.csv"))

message("Done.")
message("Regenerated tables written under:  ", tables_dir)
message("Regenerated figures written under: ", figures_dir)
