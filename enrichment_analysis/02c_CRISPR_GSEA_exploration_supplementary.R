# ============================================================
# CRISPR_GSEA_exploration_supplementary.R
#
# SUPPLEMENTARY / OPTIONAL -- not part of the paper's reported figure
# (see 02b_CRISPR_GSEA_plots.R for that). This is broader exploration
# across the full GSEA result set:
#   - both filter modes: p.adjust <= 0.30 AND pvalue < 0.05
#   - all collections (Hallmark, C2, C5, Reactome), not just Hallmark+C5
#   - per-result dotplots (one per collection x arm)
#   - paired resistance-vs-sensitisation NES scatter (pathways enriched
#     with NES > 0 in BOTH arms)
#
# Input : processed/CRISPR_GSEA/<run_label>/GSEA_summary_compact.csv
#         processed/CRISPR_GSEA/<run_label>/{MSigDB_*,Reactome}/*.rds
# Output: output/CRISPR_GSEA/tables/supplementary/...
#         output/CRISPR_GSEA/figures/supplementary/...
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(enrichplot)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
step_root <- file.path("processed", "CRISPR_GSEA")
run_label <- readLines(file.path(step_root, "_LATEST_RUN_LABEL.txt"), n = 1)
run_dir   <- file.path(step_root, run_label)

summary_file <- file.path(run_dir, "GSEA_summary_compact.csv")
if (!file.exists(summary_file)) {
  stop("Could not find summary file: ", summary_file, " -- run 02a_CRISPR_GSEA_compute.R first.")
}
summary_tbl_small <- readr::read_csv(summary_file, show_col_types = FALSE)

tables_dir  <- file.path("output", "CRISPR_GSEA", "tables",  "supplementary")
figures_dir <- file.path("output", "CRISPR_GSEA", "figures", "supplementary")
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(readr::write_csv(x, path), error = function(e) message("WARN: write_csv failed: ", path, " | ", e$message))
  invisible(TRUE)
}

safe_ggsave <- function(path, plot, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi),
           error = function(e) message("WARN: ggsave failed: ", path, " | ", e$message))
  invisible(TRUE)
}

# ------------------------------------------------------------
# 1) per-result dotplots (one per collection x arm), all collections
# ------------------------------------------------------------
result_dirs <- c(
  file.path(run_dir, "MSigDB_Hallmark"),
  file.path(run_dir, "MSigDB_C2"),
  file.path(run_dir, "MSigDB_C5"),
  file.path(run_dir, "Reactome")
)

for (d in result_dirs) {
  if (!dir.exists(d)) next
  rds_files <- list.files(d, pattern = "\\.rds$", full.names = TRUE)

  for (f in rds_files) {
    res <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(res) || !inherits(res, c("gseaResult", "enrichResult"))) next

    res_tbl <- tryCatch(tibble::as_tibble(res@result), error = function(e) tibble::tibble())
    if (nrow(res_tbl) == 0) next

    prefix <- sub("\\.rds$", "", basename(f))

    p <- tryCatch(
      enrichplot::dotplot(res, showCategory = min(20, nrow(res_tbl))) + ggplot2::ggtitle(prefix),
      error = function(e) NULL
    )
    if (!is.null(p)) {
      safe_ggsave(file.path(figures_dir, paste0(prefix, "_dotplot.png")), p)
      safe_ggsave(file.path(figures_dir, paste0(prefix, "_dotplot.pdf")), p)
    }
  }
}

# ------------------------------------------------------------
# helper for leading edge count
# ------------------------------------------------------------
count_core_genes <- function(core_enrichment_string) {
  if (is.na(core_enrichment_string) || core_enrichment_string == "") return(NA_integer_)
  length(strsplit(core_enrichment_string, "/", fixed = TRUE)[[1]])
}

make_filtered_view <- function(summary_tbl_small,
                               filter_mode = c("padj_0.3", "pval_0.05"),
                               min_leading_edge_n = 3,
                               min_set_size_for_plot = 10) {
  filter_mode <- match.arg(filter_mode)

  x <- summary_tbl_small %>%
    dplyr::select(any_of(c(
      "ID", "Description", "setSize", "NES", "pvalue", "p.adjust",
      "collection", "arm", "core_enrichment"
    ))) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      arm = dplyr::recode(
        arm,
        "positive_selection" = "Resistance arm",
        "negative_selection" = "Sensitisation arm"
      ),
      leading_edge_n = vapply(core_enrichment, count_core_genes, integer(1))
    ) %>%
    dplyr::filter(
      !is.na(arm), !is.na(NES), NES > 0,
      !is.na(setSize), setSize >= min_set_size_for_plot,
      !is.na(leading_edge_n), leading_edge_n >= min_leading_edge_n
    )

  if (filter_mode == "padj_0.3") {
    x <- x %>% dplyr::filter(!is.na(p.adjust), p.adjust <= 0.30)
    label <- paste0("padj_le_0.3_LE", min_leading_edge_n)
    sort_col <- "p.adjust"
  } else {
    x <- x %>% dplyr::filter(!is.na(pvalue), pvalue < 0.05)
    label <- paste0("pvalue_lt_0.05_LE", min_leading_edge_n)
    sort_col <- "pvalue"
  }

  list(data = x, label = label, sort_col = sort_col)
}

# ------------------------------------------------------------
# 2) full exploration for one filter mode: NES compare + paired scatter
# ------------------------------------------------------------
run_exploration_for_mode <- function(summary_tbl_small, filter_mode, top_n = 15,
                                     min_leading_edge_n = 3, min_set_size_for_plot = 10) {

  view <- make_filtered_view(summary_tbl_small, filter_mode, min_leading_edge_n, min_set_size_for_plot)
  x <- view$data
  label <- view$label
  sort_col <- view$sort_col

  safe_write_csv(x, file.path(tables_dir, paste0("GSEA_summary_compact_filtered_", label, ".csv")))

  if (nrow(x) == 0) {
    message("No enriched terms passed final filter: ", label)
    return(invisible(NULL))
  }

  # -----------------------------------------
  # top terms per arm/collection (all collections)
  # -----------------------------------------
  nes_compare <- x %>%
    dplyr::group_by(collection, arm) %>%
    dplyr::arrange(.data[[sort_col]], dplyr::desc(NES), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  safe_write_csv(nes_compare, file.path(tables_dir, paste0("NES_compare_top_terms_", label, ".csv")))

  p_nes <- ggplot(
    nes_compare,
    aes(x = reorder(Description, NES), y = NES, size = setSize)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey60") +
    geom_point() +
    coord_flip() +
    facet_grid(collection ~ arm, scales = "free_y", space = "free_y") +
    labs(
      title = paste0(
        "Top enriched pathways by selection arm (", label, ")\n",
        "Filtered for leading-edge genes >= ", min_leading_edge_n
      ),
      x = NULL, y = "Normalized Enrichment Score (NES)", size = "Gene set size"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  safe_ggsave(file.path(figures_dir, paste0("NES_compare_top_terms_", label, ".png")), p_nes, width = 11, height = 12, dpi = 300)
  safe_ggsave(file.path(figures_dir, paste0("NES_compare_top_terms_", label, ".pdf")), p_nes, width = 11, height = 12, dpi = 300)

  # -----------------------------------------
  # paired comparison between arms (NES > 0 in BOTH arms)
  # -----------------------------------------
  nes_paired <- x %>%
    dplyr::mutate(
      arm_short = dplyr::case_when(
        arm == "Resistance arm" ~ "resistance",
        arm == "Sensitisation arm" ~ "sensitisation",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(arm_short)) %>%
    dplyr::select(ID, Description, collection, setSize, leading_edge_n, arm_short, NES, pvalue, p.adjust) %>%
    tidyr::pivot_wider(
      names_from = arm_short,
      values_from = c(NES, pvalue, p.adjust, leading_edge_n),
      values_fn = dplyr::first
    )

  required_cols <- c("NES_resistance", "NES_sensitisation")
  if (!all(required_cols %in% colnames(nes_paired))) {
    nes_paired <- tibble::tibble()
  } else {
    nes_paired <- nes_paired %>%
      dplyr::filter(!is.na(NES_resistance), !is.na(NES_sensitisation), NES_resistance > 0, NES_sensitisation > 0)
  }

  safe_write_csv(nes_paired, file.path(tables_dir, paste0("NES_paired_between_arms_", label, ".csv")))

  if (nrow(nes_paired) > 0) {
    if (sort_col == "p.adjust") {
      nes_paired_plot <- nes_paired %>%
        dplyr::mutate(sort_metric = p.adjust_resistance + p.adjust_sensitisation) %>%
        dplyr::group_by(collection) %>%
        dplyr::arrange(sort_metric, .by_group = TRUE) %>%
        dplyr::slice_head(n = 20) %>%
        dplyr::ungroup()
    } else {
      nes_paired_plot <- nes_paired %>%
        dplyr::mutate(sort_metric = pvalue_resistance + pvalue_sensitisation) %>%
        dplyr::group_by(collection) %>%
        dplyr::arrange(sort_metric, .by_group = TRUE) %>%
        dplyr::slice_head(n = 20) %>%
        dplyr::ungroup()
    }

    p_scatter <- ggplot(
      nes_paired_plot,
      aes(x = NES_sensitisation, y = NES_resistance, size = setSize)
    ) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey60") +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey60") +
      geom_point() +
      facet_wrap(~ collection, scales = "free") +
      labs(
        title = paste0(
          "Pathways detected in both arms (", label, ")\n",
          "Filtered for leading-edge genes >= ", min_leading_edge_n
        ),
        x = "NES (Sensitisation arm)", y = "NES (Resistance arm)", size = "Gene set size"
      ) +
      theme_bw(base_size = 11) +
      theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

    safe_ggsave(file.path(figures_dir, paste0("NES_scatter_pos_vs_neg_", label, ".png")), p_scatter, width = 10, height = 8, dpi = 300)
    safe_ggsave(file.path(figures_dir, paste0("NES_scatter_pos_vs_neg_", label, ".pdf")), p_scatter, width = 10, height = 8, dpi = 300)
  }

  invisible(list(filtered_table = x, nes_compare = nes_compare, nes_paired = nes_paired))
}

run_exploration_for_mode(summary_tbl_small, filter_mode = "padj_0.3", top_n = 15)
run_exploration_for_mode(summary_tbl_small, filter_mode = "pval_0.05", top_n = 15)

message("Done (supplementary exploration).")
message("Tables written under:  ", tables_dir)
message("Figures written under: ", figures_dir)
