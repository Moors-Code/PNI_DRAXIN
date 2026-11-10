# ============================================================
# 02b_CRISPR_GSEA_plots.R
#
# READ-AND-PLOT ONLY: this is the PAPER output. Loads the GSEA summary
# table written by 02a_CRISPR_GSEA_compute.R and produces:
#   - the full, unfiltered GSEA results table (all collections, both
#     arms -- no significance/leading-edge filtering), as a
#     supplementary table
#   - the NES comparison plot (Hallmark + C5 collections only,
#     padj <= 0.30) used as the paper's figure
#
# For broader (non-paper) exploration -- the pvalue < 0.05 filter mode,
# per-arm paired scatter plots, per-result dotplots -- see
# CRISPR_GSEA_exploration_supplementary.R.
#
# Input : processed/CRISPR_GSEA/<run_label>/GSEA_summary_compact.csv
# Output: output/CRISPR_GSEA/tables/...
#         output/CRISPR_GSEA/figures/...
# (flat -- no run_label or per-arm subfolders, just the paper's files)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(grid)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
step_root <- file.path("processed", "CRISPR_GSEA")

# Pick up the most recent compute run automatically (written by 02a).
# To pin an older run instead, just hardcode run_label.
run_label <- readLines(file.path(step_root, "_LATEST_RUN_LABEL.txt"), n = 1)

summary_file <- file.path(step_root, run_label, "GSEA_summary_compact.csv")
if (!file.exists(summary_file)) {
  stop("Could not find summary file: ", summary_file, " -- run 02a_CRISPR_GSEA_compute.R first.")
}

tables_dir  <- file.path("output", "CRISPR_GSEA", "tables")
figures_dir <- file.path("output", "CRISPR_GSEA", "figures")
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

summary_tbl_small <- readr::read_csv(summary_file, show_col_types = FALSE)

# ------------------------------------------------------------
# full, unfiltered results table -- every collection, both arms,
# no NES/significance/leading-edge filtering. Supplementary table.
# ------------------------------------------------------------
readr::write_csv(
  summary_tbl_small,
  file.path(tables_dir, "GSEA_all_results.csv")
)

# ------------------------------------------------------------
# settings
# ------------------------------------------------------------
plot_collections <- c("Hallmark", "C5")
top_n_use <- 15
plot_width <- 10
plot_height <- 4
y_text_size <- 12
panel_spacing_y_lines <- 0.05
panel_spacing_x_lines <- 0.08

# ------------------------------------------------------------
# helper for leading edge count
# ------------------------------------------------------------
count_core_genes <- function(core_enrichment_string) {
  if (is.na(core_enrichment_string) || core_enrichment_string == "") return(NA_integer_)
  length(strsplit(core_enrichment_string, "/", fixed = TRUE)[[1]])
}

# ------------------------------------------------------------
# filtered view: Hallmark + C5, NES > 0, setSize >= 10,
# leading_edge_n >= 3, padj <= 0.30
# ------------------------------------------------------------
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
    collection %in% plot_collections,
    !is.na(arm),
    !is.na(NES),
    NES > 0,
    !is.na(setSize),
    setSize >= 10,
    !is.na(leading_edge_n),
    leading_edge_n >= 3,
    !is.na(p.adjust),
    p.adjust <= 0.30
  )

if (nrow(x) == 0) {
  stop("No rows passed the Hallmark/C5 padj <= 0.3 filter.")
}

# ------------------------------------------------------------
# top terms per arm/collection
# ------------------------------------------------------------
nes_compare <- x %>%
  dplyr::group_by(collection, arm) %>%
  dplyr::arrange(p.adjust, dplyr::desc(NES), .by_group = TRUE) %>%
  dplyr::slice_head(n = top_n_use) %>%
  dplyr::ungroup()

readr::write_csv(
  nes_compare,
  file.path(tables_dir, "NES_compare_top_terms_padj_0.3_Hallmark_C5_only.csv")
)

# ------------------------------------------------------------
# plot
# ------------------------------------------------------------
p_nes <- ggplot(
  nes_compare,
  aes(x = reorder(Description, NES), y = NES, size = setSize)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey60") +
  geom_point() +
  coord_flip() +
  facet_grid(collection ~ arm, scales = "free_y", space = "free_y") +
  labs(
    title = "Top enriched pathways by selection arm\npadj <= 0.3, Hallmark + C5 only",
    x = NULL,
    y = "Normalized Enrichment Score (NES)",
    size = "Gene set size"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = y_text_size),
    panel.spacing.y = unit(panel_spacing_y_lines, "lines"),
    panel.spacing.x = unit(panel_spacing_x_lines, "lines")
  )

ggsave(
  file.path(figures_dir, "NES_compare_top_terms_padj_0.3_Hallmark_C5_only.png"),
  p_nes, width = plot_width, height = plot_height, dpi = 300
)

ggsave(
  file.path(figures_dir, "NES_compare_top_terms_padj_0.3_Hallmark_C5_only.pdf"),
  p_nes, width = plot_width, height = plot_height, dpi = 300
)

message("Done. Tables written to:  ", tables_dir, " (GSEA_all_results.csv + filtered NES table)")
message("Figure written to: ", figures_dir)

