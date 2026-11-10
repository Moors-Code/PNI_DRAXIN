# ============================================================
# 03b_Figure2_GSEA_plots.R
#
# READ-AND-PLOT ONLY: loads every GSEA result .rds written by
# 03a_Figure2_GSEA_compute.R and produces a dotplot for each, plus
# regenerates the two rank-list QC barplots (ENTREZ, SYMBOL).
#
# Input : processed/Figure2_GSEA/GSEA/*.rds
#         processed/Figure2_GSEA/qc/ranked_genes_{entrez,symbol}.rds
# Output: output/Figure2_GSEA/figures/*_dotplot.pdf
#         output/Figure2_GSEA/figures/ranked_genes_{entrez,symbol}_barplot.pdf
# ============================================================

suppressPackageStartupMessages({
  library(enrichplot)
  library(ggplot2)
  library(readr)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
processed_root <- file.path("processed", "Figure2_GSEA")
rds_dir <- file.path(processed_root, "GSEA")
qc_dir  <- file.path(processed_root, "qc")

figures_dir <- file.path("output", "Figure2_GSEA", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

show_categories <- 30

# ------------------------------------------------------------
# dotplot per GSEA result
# ------------------------------------------------------------
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

for (f in rds_files) {
  analysis_name <- sub("\\.rds$", "", basename(f))
  gsea_obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(gsea_obj)) next

  res <- tryCatch(gsea_obj@result, error = function(e) NULL)
  pdf_file <- file.path(figures_dir, paste0(analysis_name, "_dotplot.pdf"))

  pdf(pdf_file, width = 10, height = 8, useDingbats = FALSE)
  if (!is.null(res) && nrow(res) > 0) {
    p <- tryCatch(
      dotplot(gsea_obj, x = "NES", showCategory = show_categories) +
        ggtitle(analysis_name) +
        theme_bw(base_size = 12),
      error = function(e) NULL
    )
    if (!is.null(p)) print(p)
  } else {
    plot.new()
    title(main = paste("No enriched terms for", analysis_name))
    text(0.5, 0.5, "No terms passed the enrichment threshold.")
  }
  dev.off()
}

# ------------------------------------------------------------
# rank QC barplots
# ------------------------------------------------------------
save_rank_barplot <- function(ranks, file, main_title) {
  pdf(file, width = 12, height = 6, useDingbats = FALSE)
  barplot(ranks, border = NA, las = 2, main = main_title, cex.names = 0.4)
  dev.off()
}

# Load ranks from .rds if available (03a writes these), otherwise fall
# back to reconstructing from the .csv (all a migrated pre-reorg run has,
# since the original script only ever saved ranks as .csv).
load_ranks <- function(rds_path, csv_path, id_col) {
  if (file.exists(rds_path)) return(readRDS(rds_path))
  if (file.exists(csv_path)) {
    tab <- readr::read_csv(csv_path, show_col_types = FALSE)
    ranks <- tab$rank_score
    names(ranks) <- tab[[id_col]]
    return(ranks)
  }
  NULL
}

ranks_entrez <- load_ranks(
  file.path(qc_dir, "ranked_genes_entrez.rds"),
  file.path(qc_dir, "ranked_genes_entrez.csv"),
  "ENTREZID"
)
ranks_symbol <- load_ranks(
  file.path(qc_dir, "ranked_genes_symbol.rds"),
  file.path(qc_dir, "ranked_genes_symbol.csv"),
  "Gene"
)

if (!is.null(ranks_entrez)) {
  save_rank_barplot(ranks_entrez, file.path(figures_dir, "ranked_genes_entrez_barplot.pdf"), "Ranked genes (ENTREZ)")
}

if (!is.null(ranks_symbol)) {
  save_rank_barplot(ranks_symbol, file.path(figures_dir, "ranked_genes_symbol_barplot.pdf"), "Ranked genes (SYMBOL)")
}

message("Done. Figures written to: ", figures_dir)
