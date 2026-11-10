# ============================================================
# 00_annotate_gene_symbols.R
#
# Add HGNC gene symbols to edgeR DEG tables (ENSEMBL -> SYMBOL).
# This is a lightweight annotation utility, not part of the GSEA
# ranking (GSEA does its own ENSEMBL -> ENTREZ mapping internally).
#
# Input : input/bulkRNAseq/edgeR_expBlocked_*_ALL.csv
# Output: output/tables/bulkRNAseq/00_gene_symbol_annotation/*_with_symbols.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

in_dir  <- file.path("input", "bulkRNAseq")
out_dir <- file.path("output", "bulkRNAseq", "tables", "00_gene_symbol_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  "edgeR_expBlocked_diluted_vs_ctrl_ALL.csv",
  "edgeR_expBlocked_DRAXIN_vs_ctrl_ALL.csv",
  "edgeR_expBlocked_DRAXIN_vs_diluted_ALL.csv"
)

add_gene_symbols <- function(df, ensembl_col = "gene_id") {
  df %>%
    mutate(
      ensembl_clean = sub("\\..*$", "", .data[[ensembl_col]]),
      gene_symbol = mapIds(
        org.Hs.eg.db,
        keys = ensembl_clean,
        keytype = "ENSEMBL",
        column = "SYMBOL",
        multiVals = "first"
      )
    ) %>%
    relocate(gene_symbol, .after = all_of(ensembl_col))
}

for (f in files) {
  in_path  <- file.path(in_dir, f)
  out_path <- file.path(out_dir, sub("\\.csv$", "_with_symbols.csv", f))

  if (!file.exists(in_path)) {
    warning("Missing input file, skipping: ", in_path)
    next
  }

  df       <- read.csv(in_path)
  df_annot <- add_gene_symbols(df, ensembl_col = "gene_id")
  write.csv(df_annot, out_path, row.names = FALSE)
  message("Wrote: ", out_path)
}
