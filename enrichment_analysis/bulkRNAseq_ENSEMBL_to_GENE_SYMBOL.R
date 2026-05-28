library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

files <- c(
  "input lists/bulkRNAseq/edgeR_expBlocked_diluted_vs_ctrl_ALL.csv",
  "input lists/bulkRNAseq/edgeR_expBlocked_DRAXIN_vs_ctrl_ALL.csv",
  "input lists/bulkRNAseq/edgeR_expBlocked_DRAXIN_vs_diluted_ALL.csv"
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

for (file in files) {
  df <- read.csv(file)
  
  df_annot <- add_gene_symbols(df, ensembl_col = "gene_id")
  
  out_file <- sub("\\.csv$", "_with_symbols.csv", file)
  write.csv(df_annot, out_file, row.names = FALSE)
}
