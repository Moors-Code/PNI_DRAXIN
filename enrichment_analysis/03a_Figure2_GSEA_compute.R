# ============================================================
# 03a_Figure2_GSEA_compute.R
#
# GSEA analysis for invaded vs uninvaded nerve DE results.
#
# COMPUTE ONLY: builds ranked gene lists (ENTREZ + SYMBOL) from a DE
# table and runs a broad GSEA suite:
#   - GO (CC, BP, MF), Reactome, WikiPathways, Disease Ontology,
#     Network of Cancer Genes, KEGG  (ENTREZ-based)
#   - a custom cell-state GMT signature set                (SYMBOL-based)
#   - MSigDB collections H, C1-C8                           (ENTREZ-based)
# Saves .rds / .csv results only. No plotting here -- see
# 03b_Figure2_GSEA_plots.R for dotplots and rank QC barplots.
#
# Ranking formula: sign(logFC) * -log10(FDR)
#
# Input : input/Figure2_GSEA/FullDE_List.csv          (columns: Gene, logFC, FDR)
#         input/single_nuclei_PDAC/combined_signatures.gmt
# Output: processed/Figure2_GSEA/GSEA/*.rds, *.csv
#         processed/Figure2_GSEA/GSEA/ranks_*.rds, *.csv
#         processed/Figure2_GSEA/logs/, processed/Figure2_GSEA/qc/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(DOSE)
  library(msigdbr)
  library(fgsea)
})

# ------------------------------------------------------------
# parameters
# ------------------------------------------------------------
seed_value <- 1234

n_perm          <- 100000
min_gs_size     <- 10
max_gs_size     <- 800
p_cutoff        <- 0.05
adjust_method   <- "BH"

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
de_file           <- file.path("input", "Figure2_GSEA", "FullDE_List.csv")
combined_gmt_file <- file.path("input", "single_nuclei_PDAC", "combined_signatures.gmt")

processed_root <- file.path("processed", "Figure2_GSEA")
rds_dir <- file.path(processed_root, "GSEA")
csv_dir <- file.path(processed_root, "GSEA")   # kept flat, same folder as .rds
qc_dir  <- file.path(processed_root, "qc")
log_dir <- file.path(processed_root, "logs")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("GSEA_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  write(msg, file = log_file, append = TRUE)
}

set.seed(seed_value)
log_message("Starting GSEA compute.")
log_message("Seed set to: ", seed_value)
log_message("nPermSimple set to: ", n_perm)

if (!file.exists(de_file)) {
  stop("DE table not found: ", de_file, " -- place FullDE_List.csv in input/Figure2_GSEA/")
}
if (!file.exists(combined_gmt_file)) {
  stop("Custom GMT file not found: ", combined_gmt_file, " -- place combined_signatures.gmt in input/single_nuclei_PDAC/")
}

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------
check_required_columns <- function(df, required_cols = c("Gene", "logFC", "FDR")) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Input file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

make_rank_vector <- function(df, gene_col, logfc_col = "logFC", fdr_col = "FDR") {
  df_clean <- df %>%
    dplyr::filter(!is.na(.data[[gene_col]]), !is.na(.data[[logfc_col]]), !is.na(.data[[fdr_col]])) %>%
    dplyr::filter(.data[[fdr_col]] > 0) %>%
    dplyr::mutate(rank_score = sign(.data[[logfc_col]]) * -log10(.data[[fdr_col]])) %>%
    dplyr::distinct(.data[[gene_col]], .keep_all = TRUE)

  ranks <- df_clean$rank_score
  names(ranks) <- df_clean[[gene_col]]
  sort(ranks, decreasing = TRUE)
}

save_rank_table <- function(ranks, gene_id_name, outfile) {
  rank_df <- tibble(!!gene_id_name := names(ranks), rank_score = as.numeric(ranks))
  readr::write_csv(rank_df, outfile)
}

# Export GSEA result table -- COMPUTE ONLY, no plotting
save_gsea_table <- function(gsea_obj, analysis_name, csv_dir, rds_dir) {
  res <- gsea_obj@result %>% as.data.frame()

  if (nrow(res) == 0) log_message("No enriched terms found for ", analysis_name)
  if (!"Description" %in% colnames(res)) res$Description <- NA_character_
  if (!"core_enrichment" %in% colnames(res)) res$core_enrichment <- NA_character_

  preferred_order <- c(
    "ID", "Description", "setSize", "enrichmentScore", "NES",
    "pvalue", "p.adjust", "qvalues", "rank", "leading_edge", "core_enrichment"
  )
  existing_order <- preferred_order[preferred_order %in% colnames(res)]
  remaining_cols <- setdiff(colnames(res), existing_order)
  res <- res[, c(existing_order, remaining_cols), drop = FALSE]

  readr::write_csv(res, file.path(csv_dir, paste0(analysis_name, "_results.csv")))
  saveRDS(gsea_obj, file = file.path(rds_dir, paste0(analysis_name, ".rds")))

  res
}

# Add readable gene symbols to exported results (data enrichment, not a plot)
#
# GSEA-type results (from gseGO/gsePathway/gseWP/gseDO/gseNCG/gseKEGG/GSEA())
# store their leading-edge genes in a column called `core_enrichment` -- NOT
# `geneID` (that column name belongs to a different result type, over-
# representation analysis e.g. enrichGO/enricher, which none of these are).
#
# NOTE: by the time this runs, `core_enrichment` may ALREADY be gene symbols
# -- setReadable() (called upstream, for analyses with make_readable = TRUE)
# converts ENTREZ -> SYMBOL in place before save_all_outputs() is called.
# Re-running an ENTREZ->SYMBOL bitr() lookup on IDs that are already
# symbols finds zero matches and used to silently produce an all-NA
# gene_symbols column. Auto-detect instead of assuming either way.
add_readable_gene_symbols <- function(result_df) {
  if (!"core_enrichment" %in% colnames(result_df)) {
    result_df$gene_symbols <- NA_character_
    return(result_df)
  }

  unique_ids <- unique(unlist(strsplit(result_df$core_enrichment, "/")))
  unique_ids <- unique_ids[!is.na(unique_ids) & unique_ids != ""]

  if (length(unique_ids) == 0) {
    result_df$gene_symbols <- NA_character_
    return(result_df)
  }

  # ENTREZ IDs are purely numeric strings; gene symbols are not.
  # If core_enrichment is already symbols (e.g. setReadable() already ran,
  # or this is a SYMBOL-based custom GMT analysis), there's nothing left to map.
  looks_like_entrez <- all(grepl("^[0-9]+$", unique_ids))

  if (!looks_like_entrez) {
    result_df$gene_symbols <- result_df$core_enrichment
    return(result_df)
  }

  suppressWarnings({
    id_map <- bitr(unique_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  })

  if (is.null(id_map) || nrow(id_map) == 0) {
    result_df$gene_symbols <- NA_character_
    return(result_df)
  }

  id_lookup <- stats::setNames(id_map$SYMBOL, id_map$ENTREZID)

  result_df$gene_symbols <- vapply(
    strsplit(result_df$core_enrichment, "/"),
    function(id_vec) {
      syms <- id_lookup[id_vec]
      syms <- syms[!is.na(syms)]
      paste(unique(syms), collapse = "/")
    },
    character(1)
  )

  result_df
}

save_all_outputs <- function(gsea_obj, analysis_name, csv_dir, rds_dir) {
  res <- save_gsea_table(gsea_obj, analysis_name, csv_dir, rds_dir)

  if ("core_enrichment" %in% colnames(res)) {
    res <- add_readable_gene_symbols(res)
    readr::write_csv(res, file.path(csv_dir, paste0(analysis_name, "_results.csv")))
  }

  invisible(res)
}

run_gsea_entrez <- function(fun, analysis_name, gene_list, make_readable = TRUE, ...) {
  log_message("Running analysis: ", analysis_name)
  set.seed(seed_value)

  gsea_obj <- fun(
    geneList = gene_list, nPermSimple = n_perm,
    minGSSize = min_gs_size, maxGSSize = max_gs_size,
    pvalueCutoff = p_cutoff, pAdjustMethod = adjust_method,
    verbose = TRUE, eps = 0, seed = TRUE, by = "fgsea", ...
  )

  if (make_readable) gsea_obj <- setReadable(gsea_obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

  save_all_outputs(gsea_obj, analysis_name, csv_dir, rds_dir)
  gsea_obj
}

run_gsea_term2gene <- function(analysis_name, gene_list, term2gene,
                               gene_id_type = c("ENTREZID", "SYMBOL"), make_readable = TRUE) {
  gene_id_type <- match.arg(gene_id_type)
  log_message("Running analysis: ", analysis_name)
  set.seed(seed_value)

  gsea_obj <- GSEA(
    geneList = gene_list, TERM2GENE = term2gene, nPermSimple = n_perm,
    minGSSize = min_gs_size, maxGSSize = max_gs_size,
    pvalueCutoff = p_cutoff, pAdjustMethod = adjust_method,
    eps = 0, seed = TRUE, verbose = TRUE, by = "fgsea"
  )

  if (make_readable && gene_id_type == "ENTREZID") {
    gsea_obj <- setReadable(gsea_obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }

  save_all_outputs(gsea_obj, analysis_name, csv_dir, rds_dir)
  gsea_obj
}

# ------------------------------------------------------------
# read input data
# ------------------------------------------------------------
log_message("Reading input DE table: ", de_file)
background <- readr::read_csv(de_file, show_col_types = FALSE)
check_required_columns(background)
log_message("Number of rows in DE table: ", nrow(background))

# ------------------------------------------------------------
# prepare ranked gene lists
# ------------------------------------------------------------
log_message("Converting gene symbols to ENTREZ IDs.")

background_entrez <- bitr(background$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

background_with_id <- background %>%
  dplyr::left_join(background_entrez, by = c("Gene" = "SYMBOL")) %>%
  dplyr::distinct(Gene, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(ENTREZID))

log_message("Genes with valid ENTREZ IDs: ", nrow(background_with_id))

ranks_entrez <- make_rank_vector(df = background_with_id, gene_col = "ENTREZID")
ranks_symbol <- make_rank_vector(df = background, gene_col = "Gene")

log_message("ENTREZ ranked list length: ", length(ranks_entrez))
log_message("SYMBOL ranked list length: ", length(ranks_symbol))

save_rank_table(ranks_entrez, "ENTREZID", file.path(qc_dir, "ranked_genes_entrez.csv"))
save_rank_table(ranks_symbol, "Gene",     file.path(qc_dir, "ranked_genes_symbol.csv"))
saveRDS(ranks_entrez, file.path(qc_dir, "ranked_genes_entrez.rds"))
saveRDS(ranks_symbol, file.path(qc_dir, "ranked_genes_symbol.rds"))

log_message("Ranked gene lists saved.")

# ------------------------------------------------------------
# GO and pathway enrichment analyses (ENTREZ-based)
# ------------------------------------------------------------
gse_cc <- run_gsea_entrez(gseGO, "gse_go_cc", ranks_entrez, make_readable = TRUE, ont = "CC", OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gse_bp <- run_gsea_entrez(gseGO, "gse_go_bp", ranks_entrez, make_readable = TRUE, ont = "BP", OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gse_mf <- run_gsea_entrez(gseGO, "gse_go_mf", ranks_entrez, make_readable = TRUE, ont = "MF", OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

gse_reactome <- run_gsea_entrez(gsePathway, "gse_reactome", ranks_entrez, make_readable = TRUE, organism = "human")
gse_wikipathways <- run_gsea_entrez(gseWP, "gse_wikipathways", ranks_entrez, make_readable = TRUE, organism = "Homo sapiens")
gse_do  <- run_gsea_entrez(gseDO,  "gse_disease_ontology", ranks_entrez, make_readable = TRUE, organism = "hsa")
gse_ncg <- run_gsea_entrez(gseNCG, "gse_ncg", ranks_entrez, make_readable = TRUE, exponent = 1)
gse_kegg <- run_gsea_entrez(gseKEGG, "gse_kegg", ranks_entrez, make_readable = TRUE)

# ------------------------------------------------------------
# custom GMT-based enrichment analysis (SYMBOL-based)
# ------------------------------------------------------------
log_message("Reading custom GMT file: ", combined_gmt_file)
m_t2g_cells_combined <- read.gmt(combined_gmt_file)

gse_cells_combined <- run_gsea_term2gene(
  "gse_cells_combined", ranks_symbol, m_t2g_cells_combined,
  gene_id_type = "SYMBOL", make_readable = FALSE
)

# ------------------------------------------------------------
# MSigDB collection analyses (ENTREZ-based)
# ------------------------------------------------------------
run_msigdb_collection <- function(collection_name) {
  log_message("Preparing MSigDB collection: ", collection_name)
  
  m0 <- msigdbr(species = "Homo sapiens", collection = collection_name)
  
  # msigdbr renamed entrez_gene -> ncbi_gene in v10
  entrez_col <- dplyr::case_when(
    "entrez_gene" %in% names(m0) ~ "entrez_gene",
    "ncbi_gene"   %in% names(m0) ~ "ncbi_gene",
    TRUE ~ NA_character_
  )
  if (is.na(entrez_col)) {
    stop("msigdbr output has no recognised Entrez column. Columns are:\n",
         paste(names(m0), collapse = ", "))
  }
  
  term2gene <- m0 %>%
    dplyr::transmute(gs_name = gs_name, entrez_gene = as.character(.data[[entrez_col]])) %>%
    dplyr::filter(!is.na(entrez_gene), entrez_gene != "") %>%
    dplyr::distinct()
  
  run_gsea_term2gene(
    paste0("gse_msigdb_", collection_name), ranks_entrez, term2gene,
    gene_id_type = "ENTREZID", make_readable = TRUE
  )
}

gse_msigdb_h  <- run_msigdb_collection("H")
gse_msigdb_c1 <- run_msigdb_collection("C1")
gse_msigdb_c2 <- run_msigdb_collection("C2")
gse_msigdb_c3 <- run_msigdb_collection("C3")
gse_msigdb_c4 <- run_msigdb_collection("C4")
gse_msigdb_c5 <- run_msigdb_collection("C5")
gse_msigdb_c6 <- run_msigdb_collection("C6")
gse_msigdb_c7 <- run_msigdb_collection("C7")
gse_msigdb_c8 <- run_msigdb_collection("C8")

# ------------------------------------------------------------
# reproducibility info
# ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), con = file.path(qc_dir, "sessionInfo.txt"))

params_df <- tibble(
  parameter = c("seed_value", "n_perm", "min_gs_size", "max_gs_size", "p_cutoff",
                "adjust_method", "de_file", "combined_gmt_file"),
  value = c(as.character(seed_value), as.character(n_perm), as.character(min_gs_size),
            as.character(max_gs_size), as.character(p_cutoff), adjust_method,
            de_file, combined_gmt_file)
)
readr::write_csv(params_df, file.path(qc_dir, "run_parameters.csv"))

log_message("GSEA compute completed successfully.")
log_message("Results (.rds/.csv) written to: ", rds_dir)
log_message("QC (ranks, session info, params) written to: ", qc_dir)
message("Done. Run 03b_Figure2_GSEA_plots.R next for dotplots and rank barplots.")
