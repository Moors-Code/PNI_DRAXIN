# ============================================================
# 01c_ssGSEA_compute.R
#
# ssGSEA pathway scoring across 0 nM / 17 nM / 68 nM DRAXIN.
#
# COMPUTE ONLY: this script
#   1) loads precomputed GSEA results from the 3 pairwise comparisons
#      (written by 01a_GSEA_compute.R),
#   2) selects full pathway gene sets of interest,
#   3) computes per-sample ssGSEA scores from batch-corrected logCPM data,
#   4) saves a bundle of everything 01d_ssGSEA_plots.R needs (including the
#      raw + standardized pathway score table -- 01d writes that table to
#      output/, since 01c never writes to output/, only processed/).
# No plotting, and no output/ writes, happen here -- see 01d_ssGSEA_plots.R.
#
# STANDARDIZATION (kept for downstream use):
#   score_sd_all  = z-score of raw ssGSEA values across ALL samples
#                   (general visualization; NOT "relative to ctrl")
#   score_sd_ctrl = (sample score - ctrl mean) / ctrl SD
#                   (expresses each sample relative to the ctrl group only)
#
# GENE SET CONTENT:
#   ssGSEA uses FULL pathway gene sets (not leading-edge subsets),
#   restricted to genes present in the processed expression matrix.
#
# Input : input/bulkRNAseq/Heatmaps_input/expr_logCPM_TMM.rds, meta2.csv
#         processed/bulkRNAseq/GSEA/<run_label>/<comparison>/**.rds
# Output: processed/bulkRNAseq/ssGSEA/ssGSEA_scored_data.rds  (bundle for 01d)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(tidyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(limma)
  library(GSVA)
  library(matrixStats)
  library(msigdbr)
  library(GO.db)
})

# ------------------------------------------------------------
# 0) paths + required files
# ------------------------------------------------------------
gsea_step_root <- file.path("processed", "bulkRNAseq", "GSEA")

# Pick up the most recent GSEA run automatically (written by 01a_GSEA_compute.R).
# To pin an older run instead, just hardcode run_label.
run_label <- readLines(file.path(gsea_step_root, "_LATEST_RUN_LABEL.txt"), n = 1)

comparisons <- c("DRAXIN_vs_ctrl", "diluted_vs_ctrl", "DRAXIN_vs_diluted")

gsea_base_dir <- file.path(gsea_step_root, run_label)

processed_dir <- file.path("processed", "bulkRNAseq", "ssGSEA")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path("input", "bulkRNAseq", "Heatmaps_input", "expr_logCPM_TMM.rds")
meta_path <- file.path("input", "bulkRNAseq", "Heatmaps_input", "meta2.csv")
stopifnot(file.exists(expr_path), file.exists(meta_path))

gsea_paths <- lapply(comparisons, function(comp) {
  list(
    GO_BP  = file.path(gsea_base_dir, comp, "GO",     paste0(run_label, "__", comp, "__gseGO_BP.rds")),
    GO_MF  = file.path(gsea_base_dir, comp, "GO",     paste0(run_label, "__", comp, "__gseGO_MF.rds")),
    MSIG_H = file.path(gsea_base_dir, comp, "MSigDB", paste0(run_label, "__", comp, "__msig_H.rds")),
    KEGG   = file.path(gsea_base_dir, comp, "KEGG",   paste0(run_label, "__", comp, "__gseKEGG.rds"))
  )
})
names(gsea_paths) <- comparisons

missing <- unlist(lapply(names(gsea_paths), function(comp) {
  p <- unlist(gsea_paths[[comp]])
  p[!file.exists(p)]
}))
if (length(missing) > 0) {
  stop("Missing GSEA files (run 01a_GSEA_compute.R first):\n", paste(missing, collapse = "\n"))
}

# ------------------------------------------------------------
# 1) load normalized expression matrix and sample metadata
# ------------------------------------------------------------
expr_mat <- readRDS(expr_path)

meta <- readr::read_csv(meta_path, show_col_types = FALSE) |> as.data.frame()
stopifnot("sample" %in% colnames(meta))
rownames(meta) <- meta$sample
meta$sample <- NULL

stopifnot(all(colnames(expr_mat) %in% rownames(meta)))
meta <- meta[colnames(expr_mat), , drop = FALSE]
stopifnot(all(colnames(expr_mat) == rownames(meta)))

# Remove Ensembl version suffix
rownames(expr_mat) <- sub("\\.\\d+$", "", rownames(expr_mat))

# ------------------------------------------------------------
# 2) treatment labels, batch variable, plotting colors
# ------------------------------------------------------------
meta$condition <- factor(meta$condition, levels = c("ctrl", "diluted_DRAXIN", "DRAXIN"))

cond_map <- c(ctrl = "0nM", diluted_DRAXIN = "17nM", DRAXIN = "68nM")
meta$condition_nm <- factor(cond_map[as.character(meta$condition)],
                            levels = c("0nM", "17nM", "68nM"))

batch_var <- "batch3"
if (!batch_var %in% colnames(meta)) {
  stop("batch_var not in meta: ", batch_var,
       "\nAvailable columns: ", paste(colnames(meta), collapse = ", "))
}
meta[[batch_var]] <- factor(meta[[batch_var]],
                            levels = c("exp1_altprep", "exp1_main", "exp2_main"))

cond_cols <- c("0nM" = "#4E79A7", "17nM" = "#59A14F", "68nM" = "#E15759")

# ------------------------------------------------------------
# 3) load GSEA result objects from all comparisons
# ------------------------------------------------------------
read_gsea_multi <- function(paths_by_comp,
                            slot_name = c("GO_BP", "GO_MF", "MSIG_H", "KEGG")) {
  slot_name <- match.arg(slot_name)
  out <- lapply(names(paths_by_comp), function(comp) {
    readRDS(paths_by_comp[[comp]][[slot_name]])
  })
  names(out) <- names(paths_by_comp)
  out
}

gsea_GO_BP_list  <- read_gsea_multi(gsea_paths, "GO_BP")
gsea_GO_MF_list  <- read_gsea_multi(gsea_paths, "GO_MF")
gsea_MSIG_H_list <- read_gsea_multi(gsea_paths, "MSIG_H")
gsea_KEGG_list   <- read_gsea_multi(gsea_paths, "KEGG")

# ------------------------------------------------------------
# 4) helper functions for pathway lookup and gene mapping
# ------------------------------------------------------------
map_entrez_to_symbol <- function(entrez) {
  entrez <- unique(as.character(entrez))
  mp <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez,
    keytype = "ENTREZID",
    columns = "SYMBOL"
  ) |>
    tibble::as_tibble()

  unique(mp$SYMBOL[!is.na(mp$SYMBOL) & mp$SYMBOL != ""])
}

show_gsea_terms_all <- function(gsea_list, pattern = NULL, n = 50) {
  combined <- lapply(names(gsea_list), function(comp) {
    tb <- tibble::as_tibble(gsea_list[[comp]]@result)
    tb$comparison <- comp
    tb
  }) |>
    dplyr::bind_rows()

  keep <- intersect(c("comparison", "ID", "Description", "NES", "p.adjust"),
                    names(combined))
  combined <- combined |>
    dplyr::select(dplyr::all_of(keep))

  if (!is.null(pattern)) {
    combined <- combined |>
      dplyr::filter(stringr::str_detect(
        tolower(paste(ID, Description)),
        tolower(pattern)
      ))
  }

  combined |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n = n)
}

show_kegg_all_with_sizes <- function(gsea_kegg_list, pattern = NULL, n = 50) {
  gs <- gsea_kegg_list[[1]]@geneSets
  gs_n <- tibble::tibble(ID = names(gs), n_genes = vapply(gs, length, integer(1)))

  combined <- lapply(names(gsea_kegg_list), function(comp) {
    tb <- tibble::as_tibble(gsea_kegg_list[[comp]]@result)
    tb$comparison <- comp
    tb
  }) |>
    dplyr::bind_rows()

  out <- combined |>
    dplyr::left_join(gs_n, by = "ID") |>
    dplyr::select(comparison, ID, Description, NES, p.adjust, n_genes) |>
    dplyr::distinct()

  if (!is.null(pattern)) {
    out <- out |>
      dplyr::filter(stringr::str_detect(
        tolower(paste(ID, Description)),
        tolower(pattern)
      ))
  }

  out |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n = n)
}

kegg_symbols_from_gse_fullset_by_id <- function(gsea_kegg_list, kegg_id) {
  gse <- gsea_kegg_list[[1]]
  gs_list <- gse@geneSets

  if (!kegg_id %in% names(gs_list)) {
    stop("KEGG ID not found in geneSets: ", kegg_id,
         "\nExample available IDs: ",
         paste(head(names(gs_list), 5), collapse = ", "))
  }

  entrez <- gs_list[[kegg_id]]
  map_entrez_to_symbol(entrez)
}

show_hallmark_msigdbr <- function(pattern = NULL, n = 50) {
  hm <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") |>
    dplyr::group_by(gs_name, gs_description) |>
    dplyr::summarise(n_genes = dplyr::n_distinct(gene_symbol), .groups = "drop")

  if (!is.null(pattern)) {
    hm <- hm |>
      dplyr::filter(
        stringr::str_detect(tolower(gs_name), tolower(pattern)) |
          stringr::str_detect(tolower(gs_description), tolower(pattern))
      )
  }

  hm |>
    dplyr::arrange(dplyr::desc(n_genes), gs_name) |>
    dplyr::slice_head(n = n)
}

hallmark_symbols <- function(set_name) {
  msigdbr::msigdbr(species = "Homo sapiens", collection = "H") |>
    dplyr::filter(gs_name == set_name) |>
    dplyr::pull(gene_symbol) |>
    unique()
}

show_go_terms <- function(ontology = c("BP", "MF"), pattern = NULL, n = 50) {
  ontology <- match.arg(ontology)

  tt <- AnnotationDbi::select(
    GO.db,
    keys = AnnotationDbi::keys(GO.db, keytype = "GOID"),
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(!is.na(TERM), ONTOLOGY == ontology) |>
    dplyr::distinct(GOID, TERM, ONTOLOGY)

  if (!is.null(pattern)) {
    tt <- tt |>
      dplyr::filter(stringr::str_detect(tolower(TERM), tolower(pattern)))
  }

  tt |>
    dplyr::arrange(TERM) |>
    dplyr::slice_head(n = n)
}

goid_to_symbols <- function(go_id) {
  entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go_id,
    keytype = "GOALL",
    columns = "ENTREZID"
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::pull(ENTREZID) |>
    unique()

  map_entrez_to_symbol(entrez)
}

pick_go_by_text <- function(ontology = c("BP", "MF"), pattern, prefer_short = TRUE) {
  ontology <- match.arg(ontology)
  pat <- tolower(pattern)

  tt <- AnnotationDbi::select(
    GO.db,
    keys = AnnotationDbi::keys(GO.db, keytype = "GOID"),
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(!is.na(TERM), ONTOLOGY == ontology)

  hit <- tt |>
    dplyr::filter(stringr::str_detect(tolower(TERM), pat))

  if (nrow(hit) == 0) {
    stop("No GO ", ontology, " term matched: ", pattern)
  }

  if (prefer_short) {
    hit <- hit |>
      dplyr::mutate(nchar_term = nchar(TERM)) |>
      dplyr::arrange(nchar_term)
  }

  hit[1, , drop = FALSE]
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
# 5) inspect enriched pathways across comparisons (console only)
# ------------------------------------------------------------
cat("\n--- KEGG across ALL comparisons (includes full-set sizes) ---\n")
print(show_kegg_all_with_sizes(gsea_KEGG_list, n = 25), n = 25)

cat("\n--- KEGG search across ALL comparisons: 'ferroptosis' ---\n")
print(show_kegg_all_with_sizes(gsea_KEGG_list, "ferroptosis", n = 25), n = 25)

cat("\n--- GO BP terms enriched across ALL comparisons (top by q) ---\n")
print(show_gsea_terms_all(gsea_GO_BP_list, n = 25), n = 25)

cat("\n--- GO MF terms enriched across ALL comparisons (top by q) ---\n")
print(show_gsea_terms_all(gsea_GO_MF_list, n = 25), n = 25)

cat("\n--- Hallmark terms enriched across ALL comparisons (top by q) ---\n")
print(show_gsea_terms_all(gsea_MSIG_H_list, n = 25), n = 25)

cat("\n--- Hallmark universe (msigdbr) top by size ---\n")
print(show_hallmark_msigdbr(n = 25), n = 25)

cat("\n--- GO BP browse: 'cell adhesion' ---\n")
print(show_go_terms("BP", "cell adhesion", n = 25), n = 25)

# ------------------------------------------------------------
# 6) define final full gene sets to score
# ------------------------------------------------------------
hallmark_terms <- c(
  OXPHOS = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  E2F_TARGETS = "HALLMARK_E2F_TARGETS",
  MYC = "HALLMARK_MYC_TARGETS_V1",
  DNA_REPAIR = "HALLMARK_DNA_REPAIR",
  G2M = "HALLMARK_G2M_CHECKPOINT",
  CHOLESTEROL_HOMEOSTASIS = "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  UNFOLDED = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
)

go_bp_1 <- pick_go_by_text("BP", "cell adhesion mediated by integrin")
go_bp_2 <- pick_go_by_text("BP", "regulation of cell adhesion")
go_bp_3 <- pick_go_by_text("BP", "SULFUR COMPOUND METABOLIC PROCESS")

go_mf_1 <- pick_go_by_text("MF", "oxidoreduction-driven active transmembrane transporter activity")
go_mf_2 <- pick_go_by_text("MF", "NADH dehydrogenase activity")
go_mf_3 <- pick_go_by_text("MF", "iron ion binding")
go_mf_4 <- pick_go_by_text("MF", "cell adhesion molecule binding")

gene_sets <- list()

for (nm in names(hallmark_terms)) {
  gene_sets[[paste0("HALLMARK__", nm)]] <- hallmark_symbols(hallmark_terms[[nm]])
}

gene_sets[[paste0("GOBP__", go_bp_1$GOID, "__", make.names(go_bp_1$TERM))]] <- goid_to_symbols(go_bp_1$GOID)
gene_sets[[paste0("GOBP__", go_bp_2$GOID, "__", make.names(go_bp_2$TERM))]] <- goid_to_symbols(go_bp_2$GOID)
gene_sets[[paste0("GOBP__", go_bp_3$GOID, "__", make.names(go_bp_3$TERM))]] <- goid_to_symbols(go_bp_3$GOID)
gene_sets[[paste0("GOMF__", go_mf_1$GOID, "__", make.names(go_mf_1$TERM))]] <- goid_to_symbols(go_mf_1$GOID)
gene_sets[[paste0("GOMF__", go_mf_2$GOID, "__", make.names(go_mf_2$TERM))]] <- goid_to_symbols(go_mf_2$GOID)
gene_sets[[paste0("GOMF__", go_mf_3$GOID, "__", make.names(go_mf_3$TERM))]] <- goid_to_symbols(go_mf_3$GOID)
gene_sets[[paste0("GOMF__", go_mf_4$GOID, "__", make.names(go_mf_4$TERM))]] <- goid_to_symbols(go_mf_4$GOID)

gene_sets[["KEGG__hsa04216__Ferroptosis"]] <-
  kegg_symbols_from_gse_fullset_by_id(gsea_KEGG_list, "hsa04216")

names(gene_sets) <- paste0(names(gene_sets), " (n=", sapply(gene_sets, length), ")")

gene_sets <- gene_sets[sapply(gene_sets, length) >= 10]
stopifnot(length(gene_sets) >= 2)

cat("\nSelected gene sets (PRE-intersection sizes):\n")
print(sort(sapply(gene_sets, length)))

# ------------------------------------------------------------
# 7) prepare batch-corrected expression matrix for ssGSEA
# ------------------------------------------------------------
design3 <- model.matrix(~ condition, data = meta)
expr_bc <- limma::removeBatchEffect(expr_mat,
                                    batch = meta[[batch_var]],
                                    design = design3)

ensg <- rownames(expr_bc)
mp <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensg,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) |>
  tibble::as_tibble() |>
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

expr_bc_sym <- expr_bc[mp$ENSEMBL, , drop = FALSE]
rownames(expr_bc_sym) <- mp$SYMBOL

expr_bc_sym <- expr_bc_sym[order(rowMeans(expr_bc_sym), decreasing = TRUE), , drop = FALSE]
expr_bc_sym <- expr_bc_sym[!duplicated(rownames(expr_bc_sym)), , drop = FALSE]

expr_for_gsva <- as.matrix(expr_bc_sym)

v <- matrixStats::rowVars(expr_for_gsva)
expr_for_gsva <- expr_for_gsva[v > 0, , drop = FALSE]

gene_sets_use <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr_for_gsva)))
gene_sets_use <- gene_sets_use[sapply(gene_sets_use, length) >= 10]
stopifnot(length(gene_sets_use) >= 2)

cat("\nSelected gene sets (POST-intersection sizes actually used in ssGSEA):\n")
print(sort(sapply(gene_sets_use, length)))

# ------------------------------------------------------------
# 8) compute raw ssGSEA scores
# ------------------------------------------------------------
gsva_mat <- GSVA::gsva(
  expr = expr_for_gsva,
  gset.idx.list = gene_sets_use,
  method = "ssgsea",
  kcdf = "Gaussian",
  abs.ranking = FALSE,
  min.sz = 10,
  max.sz = 5000
)

scores_raw <- as.data.frame(t(gsva_mat))
scores_raw$sample <- rownames(scores_raw)
scores_raw$condition_nm <- meta[scores_raw$sample, "condition_nm"]

df_long_raw <- scores_raw |>
  tidyr::pivot_longer(
    cols = -c(sample, condition_nm),
    names_to = "pathway",
    values_to = "ssgsea_score"
  ) |>
  dplyr::mutate(
    condition_nm = factor(condition_nm, levels = c("0nM", "17nM", "68nM")),
    pathway_label = make_pathway_label(pathway, wrap_width = 34)
  )

# ------------------------------------------------------------
# 9) two standardized versions of the raw ssGSEA scores
#
# score_sd_all:  z-score across ALL samples within each pathway
# score_sd_ctrl: (sample score - ctrl mean) / ctrl SD within each pathway
#                ctrl samples therefore center around 0 by definition
# ------------------------------------------------------------
sd_floor <- 1e-6

df_long <- df_long_raw |>
  dplyr::group_by(pathway) |>
  dplyr::mutate(
    score_sd_all = as.numeric(scale(ssgsea_score)),
    ctrl_mean = mean(ssgsea_score[condition_nm == "0nM"], na.rm = TRUE),
    ctrl_sd = sd(ssgsea_score[condition_nm == "0nM"], na.rm = TRUE),
    ctrl_sd = ifelse(is.na(ctrl_sd) | ctrl_sd < sd_floor, NA_real_, ctrl_sd),
    score_sd_ctrl = (ssgsea_score - ctrl_mean) / ctrl_sd
  ) |>
  dplyr::ungroup()

# NOTE: the raw + standardized score table (PATHWAY_SCORES_GSVA_ssGSEA_raw_allSD_and_ctrlSD.csv)
# is written by 01d_ssGSEA_plots.R, not here -- 01c never writes to output/, only processed/.

# ------------------------------------------------------------
# save everything 01d_ssGSEA_plots.R needs
# ------------------------------------------------------------
saveRDS(
  list(
    df_long = df_long,
    meta = meta,
    cond_cols = cond_cols,
    gene_sets_use = gene_sets_use,
    expr_for_gsva = expr_for_gsva,
    gsva_mat = gsva_mat,
    gsea_run_label = run_label
  ),
  file.path(processed_dir, "ssGSEA_scored_data.rds")
)

message("ssGSEA compute done.")
message("Scored data bundle: ", file.path(processed_dir, "ssGSEA_scored_data.rds"))
message("Run 01d_ssGSEA_plots.R next to write the score table and all figures to output/.")
