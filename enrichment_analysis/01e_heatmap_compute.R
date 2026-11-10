# ============================================================
# 01e_heatmap_compute.R
#
# LEADING-EDGE pathway heatmaps across 0nM / 17nM / 68nM (BATCH-CORRECTED)
# MERGED across comparisons:
#   - For each target term, union of leading-edge genes across the 3
#     comparisons (unique genes), with a membership table showing which
#     comparisons contributed each gene.
#
# COMPUTE ONLY: batch-corrects expression, reads leading-edge gene lists
# from 01a's saved GSEA .rds results, maps them onto expr matrix IDs,
# and merges across comparisons per target term. Saves one bundle only.
# No plotting here -- see 01f_heatmap_plots.R.
#
# Input : input/bulkRNAseq/Heatmaps_input/{expr_logCPM_TMM.rds, meta2.csv}
#         processed/bulkRNAseq/GSEA/<run_label>/<comparison>/**.rds (via 01a)
# Output: processed/bulkRNAseq/heatmap/merged_leadingedge_bundle.rds
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(limma)
})

# ------------------------------------------------------------
# 0) paths
# ------------------------------------------------------------
gsea_step_root <- file.path("processed", "bulkRNAseq", "GSEA")
run_label <- readLines(file.path(gsea_step_root, "_LATEST_RUN_LABEL.txt"), n = 1)
comparisons <- c("DRAXIN_vs_ctrl", "diluted_vs_ctrl", "DRAXIN_vs_diluted")

base_dir <- file.path(gsea_step_root, run_label)

expr_path <- file.path("input", "bulkRNAseq", "Heatmaps_input", "expr_logCPM_TMM.rds")
meta_path <- file.path("input", "bulkRNAseq", "Heatmaps_input", "meta2.csv")
stopifnot(file.exists(expr_path), file.exists(meta_path))

processed_dir <- file.path("processed", "bulkRNAseq", "heatmap")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1) load expression + metadata
# ------------------------------------------------------------
expr_mat <- readRDS(expr_path)

meta <- readr::read_csv(meta_path, show_col_types = FALSE) |> as.data.frame()
stopifnot("sample" %in% colnames(meta))
rownames(meta) <- meta$sample
meta$sample <- NULL

stopifnot(all(colnames(expr_mat) %in% rownames(meta)))
meta <- meta[colnames(expr_mat), , drop = FALSE]
stopifnot(all(colnames(expr_mat) == rownames(meta)))

# strip ".1" etc in gene rownames if present
rownames(expr_mat) <- sub("\\.\\d+$", "", rownames(expr_mat))

stopifnot("condition" %in% colnames(meta))
meta$condition <- factor(meta$condition, levels = c("ctrl", "diluted_DRAXIN", "DRAXIN"))
cond_map <- c(ctrl = "0nM", diluted_DRAXIN = "17nM", DRAXIN = "68nM")
meta$condition_nm <- factor(cond_map[as.character(meta$condition)], levels = c("0nM", "17nM", "68nM"))

# ------------------------------------------------------------
# 1b) batch + experiment columns
# ------------------------------------------------------------
batch_var <- "batch3"
if (!batch_var %in% colnames(meta)) {
  stop("Batch column not found in meta: '", batch_var, "'. Available: ", paste(colnames(meta), collapse = ", "))
}

exp_var <- dplyr::case_when(
  "exp" %in% colnames(meta) ~ "exp",
  "experiment" %in% colnames(meta) ~ "experiment",
  "Experiment" %in% colnames(meta) ~ "Experiment",
  "run" %in% colnames(meta) ~ "run",
  TRUE ~ NA_character_
)
if (is.na(exp_var)) {
  meta$exp <- "exp1"
  exp_var <- "exp"
}

meta[[batch_var]] <- as.factor(meta[[batch_var]])
meta[[exp_var]]   <- as.factor(meta[[exp_var]])

# ------------------------------------------------------------
# 1c) batch correction (keeps condition effects)
# ------------------------------------------------------------
expr_mat <- limma::removeBatchEffect(
  expr_mat,
  batch  = meta[[batch_var]],
  design = model.matrix(~ condition_nm, data = meta)
)

# ------------------------------------------------------------
# 2) mapping helpers (core_enrichment IDs -> expr rownames)
# ------------------------------------------------------------
detect_keytype <- function(rownms) {
  if (any(grepl("^ENSG", rownms))) return("ENSEMBL")
  if (all(grepl("^[0-9]+$", rownms))) return("ENTREZID")
  "SYMBOL"
}

map_entrez <- function(entrez) {
  entrez <- unique(as.character(entrez))
  mp <- suppressWarnings(
    AnnotationDbi::select(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", columns = c("SYMBOL", "ENSEMBL"))
  ) |> tibble::as_tibble()

  mp |>
    mutate(ENTREZID = as.character(ENTREZID)) |>
    arrange(ENTREZID, is.na(SYMBOL), is.na(ENSEMBL)) |>
    distinct(ENTREZID, .keep_all = TRUE)
}

symbols_to_expr_ids <- function(expr_mat, symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != ""]

  keytype <- detect_keytype(rownames(expr_mat))

  if (keytype == "SYMBOL") {
    ids <- intersect(symbols, rownames(expr_mat))
    return(list(ids = ids, label_map = setNames(ids, ids)))
  }

  mp <- suppressWarnings(
    AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", columns = c("SYMBOL", "ENSEMBL", "ENTREZID"))
  ) |> tibble::as_tibble() |> distinct(SYMBOL, .keep_all = TRUE)

  if (keytype == "ENSEMBL") {
    ids <- mp$ENSEMBL
    lbl <- mp$SYMBOL; names(lbl) <- mp$ENSEMBL
  } else if (keytype == "ENTREZID") {
    ids <- as.character(mp$ENTREZID)
    lbl <- mp$SYMBOL; names(lbl) <- as.character(mp$ENTREZID)
  } else {
    stop("Unexpected keytype: ", keytype)
  }

  ids <- ids[!is.na(ids) & ids != ""]
  ids <- intersect(ids, rownames(expr_mat))

  lbl <- lbl[ids]
  lbl[is.na(lbl) | lbl == ""] <- ids[is.na(lbl) | lbl == ""]
  list(ids = ids, label_map = lbl)
}

entrez_to_expr_ids <- function(expr_mat, entrez) {
  entrez <- unique(as.character(entrez))
  entrez <- entrez[!is.na(entrez) & entrez != ""]
  keytype <- detect_keytype(rownames(expr_mat))

  mp <- map_entrez(entrez)

  if (keytype == "ENTREZID") {
    ids <- intersect(mp$ENTREZID, rownames(expr_mat))
    lbl <- mp$SYMBOL; names(lbl) <- mp$ENTREZID
    return(list(ids = ids, label_map = lbl[ids]))
  }
  if (keytype == "ENSEMBL") {
    ids <- intersect(mp$ENSEMBL, rownames(expr_mat))
    lbl <- mp$SYMBOL; names(lbl) <- mp$ENSEMBL
    return(list(ids = ids, label_map = lbl[ids]))
  }
  ids <- intersect(mp$SYMBOL, rownames(expr_mat))
  lbl <- mp$SYMBOL; names(lbl) <- mp$SYMBOL
  list(ids = ids, label_map = lbl[ids])
}

guess_id_space <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return("UNKNOWN")
  if (all(grepl("^[0-9]+$", x))) return("ENTREZID")
  if (any(grepl("^ENSG", x))) return("ENSEMBL")
  "SYMBOL"
}

map_le_to_expr <- function(expr_mat, le_genes) {
  le_genes <- unique(as.character(le_genes))
  le_genes <- le_genes[!is.na(le_genes) & le_genes != ""]
  if (length(le_genes) == 0) return(list(ids = character(0), label_map = character(0)))

  space <- guess_id_space(le_genes)
  if (space == "ENTREZID") return(entrez_to_expr_ids(expr_mat, le_genes))
  if (space == "SYMBOL")   return(symbols_to_expr_ids(expr_mat, le_genes))

  keytype <- detect_keytype(rownames(expr_mat))
  mp <- suppressWarnings(
    AnnotationDbi::select(org.Hs.eg.db, keys = le_genes, keytype = "ENSEMBL", columns = c("SYMBOL", "ENSEMBL", "ENTREZID"))
  ) |> tibble::as_tibble() |> distinct(ENSEMBL, .keep_all = TRUE)

  if (keytype == "ENSEMBL") {
    ids <- intersect(mp$ENSEMBL, rownames(expr_mat))
    lbl <- mp$SYMBOL; names(lbl) <- mp$ENSEMBL
    return(list(ids = ids, label_map = lbl[ids]))
  } else if (keytype == "ENTREZID") {
    ids <- intersect(as.character(mp$ENTREZID), rownames(expr_mat))
    lbl <- mp$SYMBOL; names(lbl) <- as.character(mp$ENTREZID)
    return(list(ids = ids, label_map = lbl[ids]))
  } else {
    ids <- intersect(mp$SYMBOL, rownames(expr_mat))
    lbl <- mp$SYMBOL; names(lbl) <- mp$SYMBOL
    return(list(ids = ids, label_map = lbl[ids]))
  }
}

# ------------------------------------------------------------
# 3) leading-edge extraction from a saved GSEA result object
# ------------------------------------------------------------
get_leading_edge <- function(gsea_obj, term_id) {
  res <- gsea_obj@result |> as.data.frame()
  hit <- res[res$ID == term_id, , drop = FALSE]
  if (nrow(hit) == 0) return(NULL)  # term not in this comparison's result
  core <- as.character(hit$core_enrichment[1])
  if (is.na(core) || core == "") return(character(0))
  strsplit(core, "/", fixed = TRUE)[[1]]
}

# ------------------------------------------------------------
# 4) RDS path builders (point at 01a's saved GSEA results)
# ------------------------------------------------------------
rds_go_bp <- function(comp) file.path(base_dir, comp, "GO", paste0(run_label, "__", comp, "__gseGO_BP.rds"))
rds_go_mf <- function(comp) file.path(base_dir, comp, "GO", paste0(run_label, "__", comp, "__gseGO_MF.rds"))
rds_kegg  <- function(comp) file.path(base_dir, comp, "KEGG", paste0(run_label, "__", comp, "__gseKEGG.rds"))
rds_hallmark <- function(comp) {
  p1 <- file.path(base_dir, comp, "MSigDB", paste0(run_label, "__", comp, "__msig_H.rds"))
  p2 <- file.path(base_dir, comp, "MSigDB_ALL", "H", paste0(run_label, "__", comp, "__msig_H.rds"))
  if (file.exists(p1)) return(p1)
  if (file.exists(p2)) return(p2)
  return(NA_character_)
}

# ------------------------------------------------------------
# 5) target terms (one heatmap per term, merged across comparisons)
# ------------------------------------------------------------
targets <- list(
  list(collection = "HALLMARK", rds_fun = rds_hallmark, term_id = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", label = "HM_OXPHOS_LE"),
  list(collection = "HALLMARK", rds_fun = rds_hallmark, term_id = "HALLMARK_G2M_CHECKPOINT", label = "HM_G2M_LE"),
  list(collection = "HALLMARK", rds_fun = rds_hallmark, term_id = "HALLMARK_DNA_REPAIR", label = "HM_DNA_REPAIR_LE"),
  list(collection = "HALLMARK", rds_fun = rds_hallmark, term_id = "HALLMARK_CHOLESTEROL_HOMEOSTASIS", label = "HM_CHOLESTEROL_HOMEOSTASIS_LE"),

  list(collection = "GO_MF", rds_fun = rds_go_mf, term_id = "GO:0003954", label = "GOMF_NADH_DEHYDROGENASE_ACTIVITY_LE"),
  list(collection = "GO_MF", rds_fun = rds_go_mf, term_id = "GO:0005506", label = "GOMF_IRON_ION_BINDING_LE"),
  list(collection = "GO_MF", rds_fun = rds_go_mf, term_id = "GO:0015453", label = "GOMF_OXIDOREDUCTION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY_LE"),
  list(collection = "GO_MF", rds_fun = rds_go_mf, term_id = "GO:0050839", label = "GOMF_CELL_ADHESION_MOLECULE_BINDING"),

  list(collection = "GO_BP", rds_fun = rds_go_bp, term_id = "GO:0033627", label = "GOBP_INTEGRIN_MEDIATED_CELL_ADHESION_LE"),

  list(collection = "KEGG", rds_fun = rds_kegg, term_id = "hsa04216", label = "KEGG_FERROPTOSIS_LE")
)

# ------------------------------------------------------------
# 6) merge leading-edge genes across comparisons (union + membership)
# ------------------------------------------------------------
merge_leadingedge_across_comparisons <- function(expr_mat, comparisons, tgt, require_present_in_at_least = 1) {
  mem_list <- list()
  label_map_all <- c()

  for (comp in comparisons) {
    rds_path <- tgt$rds_fun(comp)
    if (is.na(rds_path) || !file.exists(rds_path)) {
      mem_list[[comp]] <- character(0)
      next
    }

    g_obj <- readRDS(rds_path)
    le_raw <- get_leading_edge(g_obj, tgt$term_id)
    if (is.null(le_raw)) {
      mem_list[[comp]] <- character(0)
      next
    }

    mp <- map_le_to_expr(expr_mat, le_raw)
    mem_list[[comp]] <- unique(mp$ids)

    if (length(mp$label_map)) {
      label_map_all <- c(label_map_all, mp$label_map)
    }
  }

  all_ids <- unique(unlist(mem_list, use.names = FALSE))
  if (!length(all_ids)) {
    return(list(ids = character(0), label_map = character(0), mem = data.frame()))
  }

  mem_df <- data.frame(expr_id = all_ids, stringsAsFactors = FALSE)
  for (comp in comparisons) {
    mem_df[[comp]] <- mem_df$expr_id %in% mem_list[[comp]]
  }

  mem_df$n_comp <- rowSums(mem_df[, comparisons, drop = FALSE])
  mem_df <- mem_df[mem_df$n_comp >= require_present_in_at_least, , drop = FALSE]
  ids <- mem_df$expr_id

  if (length(label_map_all)) {
    label_map_all <- label_map_all[!is.na(names(label_map_all)) & names(label_map_all) != ""]
    label_map_all <- label_map_all[!duplicated(names(label_map_all))]
    lbl <- label_map_all[ids]
  } else {
    lbl <- setNames(ids, ids)
    lbl <- lbl[ids]
  }
  lbl[is.na(lbl) | lbl == ""] <- ids[is.na(lbl) | lbl == ""]
  list(ids = ids, label_map = lbl, mem = mem_df)
}

# ------------------------------------------------------------
# 7) run merge for every target, save one bundle
# ------------------------------------------------------------
merged_results <- list()

for (tgt in targets) {
  message("Merging leading-edge genes for: ", tgt$label, " | ", tgt$collection)

  merged <- tryCatch(
    merge_leadingedge_across_comparisons(expr_mat, comparisons, tgt, require_present_in_at_least = 1),
    error = function(e) {
      message("  FAILED merge: ", tgt$label, " | ", conditionMessage(e))
      NULL
    }
  )

  merged_results[[tgt$label]] <- merged
}

saveRDS(
  list(
    expr_mat_bc = expr_mat,
    meta = meta,
    batch_var = batch_var,
    exp_var = exp_var,
    comparisons = comparisons,
    targets = targets,
    merged_results = merged_results,
    gsea_run_label = run_label
  ),
  file.path(processed_dir, "merged_leadingedge_bundle.rds")
)

message("Done. Bundle written to: ", file.path(processed_dir, "merged_leadingedge_bundle.rds"))
message("Run 01f_heatmap_plots.R next for the heatmap PDFs.")
