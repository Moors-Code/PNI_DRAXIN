# ============================================================
# 01a_GSEA_compute.R
#
# edgeR DEG -> signed-logP ranked gene lists -> GSEA
#   (GO x3 ontologies, Reactome, DO, KEGG, MSigDB Hallmark,
#    then an extra pass across ALL MSigDB collections)
#
# COMPUTE ONLY: this script saves .rds / .csv results and does
# NOT make any plots. See 01b_GSEA_plots.R for figures, which
# reads the .rds objects written here.
#
# Input : input/bulkRNAseq/edgeR_expBlocked_*_ALL.csv
# Output: processed/bulkRNAseq/GSEA/<run_label>/<comparison>/...
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(DOSE)
  library(msigdbr)
})

set.seed(1234)

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
run_date  <- format(Sys.Date(), "%Y_%m_%d")
run_label <- paste0(run_date, "_GSEA_edgeR_3comparisons_signedlogP")

in_dir     <- file.path("input", "bulkRNAseq")
step_root  <- file.path("processed", "bulkRNAseq", "GSEA")
out_dir    <- file.path(step_root, run_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Record which run is "latest" so 01b_GSEA_plots.R (and downstream
# scripts like ssGSEA/heatmap) can find these results without the
# date being hardcoded everywhere.
writeLines(run_label, file.path(step_root, "_LATEST_RUN_LABEL.txt"))

deg_files <- list(
  DRAXIN_vs_ctrl    = file.path(in_dir, "edgeR_expBlocked_DRAXIN_vs_ctrl_ALL.csv"),
  diluted_vs_ctrl   = file.path(in_dir, "edgeR_expBlocked_diluted_vs_ctrl_ALL.csv"),
  DRAXIN_vs_diluted = file.path(in_dir, "edgeR_expBlocked_DRAXIN_vs_diluted_ALL.csv")
)

missing_files <- names(deg_files)[!file.exists(unlist(deg_files))]
if (length(missing_files) > 0) {
  warning("Some DEG files are missing:\n",
          paste0(" - ", missing_files, ": ", unlist(deg_files)[missing_files], collapse = "\n"))
}

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------
guess_col <- function(df, candidates, what) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) stop("Could not find ", what, " column. Tried: ", paste(candidates, collapse = ", "))
  hit
}

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(readr::write_csv(x, path), error = function(e) message("WARN: write_csv failed: ", path, " | ", e$message))
  invisible(TRUE)
}

safe_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(saveRDS(x, path), error = function(e) message("WARN: saveRDS failed: ", path, " | ", e$message))
  invisible(TRUE)
}

# ------------------------------------------------------------
# mapping + ranks
# ------------------------------------------------------------
strip_ensg_version <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^gene:", "", x)
  x <- sub("\\|.*$", "", x)
  sub("\\.\\d+$", "", x)
}

map_ensg_to_entrez <- function(ensg_ids, OrgDb = org.Hs.eg.db) {
  ensg_ids <- unique(stats::na.omit(strip_ensg_version(ensg_ids)))
  if (length(ensg_ids) == 0) return(tibble::tibble(ENSEMBL = character(), ENTREZID = character()))

  AnnotationDbi::select(
    OrgDb, keys = ensg_ids, keytype = "ENSEMBL", columns = c("ENTREZID")
  ) |>
    tibble::as_tibble() |>
    dplyr::transmute(ENSEMBL = .data$ENSEMBL, ENTREZID = as.character(.data$ENTREZID)) |>
    dplyr::filter(!is.na(.data$ENTREZID)) |>
    dplyr::distinct(.data$ENSEMBL, .keep_all = TRUE)
}

make_ranks_entrez_signedlogp <- function(df, label, out_dir_local = NULL, use_fdr = FALSE, p_floor = 1e-300) {

  gene_col <- guess_col(df, c("ENSEMBL", "x", "X", "...1", "gene", "Gene", "ENSG", "gene_id", "Geneid", "id"), "ENSG/ENSEMBL")
  lfc_col  <- guess_col(df, c("logFC", "log2FoldChange", "log2FC"), "logFC")
  p_col    <- if (use_fdr) {
    guess_col(df, c("FDR", "padj", "adj.P.Val", "qvalue"), "FDR/adj p")
  } else {
    guess_col(df, c("PValue", "P.Value", "pvalue", "p_value", "p", "P"), "PValue")
  }

  df2 <- df |>
    dplyr::transmute(
      ENSEMBL = strip_ensg_version(.data[[gene_col]]),
      logFC   = suppressWarnings(as.numeric(.data[[lfc_col]])),
      p       = suppressWarnings(as.numeric(.data[[p_col]]))
    ) |>
    dplyr::filter(!is.na(.data$ENSEMBL), !is.na(.data$logFC), !is.na(.data$p)) |>
    dplyr::mutate(p = pmax(.data$p, p_floor)) |>
    dplyr::group_by(.data$ENSEMBL) |>
    dplyr::slice_min(order_by = .data$p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(score = sign(.data$logFC) * -log10(.data$p))

  map <- map_ensg_to_entrez(df2$ENSEMBL)

  msg <- sprintf("[%s] ENSEMBL rows=%d | mapped=%d | map_rate=%.1f%%",
                 label, nrow(df2), nrow(map), 100 * (nrow(map) / max(1, nrow(df2))))
  message(msg)

  if (!is.null(out_dir_local)) {
    tryCatch(readr::write_lines(msg, file.path(out_dir_local, "MAPPING_RATE.txt")),
             error = function(e) message("WARN: write_lines failed: ", e$message))
    safe_write_csv(map, file.path(out_dir_local, "ENSEMBL_to_ENTREZ_mapping.csv"))
  }

  ranks <- df2$score[match(map$ENSEMBL, df2$ENSEMBL)]
  names(ranks) <- as.character(map$ENTREZID)

  ranks <- ranks[!is.na(ranks) & !is.na(names(ranks))]
  ranks <- ranks[is.finite(ranks)]
  names(ranks) <- as.character(names(ranks))
  ranks <- ranks[!duplicated(names(ranks))]

  sort(ranks, decreasing = TRUE)
}

# ------------------------------------------------------------
# run_and_save_gsea (COMPUTE ONLY -- no plotting)
# ------------------------------------------------------------
run_and_save_gsea <- function(fun, prefix, out_dir_local, readable = TRUE, ...) {
  dir.create(out_dir_local, recursive = TRUE, showWarnings = FALSE)

  status <- tibble::tibble(prefix = prefix, ok = FALSE, n_terms = NA_integer_, error = NA_character_)

  res <- tryCatch(
    suppressMessages(suppressWarnings(fun(...))),
    error = function(e) { status$error <<- e$message; NULL }
  )

  if (is.null(res)) return(list(res = NULL, status = status))

  safe_save_rds(res, file.path(out_dir_local, paste0(prefix, ".rds")))

  res_tbl <- tryCatch(tibble::as_tibble(res@result), error = function(e) tibble::tibble())
  status$ok <- TRUE
  status$n_terms <- nrow(res_tbl)
  safe_write_csv(res_tbl, file.path(out_dir_local, paste0(prefix, ".csv")))

  if (isTRUE(readable) && inherits(res, "gseaResult")) {
    res_readable <- tryCatch(
      clusterProfiler::setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) NULL
    )
    if (!is.null(res_readable)) safe_save_rds(res_readable, file.path(out_dir_local, paste0(prefix, "__readable.rds")))
  }

  list(res = res, status = status)
}

# ------------------------------------------------------------
# MSigDB TERM2GENE builder (version-proof)
# ------------------------------------------------------------
msig_term2gene <- function(collection, subcollection = NULL, species = "Homo sapiens") {

  m0 <- msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection)

  entrez_col <- dplyr::case_when(
    "entrez_gene" %in% names(m0) ~ "entrez_gene",
    "ncbi_gene"   %in% names(m0) ~ "ncbi_gene",
    TRUE ~ NA_character_
  )

  if (is.na(entrez_col)) {
    stop("msigdbr output has no recognized Entrez column. Columns are:\n",
         paste(names(m0), collapse = ", "))
  }

  m0 |>
    dplyr::transmute(term = .data$gs_name, gene = as.character(.data[[entrez_col]])) |>
    dplyr::filter(!is.na(.data$gene), .data$gene != "") |>
    dplyr::distinct()
}

# ------------------------------------------------------------
# GSEA suite (GO x3, Reactome, DO, KEGG, MSigDB Hallmark)
# ------------------------------------------------------------
gsea_params <- list(
  nPermSimple   = 10000,
  minGSSize     = 10,
  maxGSSize     = 800,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = TRUE,
  eps           = 1e-10,
  by            = "fgsea"
)

run_suite_for <- function(ranks, analysis_label, out_dir_local) {
  set.seed(1234)
  pref <- function(x) paste0(run_label, "__", analysis_label, "__", x)
  all_status <- list()

  # GO
  for (ont in c("CC", "BP", "MF")) {
    x <- run_and_save_gsea(
      fun = clusterProfiler::gseGO,
      prefix = pref(paste0("gseGO_", ont)),
      out_dir_local = file.path(out_dir_local, "GO"),
      geneList = ranks,
      ont = ont, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
      nPermSimple = gsea_params$nPermSimple, minGSSize = gsea_params$minGSSize, maxGSSize = gsea_params$maxGSSize,
      pvalueCutoff = gsea_params$pvalueCutoff, pAdjustMethod = gsea_params$pAdjustMethod,
      verbose = gsea_params$verbose, eps = gsea_params$eps, by = gsea_params$by,
      readable = TRUE
    )
    all_status[[length(all_status) + 1]] <- x$status
  }

  # Reactome
  all_status[[length(all_status) + 1]] <- run_and_save_gsea(
    fun = ReactomePA::gsePathway,
    prefix = pref("gseReactome"),
    out_dir_local = file.path(out_dir_local, "Reactome"),
    geneList = ranks,
    organism = "human",
    nPermSimple = gsea_params$nPermSimple, minGSSize = gsea_params$minGSSize, maxGSSize = gsea_params$maxGSSize,
    pvalueCutoff = gsea_params$pvalueCutoff, pAdjustMethod = gsea_params$pAdjustMethod,
    verbose = gsea_params$verbose, eps = gsea_params$eps, by = gsea_params$by,
    readable = TRUE
  )$status

  # DO
  all_status[[length(all_status) + 1]] <- run_and_save_gsea(
    fun = DOSE::gseDO,
    prefix = pref("gseDO"),
    out_dir_local = file.path(out_dir_local, "DO"),
    geneList = ranks,
    organism = "hsa",
    nPermSimple = gsea_params$nPermSimple, minGSSize = gsea_params$minGSSize, maxGSSize = gsea_params$maxGSSize,
    pvalueCutoff = gsea_params$pvalueCutoff, pAdjustMethod = gsea_params$pAdjustMethod,
    verbose = gsea_params$verbose, eps = gsea_params$eps, by = gsea_params$by,
    readable = TRUE
  )$status

  # KEGG (optional; may fail quietly)
  all_status[[length(all_status) + 1]] <- run_and_save_gsea(
    fun = clusterProfiler::gseKEGG,
    prefix = pref("gseKEGG"),
    out_dir_local = file.path(out_dir_local, "KEGG"),
    geneList = ranks,
    organism = "hsa",
    nPermSimple = gsea_params$nPermSimple, minGSSize = gsea_params$minGSSize, maxGSSize = gsea_params$maxGSSize,
    pvalueCutoff = gsea_params$pvalueCutoff, pAdjustMethod = gsea_params$pAdjustMethod,
    verbose = gsea_params$verbose, eps = gsea_params$eps, by = gsea_params$by,
    readable = TRUE
  )$status

  # MSigDB Hallmark only (fast + useful)
  m <- msig_term2gene(collection = "H")

  all_status[[length(all_status) + 1]] <- run_and_save_gsea(
    fun = clusterProfiler::GSEA,
    prefix = pref("msig_H"),
    out_dir_local = file.path(out_dir_local, "MSigDB"),
    geneList = ranks,
    TERM2GENE = m,
    nPermSimple = gsea_params$nPermSimple, minGSSize = gsea_params$minGSSize, maxGSSize = gsea_params$maxGSSize,
    pvalueCutoff = gsea_params$pvalueCutoff, pAdjustMethod = gsea_params$pAdjustMethod,
    eps = gsea_params$eps, by = gsea_params$by,
    verbose = FALSE,
    readable = TRUE
  )$status

  status_df <- dplyr::bind_rows(all_status)
  safe_write_csv(status_df, file.path(out_dir_local, "STATUS.csv"))
  status_df
}

# ------------------------------------------------------------
# main compute loop: ranks + core GSEA suite, per comparison
# ------------------------------------------------------------
all_status <- list()

for (nm in names(deg_files)) {
  message("=== Processing: ", nm, " ===")

  if (!file.exists(deg_files[[nm]])) {
    warning("Missing DEG file for ", nm, ": ", deg_files[[nm]])
    next
  }

  out_comp <- file.path(out_dir, nm)
  dir.create(out_comp, recursive = TRUE, showWarnings = FALSE)

  deg_df <- readr::read_csv(deg_files[[nm]], show_col_types = FALSE)

  if ("x" %in% names(deg_df) && !"ENSEMBL" %in% names(deg_df)) deg_df <- dplyr::rename(deg_df, ENSEMBL = x)
  if ("...1" %in% names(deg_df) && !"ENSEMBL" %in% names(deg_df)) deg_df <- dplyr::rename(deg_df, ENSEMBL = ...1)

  ranks_ent <- make_ranks_entrez_signedlogp(deg_df, nm, out_dir_local = out_comp, use_fdr = FALSE)

  if (length(ranks_ent) < 50) {
    warning("Too few ranked genes after mapping for ", nm, " (n=", length(ranks_ent),
            "). Skipping GSEA. Check mapping file in: ", out_comp)
    next
  }

  safe_save_rds(ranks_ent, file.path(out_comp, "ranks_ENTREZ.rds"))
  safe_write_csv(tibble::tibble(ENTREZID = names(ranks_ent), rank = unname(ranks_ent)),
                 file.path(out_comp, "ranks_ENTREZ.csv"))

  st <- run_suite_for(ranks_ent, analysis_label = nm, out_dir_local = out_comp)
  all_status[[nm]] <- st

  safe_write_csv(st, file.path(out_comp, "STATUS.csv"))
}

status_all <- dplyr::bind_rows(all_status, .id = "comparison")
safe_write_csv(status_all, file.path(out_dir, "STATUS_ALL.csv"))

message("Core GSEA compute done. Outputs written to: ", out_dir)

# ------------------------------------------------------------
# EXTRA PASS: MSigDB-only GSEA for ALL collections/subcollections
# Uses saved ranks_ENTREZ.rds from the loop above.
# ------------------------------------------------------------
msig_collections <- function() c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")

msig_subcollections <- function(collection) {
  if (collection == "C2") return(c("CGP", "CP"))
  if (collection == "C5") return(c("BP", "CC", "MF"))
  return(NULL)
}

safe_msig_term2gene <- function(collection, subcollection = NULL, species = "Homo sapiens") {
  tryCatch(
    msig_term2gene(collection = collection, subcollection = subcollection, species = species),
    error = function(e) {
      message("WARN: msig_term2gene failed for ", collection,
              if (!is.null(subcollection)) paste0(":", subcollection) else "",
              " | ", e$message)
      NULL
    }
  )
}

run_msigdb_only_for <- function(ranks, analysis_label, out_comp,
                                 nPermSimple = gsea_params$nPermSimple,
                                 minGSSize = gsea_params$minGSSize,
                                 maxGSSize = gsea_params$maxGSSize,
                                 pvalueCutoff = gsea_params$pvalueCutoff,
                                 pAdjustMethod = gsea_params$pAdjustMethod,
                                 eps = gsea_params$eps,
                                 by = gsea_params$by) {

  pref <- function(x) paste0(run_label, "__", analysis_label, "__", x)
  out_msig_root <- file.path(out_comp, "MSigDB_ALL")
  dir.create(out_msig_root, recursive = TRUE, showWarnings = FALSE)

  all_status <- list()

  for (coll in msig_collections()) {
    subs <- msig_subcollections(coll)

    if (is.null(subs)) {
      m <- safe_msig_term2gene(collection = coll)
      if (is.null(m) || nrow(m) == 0) next

      st <- run_and_save_gsea(
        fun = clusterProfiler::GSEA,
        prefix = pref(paste0("msig_", coll)),
        out_dir_local = file.path(out_msig_root, coll),
        geneList = ranks,
        TERM2GENE = m,
        nPermSimple = nPermSimple, minGSSize = minGSSize, maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
        eps = eps, by = by,
        verbose = FALSE,
        readable = TRUE
      )$status

      all_status[[length(all_status) + 1]] <- st

    } else {
      for (sub in subs) {
        m <- safe_msig_term2gene(collection = coll, subcollection = sub)
        if (is.null(m) || nrow(m) == 0) next

        st <- run_and_save_gsea(
          fun = clusterProfiler::GSEA,
          prefix = pref(paste0("msig_", coll, "_", sub)),
          out_dir_local = file.path(out_msig_root, coll, sub),
          geneList = ranks,
          TERM2GENE = m,
          nPermSimple = nPermSimple, minGSSize = minGSSize, maxGSSize = maxGSSize,
          pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
          eps = eps, by = by,
          verbose = FALSE,
          readable = TRUE
        )$status

        all_status[[length(all_status) + 1]] <- st
      }
    }
  }

  status_df <- dplyr::bind_rows(all_status)
  safe_write_csv(status_df, file.path(out_msig_root, "STATUS_MSigDB_ALL.csv"))
  status_df
}

msig_all_status <- list()

for (nm in names(deg_files)) {
  out_comp   <- file.path(out_dir, nm)
  ranks_path <- file.path(out_comp, "ranks_ENTREZ.rds")

  if (!file.exists(ranks_path)) {
    warning("No ranks_ENTREZ.rds found for ", nm, " at: ", ranks_path, " | skipping MSigDB-only pass")
    next
  }

  ranks_ent <- readRDS(ranks_path)

  if (length(ranks_ent) < 50) {
    warning("Too few ranks for ", nm, " (n=", length(ranks_ent), ") | skipping MSigDB-only pass")
    next
  }

  message("=== MSigDB-only pass: ", nm, " ===")
  msig_all_status[[nm]] <- run_msigdb_only_for(ranks_ent, analysis_label = nm, out_comp = out_comp)
}

msig_status_all <- dplyr::bind_rows(msig_all_status, .id = "comparison")
safe_write_csv(msig_status_all, file.path(out_dir, "STATUS_MSigDB_ALL_COMPARISONS.csv"))

message("MSigDB-only pass complete. Written under each comparison: <comparison>/MSigDB_ALL/")
message("Run label for this run (used by 01b_GSEA_plots.R): ", run_label)
