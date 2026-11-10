# ============================================================
# 02a_CRISPR_GSEA_compute.R
#
# MAGeCK gene_summary -> p-value-based ranking -> GSEA
#   (Hallmark, C2, C5, Reactome; separate positive/negative selection arms)
#
# COMPUTE ONLY: saves .rds / .csv only. No plotting here -- see
# 02b_CRISPR_GSEA_plots.R for the paper figure, and
# CRISPR_GSEA_exploration_supplementary.R for broader (non-paper)
# exploratory plots.
#
# This is the p-value-based scoring variant (rank_metric = lfc * -log10(p)),
# run on the PLUS_CONTROL_PAIRED_TOTAL MAGeCK filtering scheme -- the
# variant that goes into the paper.
#
# RANKING LOGIC:
#   positive-selection (resistance arm):   rank_metric = pos_lfc * -log10(pos_p)
#   negative-selection (sensitisation arm): rank_metric = (-neg_lfc) * -log10(neg_p)
#   (MAGeCK negative-selection hits often have negative LFC; inverting neg_lfc
#   makes the strongest sensitising hits rank at the TOP of the list for GSEA.)
#
# MULTI-TARGET / NON-GENE DETECTION (fixed 2026-09-10):
#   Guides/genes excluded from GSEA ranking come from THREE separate
#   sources, and all must be checked against the pipeline's OWN
#   gene-naming convention, not the raw annotation file's formatting:
#
#   1) By-design multi-target guides: the official annotation
#      (inzolia_annotation_duplicates_removed.csv) records these with a
#      semicolon-joined Kept_Target, e.g. "KRT17;KRT19". BUT this pipeline's
#      own Gene column (built upstream in 01_build_Library_and_Screen.R via
#      `sub("_[0-9]+$", "", sgRNA)`) never produces a semicolon -- for that
#      same guide it instead produces the underscore-joined "KRT17_KRT19".
#      A naive `grepl(";", gene)` check on gene_summary.txt's gene column
#      can therefore NEVER fire on real data, silently letting ~5,282
#      by-design multi-target gene-pairs/triples through as if they were
#      single-gene hits.
#
#   2) Dedup-driven shared-sequence ambiguity: genes where deduplication of
#      the guide library collapsed two originally-independent single-target
#      guide designs into one because they happened to share an identical
#      140bp sequence (e.g. ABCF2 / ABCF2-H2BE1, a readthrough gene pair).
#      These are recorded via Other_Guides/Other_Targets in the same
#      annotation file (~252 genes), and were not checked for at all before.
#
#   3) Control/non-gene pseudo-genes: the library's negative control guides
#      (EGFP_*, INTERGENIC_SITE_*, NO_SITE_*) also end in "_<number>", so
#      the same naive suffix-stripping collapses ALL individual control
#      sites into 3 pseudo-genes -- "EGFP" (20 sites merged),
#      "INTERGENIC_SITE" (500 sites), "NO_SITE" (50 sites). These aren't
#      real genes; "EGFP" in particular is a real gene symbol in some
#      contexts and could otherwise be wrongly pulled into a genuine
#      pathway enrichment result.
#
#   Fix: build the exclusion set by reconstructing category 1 using the
#   SAME naive suffix-stripping logic the pipeline itself uses (applied to
#   GuideID, not Kept_Target), adding category 2 directly from
#   Other_Targets, and hardcoding category 3's 3 known pseudo-gene names.
#   Verified against real examples: KRT17_KRT19, AAK1_BMP2K,
#   H2AC6_H2AJ_H2AX (3-way), ABCF2/ABCF2-H2BE1, EGFP/INTERGENIC_SITE/NO_SITE.
#
# Excluded entries (per the corrected detection above) are removed from
# GSEA ranking but kept in separate exported tables for reference.
#
# Input : input/CRISPR_GSEA/<mageck gene_summary.txt for PLUS_CONTROL_PAIRED_TOTAL>
#         input/CRISPR_GSEA/inzolia_annotation_duplicates_removed.csv
# Output: processed/CRISPR_GSEA/<run_label>/...
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(clusterProfiler)
  library(msigdbr)
  library(ReactomePA)
  library(org.Hs.eg.db)
})

set.seed(1234)

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
# Point this at whatever you named the MAGeCK gene_summary.txt file you
# placed in input/CRISPR_GSEA/ for the PLUS_CONTROL_PAIRED_TOTAL screen.
mageck_file <- file.path("input", "CRISPR", "T0_PLUS_CONTROL_PAIRED_TOTAL.gene_summary.txt")

if (!file.exists(mageck_file)) {
  stop(
    "MAGeCK gene_summary.txt not found: ", mageck_file,
    "\nPlace the PLUS_CONTROL_PAIRED_TOTAL MAGeCK gene_summary.txt in input/CRISPR_GSEA/",
    " and update `mageck_file` above if you used a different filename."
  )
}

# Fill in: path to inzolia_annotation_duplicates_removed.csv (needed for
# correct multi-target detection, see note above).
dup_annotation_file <- file.path("input", "CRISPR", "inzolia_annotation_duplicates_removed.csv")

if (!file.exists(dup_annotation_file)) {
  stop(
    "Duplicate-guide annotation file not found: ", dup_annotation_file,
    "\nThis file is required to correctly identify multi-target and",
    " dedup-driven shared-sequence-ambiguous genes. Place it in input/CRISPR_GSEA/",
    " and update `dup_annotation_file` above if you used a different filename."
  )
}

run_date  <- format(Sys.Date(), "%Y_%m_%d")
run_label <- paste0(run_date, "_GSEA_MAGeCK_LFCp")

step_root <- file.path("processed", "CRISPR_GSEA")
out_dir   <- file.path(step_root, run_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Record which run is "latest" so 02b_CRISPR_GSEA_plots.R (and the
# supplementary script) can find these results automatically.
writeLines(run_label, file.path(step_root, "_LATEST_RUN_LABEL.txt"))

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------
guess_col <- function(df, candidates, what) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) {
    stop("Could not find ", what, " column. Tried: ", paste(candidates, collapse = ", "))
  }
  hit
}

find_optional <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) return(NA_character_)
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

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  sub("^gene:", "", x)
}

to_num <- function(x) suppressWarnings(as.numeric(x))

pull_optional_num <- function(df, colname) {
  if (is.na(colname)) return(rep(NA_real_, nrow(df)))
  to_num(df[[colname]])
}

safe_result_tbl <- function(res) {
  tryCatch(tibble::as_tibble(res@result), error = function(e) tibble::tibble())
}

# ------------------------------------------------------------
# multi-target detection (fixed 2026-09-10 -- see header note)
# ------------------------------------------------------------
dup_annotation <- readr::read_csv(dup_annotation_file, show_col_types = FALSE)
colnames(dup_annotation)[1] <- "GuideID"

# Reproduce the pipeline's own Gene-naming logic exactly (from
# 01_build_Library_and_Screen.R), so the lookup set matches whatever
# actually ends up in gene_summary.txt's gene/id column.
pipeline_gene_name <- function(guide_id) {
  sub("_[0-9]+$", "", guide_id)
}

# Category 1: by-design multi-target guides (official Kept_Target has a
# semicolon -- i.e. this GuideID's underlying design targets >1 gene).
# Flagged using the PIPELINE'S version of the name (underscore-joined),
# not the official semicolon-joined Kept_Target, since that's what will
# actually appear in gene_summary.txt.
cat1_tbl <- dup_annotation %>%
  filter(grepl(";", Kept_Target, fixed = TRUE)) %>%
  transmute(
    gene = pipeline_gene_name(GuideID),
    official_kept_target = Kept_Target,
    category = "by_design_multi_target"
  ) %>%
  distinct(gene, .keep_all = TRUE)

# Category 2: dedup-driven shared-sequence ambiguity (Other_Targets
# populated). Kept_Target here is already a single plain name with no
# underscore/semicolon mismatch issue, so no reconstruction needed.
cat2_tbl <- dup_annotation %>%
  filter(!is.na(Other_Targets), trimws(Other_Targets) != "") %>%
  transmute(
    gene = Kept_Target,
    official_kept_target = Kept_Target,
    category = "dedup_shared_sequence"
  ) %>%
  distinct(gene, .keep_all = TRUE)

ambiguous_genes_tbl <- bind_rows(cat1_tbl, cat2_tbl) %>%
  distinct(gene, .keep_all = TRUE)

ambiguous_gene_set <- unique(ambiguous_genes_tbl$gene)

message("Category 1 (by-design multi-target, reconstructed to pipeline naming): ", nrow(cat1_tbl))
message("Category 2 (dedup-driven shared-sequence): ", nrow(cat2_tbl))

# Category 3: control/non-gene pseudo-genes. The library's negative
# control guides (EGFP_*, INTERGENIC_SITE_*, NO_SITE_*) end in "_<number>"
# just like real genes, so the same naive suffix-stripping used to build
# the Gene column collapses ALL individual control sites into just 3
# pseudo-genes: "EGFP" (20 sites merged), "INTERGENIC_SITE" (500 sites),
# "NO_SITE" (50 sites). These are not real genes and must be excluded --
# "EGFP" in particular is a real gene symbol in some contexts and could
# otherwise be wrongly pulled into a genuine pathway enrichment result.
control_pseudogene_set <- c("EGFP", "INTERGENIC_SITE", "NO_SITE")

message("Category 3 (control/non-gene pseudo-genes, e.g. pooled EGFP/intergenic/no-site controls): ",
        length(control_pseudogene_set))

excluded_gene_set <- union(ambiguous_gene_set, control_pseudogene_set)

message("Total unique genes excluded from GSEA (multi-target + control pseudo-genes): ", length(excluded_gene_set))

safe_write_csv(ambiguous_genes_tbl, file.path(out_dir, "genes_flagged_multi_target_all_categories.csv"))
safe_write_csv(tibble(gene = control_pseudogene_set, category = "control_pseudogene"),
               file.path(out_dir, "genes_flagged_control_pseudogenes.csv"))

# Single, unified check covering all three categories, using names in the
# SAME format the pipeline actually produces.
is_multi_target <- function(x) as.character(x) %in% excluded_gene_set

# ------------------------------------------------------------
# ranking helpers
# ------------------------------------------------------------
compute_lfc_p_metric <- function(lfc, p, p_floor = 1e-300, flip_sign = FALSE) {
  p_use <- pmax(p, p_floor)
  lfc_use <- if (flip_sign) -lfc else lfc
  lfc_use * (-log10(p_use))
}

dedup_for_lfc_p <- function(df, lfc_col, p_col, p_floor = 1e-300, flip_sign = FALSE) {
  df %>%
    filter(
      !is.na(gene), gene != "",
      !is.na(.data[[lfc_col]]), is.finite(.data[[lfc_col]]),
      !is.na(.data[[p_col]]), is.finite(.data[[p_col]]), .data[[p_col]] > 0
    ) %>%
    mutate(
      rank_metric = compute_lfc_p_metric(.data[[lfc_col]], .data[[p_col]], p_floor = p_floor, flip_sign = flip_sign)
    ) %>%
    group_by(gene) %>%
    slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
    ungroup()
}

make_rank_vector_lfc_p <- function(df, gene_col = "gene", lfc_col, p_col,
                                   p_floor = 1e-300, flip_sign = FALSE) {
  x <- df %>%
    dplyr::select(gene = all_of(gene_col), lfc = all_of(lfc_col), p = all_of(p_col)) %>%
    filter(!is.na(gene), gene != "", !is.na(lfc), is.finite(lfc), !is.na(p), is.finite(p), p > 0) %>%
    mutate(
      p_use = pmax(p, p_floor),
      lfc_use = if (flip_sign) -lfc else lfc,
      rank_score = lfc_use * -log10(p_use)
    ) %>%
    distinct(gene, .keep_all = TRUE)
  
  ranks <- x$rank_score
  names(ranks) <- x$gene
  sort(ranks, decreasing = TRUE)
}

export_rank_table_lfc_p <- function(df, lfc_col, p_col, outfile, p_floor = 1e-300, flip_sign = FALSE) {
  tab <- df %>%
    filter(
      !is.na(gene), gene != "",
      !is.na(.data[[lfc_col]]), is.finite(.data[[lfc_col]]),
      !is.na(.data[[p_col]]), is.finite(.data[[p_col]]), .data[[p_col]] > 0
    ) %>%
    transmute(
      gene = gene,
      lfc  = .data[[lfc_col]],
      p    = .data[[p_col]],
      lfc_used = if (flip_sign) -lfc else lfc,
      rank = compute_lfc_p_metric(lfc, p, p_floor = p_floor, flip_sign = flip_sign)
    ) %>%
    arrange(desc(rank))
  
  safe_write_csv(tab, outfile)
}

# ------------------------------------------------------------
# GSEA wrappers (COMPUTE ONLY -- no plotting)
# ------------------------------------------------------------
run_and_save_gsea_term2gene <- function(ranks, term2gene, prefix, out_dir_local,
                                        minGSSize = 10, maxGSSize = 500,
                                        pvalueCutoff = 1, pAdjustMethod = "BH",
                                        eps = 1e-10, nPermSimple = 10000) {
  dir.create(out_dir_local, recursive = TRUE, showWarnings = FALSE)
  
  res <- tryCatch(
    suppressMessages(suppressWarnings(
      clusterProfiler::GSEA(
        geneList = ranks, TERM2GENE = term2gene,
        minGSSize = minGSSize, maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
        eps = eps, by = "fgsea", nPermSimple = nPermSimple, verbose = FALSE
      )
    )),
    error = function(e) { message("WARN: GSEA failed for ", prefix, " | ", e$message); NULL }
  )
  
  if (is.null(res)) return(NULL)
  
  safe_save_rds(res, file.path(out_dir_local, paste0(prefix, ".rds")))
  safe_write_csv(safe_result_tbl(res), file.path(out_dir_local, paste0(prefix, ".csv")))
  
  res
}

run_and_save_reactome <- function(ranks_entrez, prefix, out_dir_local,
                                  minGSSize = 10, maxGSSize = 500,
                                  pvalueCutoff = 1, pAdjustMethod = "BH",
                                  eps = 1e-10, nPermSimple = 10000) {
  dir.create(out_dir_local, recursive = TRUE, showWarnings = FALSE)
  
  res <- tryCatch(
    suppressMessages(suppressWarnings(
      ReactomePA::gsePathway(
        geneList = ranks_entrez, organism = "human",
        minGSSize = minGSSize, maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
        eps = eps, by = "fgsea", nPermSimple = nPermSimple, verbose = FALSE
      )
    )),
    error = function(e) { message("WARN: Reactome GSEA failed for ", prefix, " | ", e$message); NULL }
  )
  
  if (is.null(res)) return(NULL)
  
  safe_save_rds(res, file.path(out_dir_local, paste0(prefix, ".rds")))
  safe_write_csv(safe_result_tbl(res), file.path(out_dir_local, paste0(prefix, ".csv")))
  
  res
}

# ------------------------------------------------------------
# read MAGeCK gene summary
# ------------------------------------------------------------
df <- readr::read_tsv(mageck_file, show_col_types = FALSE)

message("Loaded rows: ", nrow(df))
message("Columns: ", paste(colnames(df), collapse = " | "))

gene_col <- guess_col(df, c("id", "ID", "Gene", "gene", "Symbol", "symbol"), "gene identifier")

pos_score_col <- find_optional(df, c("pos|score", "pos.score", "pos_score"))
neg_score_col <- find_optional(df, c("neg|score", "neg.score", "neg_score"))

pos_p_col   <- guess_col(df, c("pos|p-value", "pos|pvalue", "pos.p.value", "pos_p_value"), "positive-selection p-value")
pos_lfc_col <- guess_col(df, c("pos|lfc", "pos.lfc", "pos_lfc", "pos|fc", "pos.fc", "pos_fc"), "positive-selection LFC")
neg_p_col   <- guess_col(df, c("neg|p-value", "neg|pvalue", "neg.p.value", "neg_p_value"), "negative-selection p-value")
neg_lfc_col <- guess_col(df, c("neg|lfc", "neg.lfc", "neg_lfc", "neg|fc", "neg.fc", "neg_fc"), "negative-selection LFC")

pos_fdr_col  <- find_optional(df, c("pos|fdr", "pos.fdr", "pos_fdr"))
pos_rank_col <- find_optional(df, c("pos|rank", "pos.rank", "pos_rank"))
neg_fdr_col  <- find_optional(df, c("neg|fdr", "neg.fdr", "neg_fdr"))
neg_rank_col <- find_optional(df, c("neg|rank", "neg.rank", "neg_rank"))

# ------------------------------------------------------------
# prep table
# ------------------------------------------------------------
df2 <- df %>%
  transmute(
    gene      = clean_gene_symbol(.data[[gene_col]]),
    pos_score = pull_optional_num(df, pos_score_col),
    neg_score = pull_optional_num(df, neg_score_col),
    pos_p     = to_num(.data[[pos_p_col]]),
    pos_fdr   = pull_optional_num(df, pos_fdr_col),
    pos_rank  = pull_optional_num(df, pos_rank_col),
    pos_lfc   = to_num(.data[[pos_lfc_col]]),
    neg_p     = to_num(.data[[neg_p_col]]),
    neg_fdr   = pull_optional_num(df, neg_fdr_col),
    neg_rank  = pull_optional_num(df, neg_rank_col),
    neg_lfc   = to_num(.data[[neg_lfc_col]])
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  mutate(target_type = ifelse(is_multi_target(gene), "multi_target", "single_target"))

safe_write_csv(df2, file.path(out_dir, "mageck_all_targets_cleaned.csv"))

df_single <- df2 %>% filter(target_type == "single_target")
df_multi  <- df2 %>% filter(target_type == "multi_target")

safe_write_csv(df_single, file.path(out_dir, "mageck_single_targets.csv"))
safe_write_csv(df_multi,  file.path(out_dir, "mageck_multi_targets_excluded_from_primary_GSEA.csv"))

message("Single-target rows: ", nrow(df_single))
message("Multi-target rows: ", nrow(df_multi))

# ------------------------------------------------------------
# deduplication per arm
# ------------------------------------------------------------
df_pos_dedup <- dedup_for_lfc_p(df_single, lfc_col = "pos_lfc", p_col = "pos_p", flip_sign = FALSE)
df_neg_dedup <- dedup_for_lfc_p(df_single, lfc_col = "neg_lfc", p_col = "neg_p", flip_sign = TRUE)

safe_write_csv(df_pos_dedup, file.path(out_dir, "mageck_single_targets_deduplicated_positive_selection.csv"))
safe_write_csv(df_neg_dedup, file.path(out_dir, "mageck_single_targets_deduplicated_negative_selection.csv"))

# ------------------------------------------------------------
# rank vectors (symbol-based, for MSigDB GSEA)
# ------------------------------------------------------------
ranks_pos <- make_rank_vector_lfc_p(df_pos_dedup, lfc_col = "pos_lfc", p_col = "pos_p", flip_sign = FALSE)
ranks_neg <- make_rank_vector_lfc_p(df_neg_dedup, lfc_col = "neg_lfc", p_col = "neg_p", flip_sign = TRUE)

safe_save_rds(ranks_pos, file.path(out_dir, "ranks_positive_selection_LFCxlog10p.rds"))
safe_save_rds(ranks_neg, file.path(out_dir, "ranks_negative_selection_LFCxlog10p.rds"))

export_rank_table_lfc_p(df_pos_dedup, lfc_col = "pos_lfc", p_col = "pos_p",
                        outfile = file.path(out_dir, "ranks_positive_selection_LFCxlog10p.csv"), flip_sign = FALSE)
export_rank_table_lfc_p(df_neg_dedup, lfc_col = "neg_lfc", p_col = "neg_p",
                        outfile = file.path(out_dir, "ranks_negative_selection_LFCxlog10p.csv"), flip_sign = TRUE)

# ------------------------------------------------------------
# MSigDB TERM2GENE builders (symbol-based)
# ------------------------------------------------------------
msig_term2gene_symbols <- function(collection, subcollection = NULL, species = "Homo sapiens") {
  m0 <- msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection)
  m0 %>%
    transmute(term = gs_name, gene = gene_symbol) %>%
    filter(!is.na(gene), gene != "") %>%
    distinct()
}

hallmark_t2g <- msig_term2gene_symbols(collection = "H")
c2_t2g       <- msig_term2gene_symbols(collection = "C2")
c5_t2g       <- msig_term2gene_symbols(collection = "C5")

safe_write_csv(hallmark_t2g, file.path(out_dir, "TERM2GENE_Hallmark.csv"))
safe_write_csv(c2_t2g,       file.path(out_dir, "TERM2GENE_C2.csv"))
safe_write_csv(c5_t2g,       file.path(out_dir, "TERM2GENE_C5.csv"))

# ------------------------------------------------------------
# map SYMBOL -> ENTREZ for Reactome (ReactomePA expects ENTREZ IDs)
# ------------------------------------------------------------
map_symbol_to_entrez <- function(symbols) {
  map <- tryCatch(
    clusterProfiler::bitr(unique(symbols), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db),
    error = function(e) NULL
  )
  if (is.null(map)) stop("SYMBOL -> ENTREZ mapping failed. Make sure org.Hs.eg.db is installed.")
  map %>% distinct(SYMBOL, .keep_all = TRUE) %>% transmute(gene = SYMBOL, entrez = ENTREZID)
}

sym2ent_pos <- map_symbol_to_entrez(df_pos_dedup$gene)
sym2ent_neg <- map_symbol_to_entrez(df_neg_dedup$gene)

safe_write_csv(sym2ent_pos, file.path(out_dir, "SYMBOL_to_ENTREZ_mapping_positive_selection.csv"))
safe_write_csv(sym2ent_neg, file.path(out_dir, "SYMBOL_to_ENTREZ_mapping_negative_selection.csv"))

make_entrez_ranks <- function(symbol_ranks, mapping_df) {
  rank_df <- tibble(gene = names(symbol_ranks), score = as.numeric(symbol_ranks)) %>%
    left_join(mapping_df, by = "gene") %>%
    filter(!is.na(entrez), !is.na(score), is.finite(score)) %>%
    distinct(entrez, .keep_all = TRUE)
  ranks <- rank_df$score
  names(ranks) <- rank_df$entrez
  sort(ranks, decreasing = TRUE)
}

ranks_pos_entrez <- make_entrez_ranks(ranks_pos, sym2ent_pos)
ranks_neg_entrez <- make_entrez_ranks(ranks_neg, sym2ent_neg)

safe_save_rds(ranks_pos_entrez, file.path(out_dir, "ranks_positive_selection_LFCxlog10p_ENTREZ.rds"))
safe_save_rds(ranks_neg_entrez, file.path(out_dir, "ranks_negative_selection_LFCxlog10p_ENTREZ.rds"))

safe_write_csv(tibble(ENTREZID = names(ranks_pos_entrez), rank = unname(ranks_pos_entrez)),
               file.path(out_dir, "ranks_positive_selection_LFCxlog10p_ENTREZ.csv"))
safe_write_csv(tibble(ENTREZID = names(ranks_neg_entrez), rank = unname(ranks_neg_entrez)),
               file.path(out_dir, "ranks_negative_selection_LFCxlog10p_ENTREZ.csv"))

# ------------------------------------------------------------
# run all GSEA analyses
# ------------------------------------------------------------
analysis_specs <- list(
  list(name = "Hallmark", term2gene = hallmark_t2g),
  list(name = "C2",       term2gene = c2_t2g),
  list(name = "C5",       term2gene = c5_t2g)
)

summary_rows <- list()

for (spec in analysis_specs) {
  coll_name <- spec$name
  coll_t2g  <- spec$term2gene
  coll_dir  <- file.path(out_dir, paste0("MSigDB_", coll_name))
  dir.create(coll_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("=== Running ", coll_name, " | positive-selection ===")
  res_pos <- run_and_save_gsea_term2gene(
    ranks = ranks_pos, term2gene = coll_t2g,
    prefix = paste0(run_label, "__", coll_name, "__positive_selection"),
    out_dir_local = coll_dir
  )
  
  message("=== Running ", coll_name, " | negative-selection ===")
  res_neg <- run_and_save_gsea_term2gene(
    ranks = ranks_neg, term2gene = coll_t2g,
    prefix = paste0(run_label, "__", coll_name, "__negative_selection"),
    out_dir_local = coll_dir
  )
  
  pos_tbl <- safe_result_tbl(res_pos)
  neg_tbl <- safe_result_tbl(res_neg)
  
  if (nrow(pos_tbl) > 0) summary_rows[[length(summary_rows) + 1]] <- pos_tbl %>% mutate(collection = coll_name, arm = "positive_selection")
  if (nrow(neg_tbl) > 0) summary_rows[[length(summary_rows) + 1]] <- neg_tbl %>% mutate(collection = coll_name, arm = "negative_selection")
}

# Reactome
reactome_dir <- file.path(out_dir, "Reactome")
dir.create(reactome_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Running Reactome | positive-selection ===")
res_reactome_pos <- run_and_save_reactome(
  ranks_entrez = ranks_pos_entrez,
  prefix = paste0(run_label, "__Reactome__positive_selection"),
  out_dir_local = reactome_dir
)

message("=== Running Reactome | negative-selection ===")
res_reactome_neg <- run_and_save_reactome(
  ranks_entrez = ranks_neg_entrez,
  prefix = paste0(run_label, "__Reactome__negative_selection"),
  out_dir_local = reactome_dir
)

reactome_pos_tbl <- safe_result_tbl(res_reactome_pos)
reactome_neg_tbl <- safe_result_tbl(res_reactome_neg)

if (nrow(reactome_pos_tbl) > 0) summary_rows[[length(summary_rows) + 1]] <- reactome_pos_tbl %>% mutate(collection = "Reactome", arm = "positive_selection")
if (nrow(reactome_neg_tbl) > 0) summary_rows[[length(summary_rows) + 1]] <- reactome_neg_tbl %>% mutate(collection = "Reactome", arm = "negative_selection")

# ------------------------------------------------------------
# top gene tables (by rank_metric, independent of GSEA)
# ------------------------------------------------------------
top_pos <- df_pos_dedup %>%
  mutate(rank_metric = compute_lfc_p_metric(pos_lfc, pos_p, flip_sign = FALSE)) %>%
  arrange(desc(rank_metric), pos_fdr, pos_p) %>%
  slice_head(n = 100)

top_neg <- df_neg_dedup %>%
  mutate(rank_metric = compute_lfc_p_metric(neg_lfc, neg_p, flip_sign = TRUE)) %>%
  arrange(desc(rank_metric), neg_fdr, neg_p) %>%
  slice_head(n = 100)

safe_write_csv(top_pos, file.path(out_dir, "TOP_positive_selection_genes.csv"))
safe_write_csv(top_neg, file.path(out_dir, "TOP_negative_selection_genes.csv"))

# ------------------------------------------------------------
# combined summary tables (input for 02b / supplementary script)
# ------------------------------------------------------------
summary_tbl <- bind_rows(summary_rows)

if (nrow(summary_tbl) > 0) {
  safe_write_csv(summary_tbl, file.path(out_dir, "GSEA_summary_all_results.csv"))
}

keep_cols <- intersect(
  c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust",
    "qvalues", "rank", "leading_edge", "core_enrichment", "collection", "arm"),
  names(summary_tbl)
)

summary_tbl_small <- summary_tbl %>% dplyr::select(all_of(keep_cols))
safe_write_csv(summary_tbl_small, file.path(out_dir, "GSEA_summary_compact.csv"))

message("Done. CRISPR GSEA compute output written to: ", out_dir)
message("Run label (used by 02b and the supplementary script): ", run_label)

