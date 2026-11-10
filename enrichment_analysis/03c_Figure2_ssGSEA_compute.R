# ============================================================
# 03c_Figure2_ssGSEA_compute.R
#
# Pathway activity profiling across nerve categories (spatial data).
#
# COMPUTE ONLY: subsets nerve spots from a SpatialExperiment object,
# scores each spot for a fixed panel of MSigDB pathways via ssGSEA,
# aggregates spot-level scores to the nerve level, and standardizes
# (z-score) across nerves. Saves one bundle only. No plotting here --
# see 03d_Figure2_ssGSEA_plots.R.
#
# NOTE: the original script had an interactive gene-set discovery step
# (searching already-computed GSEA CSVs by keyword to find exact MSigDB
# gs_name spellings) that is NOT reproduced here -- it was a one-time
# lookup aid, not needed to rerun the analysis. The gene sets below are
# the final, already-identified panel from that lookup.
#
# Input : <external_spatial_dir>/Integrated_NMF.rds        (SpatialExperiment)
#         <external_spatial_dir>/DE_SPE_annotaion.csv      (colData)
# Output: processed/Figure2_GSEA/ssGSEA/nerve_scores_data.rds
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(GSVA)
  library(matrixStats)
  library(msigdbr)
  library(Matrix)
  library(SpatialExperiment)
  library(readr)
  library(tibble)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
# This spatial data lives in a separate project, not copied into this
# project's input/ folder (it's large + shared/collaborator data).
# Update this if the external project moves.
external_spatial_dir <- file.path("../Spatial_data/paper_revision_codes/input/")

spe_path     <- file.path(external_spatial_dir, "Integrated_NMF.rds")
coldata_path <- file.path(external_spatial_dir, "DE_SPE_annotaion.csv")

if (!file.exists(spe_path))     stop("SpatialExperiment file not found: ", spe_path)
if (!file.exists(coldata_path)) stop("colData csv not found: ", coldata_path)

processed_dir <- file.path("processed", "Figure2_GSEA", "ssGSEA")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# load SpatialExperiment + metadata
# ------------------------------------------------------------
spe <- readRDS(spe_path)
colData(spe) <- DataFrame(read.csv(coldata_path))
colnames(spe) <- spe$barcode

# ------------------------------------------------------------
# subset nerve spots and define the three nerve categories
# ------------------------------------------------------------
nerve_spe <- spe[, colData(spe)$is.nerve == TRUE]

nerve_spe$types <- nerve_spe$type
colData(nerve_spe)$types[!colData(nerve_spe)$types == "noPNI"] <- "Diseased"
colData(nerve_spe)$types[colData(nerve_spe)$types == "noPNI"] <- "Adjacent"

colData(nerve_spe)$nerve_types <- nerve_spe$Invaded
colData(nerve_spe)$nerve_types <- ifelse(
  colData(nerve_spe)$Invaded == "NonInvaded",
  paste0(colData(nerve_spe)$Invaded, "_", colData(nerve_spe)$types),
  as.character(colData(nerve_spe)$nerve_types)
)

nerve_spe <- nerve_spe[, !is.na(colData(nerve_spe)$nerve_types)]

nerve_type_order <- c("NonInvaded_Adjacent", "NonInvaded_Diseased", "Invaded")

# ------------------------------------------------------------
# prepare expression and metadata for scoring
# ------------------------------------------------------------
meta <- as.data.frame(colData(nerve_spe))
meta$nerve_types <- factor(as.character(meta$nerve_types), levels = nerve_type_order)

# patient order -- edit here if the patient cohort changes
patient_order <- c("54703", "6249", "70865")
meta$patient_id <- factor(as.character(meta$patient_id), levels = patient_order)

assay_name <- "logcounts"
expr <- assay(nerve_spe, assay_name)
expr <- as.matrix(expr)
stopifnot(!is.null(rownames(expr)))

# ------------------------------------------------------------
# MSigDB gene-set lookup helper (exact match only -- the panel below
# was already identified via the original interactive search)
# ------------------------------------------------------------
msig_exact <- function(gs_name_in, species = "Homo sapiens") {
  tb <- msigdbr(species = species) %>% dplyr::filter(.data$gs_name == gs_name_in)
  if (nrow(tb) == 0) stop("msigdbr: exact gs_name not found: ", gs_name_in)
  unique(tb$gene_symbol)
}

# ------------------------------------------------------------
# fixed pathway panel
# ------------------------------------------------------------
gene_sets <- list(
  "H: INFLAMMATORY_RESPONSE" = msig_exact("HALLMARK_INFLAMMATORY_RESPONSE"),
  "H: IL6_JAK_STAT3"         = msig_exact("HALLMARK_IL6_JAK_STAT3_SIGNALING"),
  "H: IFNG_RESPONSE"         = msig_exact("HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  "H: IFNA_RESPONSE"         = msig_exact("HALLMARK_INTERFERON_ALPHA_RESPONSE"),

  "C4: GAVISH_FIBROBLASTS_COMPLEMENT" =
    msig_exact("GAVISH_3CA_METAPROGRAM_FIBROBLASTS_COMPLEMENT"),

  "C2: WP_BCELL_SARS_COV2" =
    msig_exact("WP_EXTRAFOLLICULAR_AND_FOLLICULAR_B_CELL_ACTIVATION_BY_SARSCOV2"),

  "C5: GOBP_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION" =
    msig_exact("GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"),

  "C5: GOBP_RESPONSE_TO_PROSTAGLANDIN_E" =
    msig_exact("GOBP_RESPONSE_TO_PROSTAGLANDIN_E"),

  "C2: REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT" =
    msig_exact("REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT"),

  "C5: GOBP_LEUKOCYTE_TETHERING_OR_ROLLING" =
    msig_exact("GOBP_LEUKOCYTE_TETHERING_OR_ROLLING")
)

# Intersect with detected genes, keep sets with >= 10 genes present
gene_sets_use <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr)))
gene_sets_use <- gene_sets_use[sapply(gene_sets_use, length) >= 10]

message("Gene set sizes after intersection:")
print(sort(sapply(gene_sets_use, length)))

stopifnot(length(gene_sets_use) >= 2)

# ------------------------------------------------------------
# compute ssGSEA pathway scores per spot
# ------------------------------------------------------------
ssgsea_scores <- gsva(
  expr, gene_sets_use,
  method = "ssgsea", kcdf = "Gaussian", abs.ranking = FALSE
)
# result: pathways x spots

# ------------------------------------------------------------
# aggregate spot-level scores to the nerve level
# ------------------------------------------------------------
score_df <- as.data.frame(t(ssgsea_scores))
score_df$spot_id <- rownames(score_df)

score_long <- score_df %>%
  tidyr::pivot_longer(cols = -spot_id, names_to = "pathway", values_to = "ssgsea_score") %>%
  left_join(meta %>% tibble::rownames_to_column("spot_id"), by = "spot_id")

nerve_scores <- score_long %>%
  group_by(NerveID, pathway, nerve_types, patient_id) %>%
  summarise(mean_score = mean(ssgsea_score), .groups = "drop")

# ------------------------------------------------------------
# standardize (z-score) pathway scores across nerves
# ------------------------------------------------------------
nerve_scores <- nerve_scores %>%
  group_by(pathway) %>%
  mutate(z_score = as.numeric(scale(mean_score))) %>%
  ungroup()

# ------------------------------------------------------------
# save bundle for 03d (and the stats supplementary script)
# ------------------------------------------------------------
saveRDS(
  list(
    nerve_scores = nerve_scores,
    gene_sets_use = gene_sets_use,
    meta = meta,
    nerve_type_order = nerve_type_order,
    patient_order = patient_order
  ),
  file.path(processed_dir, "nerve_scores_data.rds")
)

message("Done. Bundle written to: ", file.path(processed_dir, "nerve_scores_data.rds"))
message("Run 03d_Figure2_ssGSEA_plots.R next for the paper figures.")
