# ======================================================================
# 06_screening_candidates_in_vitro_step1.R
#
# Table S8: summary of the 15 candidate genes taken into the in vitro
# screen. For each candidate:
#   * spatial association with invaded nerves (Lee's L vs 200 relocated
#     masks, 02 step1): delta L, combined test, patients above the
#     relocated 95th percentile, core-signature membership; interface test.
#     Core signature = upregulated in InvN (DE) AND above the relocated
#     controls in all 3 patients; an LR-derived gene that is not upregulated
#     can be nerve-associated in 3/3 patients without being "core";
#   * differential expression InvN vs NI (TableS3);
#   * ligand-receptor support in the revised LR analysis (TableS10:
#     Tier 1 = filtered 130 pairs, Tier 2);
#   * mean expression per nerve category and patient;
#   * provenance: how the candidate was nominated.
# The descriptive "Summary" column is generated from these values (no
# hand-written rationale), so it always matches the numbers.
#
# GENE IDENTITY (placental growth factor): the candidate was nominated from
# an earlier LR list as "PIGF^FLT1". In the LR resource that pair is listed
# under the symbol PIGF, but the official symbol of placental growth factor
# is PGF (NCBI 5228, ENSG00000119630); PIGF (NCBI 5281) is an unrelated
# GPI-anchor biosynthesis enzyme. LR scores computed for "PIGF^FLT1"
# therefore used the expression of the PIGF (GPI) gene. The in vitro screen
# used recombinant human PlGF-1 (the PGF product). This table reports the
# candidate as PGF and records the nomination symbol in `Nominated_as`.
#
# INPUTS: input/Integrated_NMF.rds, input/DE_SPE_annotaion.csv,
#   input/TableS3_Full_DE_invaded_vs_non_invaded_nerve.csv,
#   input/TableS10_all_LR_cosine_score.xls,
#   output/02_spatial_association/tables/lee_sham_all_candidates.csv (02 step1)
# OUTPUTS (output/06_candidate_gene_screening/):
#   tables/Table_S8_candidate_summary.csv, screening_candidates_summary.csv,
#   screening_candidates_expression_by_patient.csv,
#   screening_candidates_expression_pooled.csv, figures/Table_S8_candidate_summary.pdf
# ======================================================================

library(dplyr)
library(tidyr)
library(readxl)
library(SpatialExperiment)
library(gridExtra)
library(grid)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
summarise <- dplyr::summarise; group_by <- dplyr::group_by; arrange <- dplyr::arrange

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input"))
  stop("Run this script from the project folder (the one containing input/).")

PROJECT_DIR   <- "."
INPUT_DIR     <- file.path(PROJECT_DIR, "input")
ANALYSIS_NAME <- "06_candidate_gene_screening"
TABLES_DIR    <- file.path(PROJECT_DIR, "output", ANALYSIS_NAME, "tables")
FIGURES_DIR   <- file.path(PROJECT_DIR, "output", ANALYSIS_NAME, "figures")
# Written by 02_spatial_association_with_nerve_step1_lee_sham.R
# (spatial_enrichment_prioritisation); place a copy in input/ before running.
LEE_TABLE     <- file.path(INPUT_DIR, "lee_sham_all_candidates.csv")
LR_FILE       <- file.path(INPUT_DIR, "TableS10_all_LR_cosine_score.xls")
LR_TIER1_SHEET <- "filtered_130_LR_pairs (Tier1)"
LR_TIER2_SHEET <- "Tier2"
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")

# The 15 candidates, defined once. `Gene` = official HGNC symbol used for all
# look-ups; `Nominated_as` = symbol under which the candidate was nominated.
candidates <- data.frame(
  Gene         = c("ADAM23", "BTC", "CLCF1", "DRAXIN", "FGF5", "FLT3LG", "IL11", "IL34",
                   "NLGN4X", "NRCAM", "PDGFB", "PGF", "SFRP1", "SLITRK2", "WNT5A"),
  Nominated_as = c("ADAM23", "BTC", "CLCF1", "DRAXIN", "FGF5", "FLT3LG", "IL11", "IL34",
                   "NLGN4X", "NRCAM", "PDGFB", "PIGF", "SFRP1", "SLITRK2", "WNT5A"),
  stringsAsFactors = FALSE)
provenance_note <- c(
  PGF   = paste("Nominated as 'PIGF^FLT1' in an earlier LR analysis; the LR resource lists placental growth",
                "factor under the symbol PIGF, which in the expression data is the GPI-anchor gene PIGF.",
                "Screened as recombinant PlGF-1 (PGF)."),
  WNT5A = "Nominated from an earlier LR analysis (WNT5A^FZD4 / WNT5A^FZD8).")

# ---- 1. Spot data (annotation attached by barcode) --------------------------
if (!exists("spe") || !inherits(spe, "SpatialExperiment")) spe <- readRDS(file.path(INPUT_DIR, "Integrated_NMF.rds"))
attach_annotation <- function(spe, path) {
  ann <- read.csv(path)
  if (nrow(ann) != ncol(spe)) stop("Annotation has ", nrow(ann), " rows but the object has ", ncol(spe), " spots.")
  if (anyDuplicated(ann$barcode)) stop("Duplicated barcodes in the annotation file.")
  keys <- list(colData_barcode = if ("barcode" %in% names(colData(spe))) as.character(spe$barcode),
               colnames = colnames(spe))
  keys <- Filter(function(k) !is.null(k) && !anyNA(match(k, ann$barcode)), keys)
  if (!length(keys)) stop("Annotation barcodes do not match the object's barcodes or column names.")
  ann <- ann[match(keys[[1]], ann$barcode), , drop = FALSE]
  if ("sample_id" %in% names(colData(spe)) && !all(as.character(spe$sample_id) == as.character(ann$sample_id)))
    stop("Annotation sample_id does not match the object's sample_id for the same barcodes.")
  colData(spe) <- DataFrame(ann); colnames(spe) <- ann$barcode
  spe
}
spe <- attach_annotation(spe, file.path(INPUT_DIR, "DE_SPE_annotaion.csv"))
missing_genes <- setdiff(candidates$Gene, rownames(spe))
if (length(missing_genes)) warning("Not in the expression data: ", paste(missing_genes, collapse = ", "))

# ---- 2. Spatial association (02 step1) --------------------------------------
lee <- read.csv(LEE_TABLE)
spatial <- lee %>% filter(Gene %in% candidates$Gene) %>%
  transmute(Gene, Tested_pool = ifelse(family == "DE_pass", "Upregulated (DE) genes", "LR-derived genes"),
            Nerve_deltaL = nerve_delta_L_mean, Nerve_padj = nerve_padj,
            Nerve_patients_above_q95 = nerve_n_patients_exceed_q95,
            Nerve_association = nerve_status == "pass",
            Core_signature = if ("signature_core" %in% names(lee)) signature_core %in% TRUE else
              family == "DE_pass" & nerve_status == "pass" & nerve_n_patients_exceed_q95 == 3,
            Interface_deltaL = interface_delta_L_mean, Interface_association = interface_status == "pass")
not_tested <- setdiff(candidates$Gene, spatial$Gene)
if (length(not_tested)) cat("Not in the spatial test (02 step1):", paste(not_tested, collapse = ", "), "\n")

# ---- 3. Differential expression (TableS3) ------------------------------------
full_de <- read.csv(file.path(INPUT_DIR, "TableS3_Full_DE_invaded_vs_non_invaded_nerve.csv"), check.names = FALSE)
names(full_de) <- trimws(names(full_de))
full_de <- full_de %>% transmute(Gene = trimws(Gene), DE_log2FC = `log2FC (logarithm of fold change)`,
                                 DE_padj = `adj.P.Val (FDR)`) %>% distinct(Gene, .keep_all = TRUE)

# ---- 4. Ligand-receptor support (revised LR analysis, TableS10) -------------
lr_genes <- function(pair) unlist(strsplit(unlist(strsplit(pair, "\\^")), "_"))
lr_hits <- function(sheet) {
  lr <- read_excel(LR_FILE, sheet = sheet)
  bind_rows(lapply(candidates$Gene, function(g) {
    hit <- lr$Interaction[vapply(lr$Interaction, function(p) g %in% lr_genes(p), logical(1))]
    data.frame(Gene = g, pairs = paste(hit, collapse = "; "))
  }))
}
lr <- full_join(lr_hits(LR_TIER1_SHEET) %>% rename(LR_Tier1_pairs = pairs),
                lr_hits(LR_TIER2_SHEET) %>% rename(LR_Tier2_pairs = pairs), by = "Gene") %>%
  mutate(LR_support = case_when(LR_Tier1_pairs != "" ~ "Tier 1", LR_Tier2_pairs != "" ~ "Tier 2",
                                TRUE ~ "Not in Tier 1/2 pairs"))

# ---- 5. Expression per nerve category and patient ---------------------------
spe$nerve_category <- case_when(
  !is.na(spe$NerveID) & spe$context == "Adjacent" & spe$Invaded == "NonInvaded" ~ "NI_Adjacent",
  !is.na(spe$NerveID) & spe$context == "Diseased" & spe$Invaded == "NonInvaded" ~ "NI_TumourBearing",
  !is.na(spe$NerveID) & spe$context == "Diseased" & spe$Invaded == "Invaded"    ~ "InvN",
  is.na(spe$NerveID)  & !is.na(spe$TumorDistalID) & spe$TumorDistalID != "none" ~ "Tumour_distant",
  TRUE ~ NA_character_)
cat("\nSpots per category:\n"); print(table(spe$nerve_category, useNA = "always"))
genes_expr <- intersect(candidates$Gene, rownames(spe))
ex <- as.matrix(assay(spe, "logcounts")[genes_expr, , drop = FALSE])
expression_detail <- bind_rows(lapply(genes_expr, function(g)
  data.frame(expr = ex[g, ], Category = spe$nerve_category, Patient = spe$patient_id) %>%
    filter(!is.na(Category)) %>% group_by(Category, Patient) %>%
    summarise(mean_expr = mean(expr), pct_detected = 100 * mean(expr > 0), n_spots = n(), .groups = "drop") %>%
    mutate(Gene = g)))
write.csv(expression_detail, file.path(TABLES_DIR, "screening_candidates_expression_by_patient.csv"), row.names = FALSE)
expression_pooled <- expression_detail %>% group_by(Gene, Category) %>%
  summarise(mean_expr = mean(mean_expr), .groups = "drop") %>%                       # mean of patient means
  pivot_wider(names_from = Category, values_from = mean_expr, names_prefix = "meanExpr_")
write.csv(expression_pooled, file.path(TABLES_DIR, "screening_candidates_expression_pooled.csv"), row.names = FALSE)

# ---- 6. Summary table ---------------------------------------------------------
summary_table <- candidates %>%
  left_join(spatial, by = "Gene") %>% left_join(full_de, by = "Gene") %>% left_join(lr, by = "Gene") %>%
  left_join(expression_pooled, by = "Gene") %>%
  mutate(
    # Same wording for every gene; "core signature" is added only for genes in
    # it (core = upregulated in InvN AND above the relocated controls in 3/3
    # patients, so LR-derived genes that are not upregulated can be 3/3 but
    # are never core).
    Spatial = case_when(
      is.na(Nerve_deltaL)  ~ "Not in the spatial test",
      Nerve_association    ~ paste0("Nerve-associated (", Nerve_patients_above_q95, "/3 patients above controls",
                                    ifelse(Core_signature, "; core signature", ""), ")"),
      TRUE                 ~ "Not nerve-associated"),
    DE = case_when(is.na(DE_padj) ~ "Not tested in DE",
                   DE_padj < 0.1 & DE_log2FC > 0 ~ "Up in InvN vs NI",
                   DE_padj < 0.1 & DE_log2FC < 0 ~ "Down in InvN vs NI",
                   TRUE ~ "No significant change"),
    Summary = paste0(Spatial, "; ", DE, "; LR: ", LR_support),
    Provenance = ifelse(Gene %in% names(provenance_note), provenance_note[Gene], "")) %>%
  arrange(desc(Nerve_deltaL))
if (any(summary_table$Gene == "PGF" & is.na(summary_table$Nerve_deltaL)))
  summary_table$Summary[summary_table$Gene == "PGF"] <-
  paste0("Not in the spatial test (nominated as PIGF); ", summary_table$DE[summary_table$Gene == "PGF"],
         "; LR: ", summary_table$LR_support[summary_table$Gene == "PGF"])
write.csv(summary_table, file.path(TABLES_DIR, "screening_candidates_summary.csv"), row.names = FALSE)

final_table <- summary_table %>%
  select(Gene, Nominated_as, Nerve_deltaL, Nerve_patients_above_q95, Core_signature, Interface_deltaL,
         DE_log2FC, DE_padj, any_of(c("meanExpr_NI_Adjacent", "meanExpr_NI_TumourBearing", "meanExpr_InvN",
                                      "meanExpr_Tumour_distant")),
         LR_support, LR_Tier1_pairs, LR_Tier2_pairs, Summary, Provenance)
write.csv(final_table, file.path(TABLES_DIR, "Table_S8_candidate_summary.csv"), row.names = FALSE)
cat("\nTable S8:\n"); print(as.data.frame(final_table %>% select(Gene, Nerve_deltaL, Nerve_patients_above_q95,
                                                                 DE_log2FC, DE_padj, LR_support, Summary)), digits = 3)

# ---- 7. Rendered table (generated, no hand-curated copy) ---------------------
fmt <- function(x, d = 2) ifelse(is.na(x), "-", formatC(x, format = "f", digits = d))
display <- final_table %>% transmute(
  Gene = ifelse(Gene != Nominated_as, paste0(Gene, " (as ", Nominated_as, ")"), Gene),
  `Nerve dL` = fmt(Nerve_deltaL, 3),
  `Patients > controls` = ifelse(is.na(Nerve_patients_above_q95), "-", paste0(Nerve_patients_above_q95, "/3")),
  `Interface dL` = fmt(Interface_deltaL, 3),
  `log2FC` = fmt(DE_log2FC), `FDR` = ifelse(is.na(DE_padj), "-", formatC(DE_padj, format = "g", digits = 2)),
  `NI-Adj` = fmt(meanExpr_NI_Adjacent), `NI-TB` = fmt(meanExpr_NI_TumourBearing), `InvN` = fmt(meanExpr_InvN),
  `LR` = LR_support, Summary = Summary)
tbl <- tableGrob(display, rows = NULL,
                 theme = ttheme_minimal(base_size = 8, base_family = "Helvetica",
                                        core = list(bg_params = list(col = "grey85", lwd = 0.5)),
                                        colhead = list(bg_params = list(col = "grey85", lwd = 0.5))))
w <- convertWidth(sum(tbl$widths), "in", valueOnly = TRUE) + 0.3
h <- convertHeight(sum(tbl$heights), "in", valueOnly = TRUE) + 0.3
pdf_device(file.path(FIGURES_DIR, "Table_S8_candidate_summary.pdf"), width = w, height = h)
grid.draw(tbl); invisible(dev.off())
cat("\nSaved Table S8 (csv + pdf) in", file.path(PROJECT_DIR, "output", ANALYSIS_NAME), "\n")
