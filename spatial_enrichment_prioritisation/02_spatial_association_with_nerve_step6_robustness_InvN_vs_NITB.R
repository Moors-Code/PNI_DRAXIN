# ======================================================================
# 02_spatial_association_with_nerve_step6_robustness_InvN_vs_NITB.R
#
# ROBUSTNESS CHECK, not used in the manuscript. The manuscript reports the
# differential expression of invaded nerves against ALL non-invaded nerves
# (InvN vs NI, 04_DE_analysis). This script repeats the comparison against
# non-invaded nerves OF THE SAME TISSUE CONTEXT only (InvN vs NI-TB): both
# groups then lie in tumour-bearing sections, so the comparison is within
# section and cannot be carried by the differences between the tumour-bearing
# and the adjacent block. It is a quality check on the reported contrast; its
# outputs are named "robustness_" and no figure or number from it is in the
# paper.
#
# WHAT IT DOES
#   The pseudobulk pipeline of 04_DE_analysis is reproduced exactly, with one
#   change: the area factor is refined by tissue context, so that non-invaded
#   nerves of tumour-bearing sections (NI-TB) and of adjacent sections (NI-Adj)
#   become separate levels. The model is refitted on ALL pseudobulks - keeping
#   every group improves the variance estimate - with section as a blocking
#   factor via duplicateCorrelation, and one contrast is evaluated:
#
#       InvadedNerve_Diseased  -  Nerve_Diseased          (InvN vs NI-TB)
#
#   The test is genome-wide. The spatially prioritised genes are extracted from
#   the finished result afterwards, so no gene set is redefined and the spatial
#   prioritisation is NOT repeated.
#
#   Nothing in the DE project is modified; the input is the saved bundle
#   DE_analysis_results.rds (its element `sub`, the unaggregated subset) and all
#   outputs go to the 02_spatial_association folder under new file names.
#
# OUTPUTS
#   tables/robustness_InvN_vs_NITB_all_genes.csv          every tested gene
#   tables/robustness_InvN_vs_NITB_prioritised_genes.csv  the prioritised genes only
#   tables/robustness_InvN_vs_NITB_summary.txt            counts
#   figures/robustness_InvN_vs_NITB_heatmap.pdf           per-nerve expression of the
#       prioritised genes that remain significant, z-scored per gene
# ======================================================================

library(SpatialExperiment)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(scuttle)
library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(ggplot2)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; rename <- dplyr::rename

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input"))
  stop("Run this script from the project folder (the one containing input/).")

PROJECT_DIR   <- "."

# Written by 04a_DE_analysis_COMPUTE.Rmd (spatial_data_processing); place a
# copy in input/ before running.
DE_BUNDLE   <- file.path(PROJECT_DIR, "input", "DE_analysis_results.rds")
TABLES_DIR  <- file.path(PROJECT_DIR, "output", "02_spatial_association", "tables")
FIGURES_DIR <- file.path(PROJECT_DIR, "output", "02_spatial_association", "figures")
GENE_SET_FILE <- file.path(TABLES_DIR, "lee_sham_signature_core.csv")
for (d in c(TABLES_DIR, FIGURES_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

FDR_CUT   <- 0.1    # as in the original differential expression analysis
MIN_CELLS <- 1      # pseudobulks kept if ncells > MIN_CELLS and < MAX_CELLS
MAX_CELLS <- 100

GRP_INV <- "InvadedNerve_Diseased"   # InvN  (tumour-bearing sections)
GRP_NI  <- "Nerve_Diseased"          # NI-TB (tumour-bearing sections)

pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")

# Section -> patient / context, identical to 04_DE_analysis
sample_key <- data.frame(
  sample_id = c("54703bl38_noPNI", "54703bl28_highPNI",
                "6249bl23_noPNI",  "6249bl18_noPNIandhighPNI",
                "70865bl5_noPNI",  "70865bl15_highPNI"),
  patient   = c("Patient 1", "Patient 1", "Patient 2", "Patient 2", "Patient 3", "Patient 3"),
  context   = c("Adjacent", "Diseased", "Adjacent", "Diseased", "Adjacent", "Diseased"),
  stringsAsFactors = FALSE)

# ---- 1. Input: the unaggregated subset from the DE analysis ----------------
if (!file.exists(DE_BUNDLE))
  stop("DE bundle not found: ", DE_BUNDLE,
       "\nPlace DE_analysis_results.rds (04a_DE_analysis_COMPUTE.Rmd) in input/.")
bundle <- readRDS(DE_BUNDLE)
if (!"sub" %in% names(bundle))
  stop("The DE bundle contains no element 'sub' (the unaggregated subset).")
sub <- bundle$sub
rm(bundle); invisible(gc())
cat("Spots loaded:", ncol(sub), "in", nrow(sub), "genes\n")

sub$context <- as.character(sample_key$context[match(sub$sample_id, sample_key$sample_id)])
sub$patient <- as.character(sample_key$patient[match(sub$sample_id, sample_key$sample_id)])
if (any(is.na(sub$context)))
  stop("Sections absent from sample_key: ",
       paste(unique(sub$sample_id[is.na(sub$context)]), collapse = ", "))

# Area refined by tissue context. This is the only departure from
# 04_DE_analysis: "Nerve" is split into Nerve_Diseased (NI-TB) and
# Nerve_Adjacent (NI-Adj).
sub$grp <- paste(sub$Area, sub$context, sep = "_")
cat("\nSpots per group and section:\n")
print(table(sub$grp, sub$sample_id))

# ---- 2. Pseudobulk, exactly as in 04_DE_analysis ---------------------------
sub$ID <- sub$NerveID
sub$ID[is.na(sub$ID)] <- sub$BorderID[is.na(sub$ID)]
sub$ID[is.na(sub$ID)] <- sub$TumorDistalID[is.na(sub$ID)] + 1000

ids <- DataFrame(grp = sub$grp, sample_id = sub$sample_id, id = sub$ID)
sb <- aggregateAcrossCells(sub, ids)

levs <- paste(colData(sb)$grp, colData(sb)$id, colData(sb)$sample_id, sep = ".")
colData(sb)$level_id <- levs
colnames(sb) <- levs

sb <- sb[, sb$ncells > MIN_CELLS & sb$ncells < MAX_CELLS]
sb <- computeLibraryFactors(sb)
sb <- logNormCounts(sb)

cd <- as.data.frame(colData(sb))
cat("\nPseudobulks per group:\n"); print(table(cd$grp))

if (!all(c(GRP_INV, GRP_NI) %in% cd$grp))
  stop("Contrast groups missing. Groups present: ",
       paste(sort(unique(cd$grp)), collapse = ", "))

# Groups with a single pseudobulk cannot be estimated; drop them from the model.
keep_grp <- names(which(table(cd$grp) >= 2))
if (!all(c(GRP_INV, GRP_NI) %in% keep_grp))
  stop("A contrast group has fewer than two pseudobulks.")
sb <- sb[, cd$grp %in% keep_grp]
cd <- as.data.frame(colData(sb))
cd$grp <- factor(cd$grp)
cat("\nGroups entering the model:\n"); print(table(cd$grp, cd$sample_id))
cat("\nContrast: ", GRP_INV, " (n = ", sum(cd$grp == GRP_INV), ") vs ",
    GRP_NI, " (n = ", sum(cd$grp == GRP_NI), ")\n", sep = "")

# ---- 3. Fit and contrast ---------------------------------------------------
mdl <- model.matrix(~ 0 + grp, data = cd)
colnames(mdl) <- sub("^grp", "", colnames(mdl))

keep <- edgeR::filterByExpr(logcounts(sb), design = mdl, min.count = 2)
cat("\nGenes retained after expression filtering:", sum(keep), "of", nrow(sb), "\n")
sb <- sb[keep, ]
m <- logcounts(sb)

corfit <- duplicateCorrelation(m, mdl, block = cd$sample_id)
cat("Consensus correlation:", round(corfit$consensus.correlation, 4), "\n")
fit <- lmFit(m, design = mdl, block = cd$sample_id,
             correlation = corfit$consensus.correlation)

ctr  <- makeContrasts(paste0(GRP_INV, "-", GRP_NI), levels = colnames(mdl))
fit2 <- eBayes(contrasts.fit(fit, contrasts = ctr))

res <- topTable(fit2, n = nrow(fit2))
res$Gene <- rownames(res)
res <- res %>%
  transmute(Gene,
            logFC_InvN_vs_NITB = logFC,
            AveExpr, t, P.Value, FDR = adj.P.Val, B) %>%
  arrange(FDR)
write.csv(res, file.path(TABLES_DIR, "robustness_InvN_vs_NITB_all_genes.csv"), row.names = FALSE)
cat("\nGenome-wide:", nrow(res), "genes tested;", sum(res$FDR < FDR_CUT),
    "at FDR <", FDR_CUT, "of which",
    sum(res$FDR < FDR_CUT & res$logFC_InvN_vs_NITB > 0), "higher in InvN\n")

# ---- 4. The spatially prioritised genes in this contrast -------------------
core <- unique(trimws(read.csv(GENE_SET_FILE)$Gene))
prio <- res %>%
  filter(Gene %in% core) %>%
  mutate(significant_InvN_vs_NITB = FDR < FDR_CUT & logFC_InvN_vs_NITB > 0) %>%
  arrange(FDR)
missing <- setdiff(core, res$Gene)
write.csv(prio, file.path(TABLES_DIR, "robustness_InvN_vs_NITB_prioritised_genes.csv"), row.names = FALSE)

n_sig <- sum(prio$significant_InvN_vs_NITB)
summary_lines <- c(
  paste0("InvN vs NI-TB, within tumour-bearing sections - ", Sys.Date()),
  "",
  paste0("Pseudobulks: ", sum(cd$grp == GRP_INV), " invaded nerves vs ",
         sum(cd$grp == GRP_NI), " non-invaded nerves of tumour-bearing sections"),
  paste0("Model: ~ 0 + area x context over all ", length(keep_grp),
         " estimable groups; section as blocking factor (consensus correlation ",
         round(corfit$consensus.correlation, 3), ")"),
  paste0("Genes tested: ", nrow(res), "; significant at FDR < ", FDR_CUT, ": ",
         sum(res$FDR < FDR_CUT)),
  "",
  paste0("Spatially prioritised genes: ", length(core), " (", nrow(prio),
         " tested, ", length(missing), " filtered on expression)"),
  paste0("  higher in InvN at FDR < ", FDR_CUT, ": ", n_sig,
         " (", round(100 * n_sig / max(nrow(prio), 1)), "% of those tested)"),
  paste0("  positive log2FC irrespective of significance: ",
         sum(prio$logFC_InvN_vs_NITB > 0), " of ", nrow(prio)),
  "",
  "Summary:",
  paste0("  \"", n_sig, " of the ", length(core), " spatially prioritised genes remained ",
         "significantly higher in InvN than in non-invaded nerves of the same ",
         "tumour-bearing sections (FDR < ", FDR_CUT, "; Table S--).\""))
if (length(missing))
  summary_lines <- c(summary_lines, "",
                     paste0("Not tested: ", paste(missing, collapse = ", ")))
writeLines(summary_lines, file.path(TABLES_DIR, "robustness_InvN_vs_NITB_summary.txt"))
cat("\n"); cat(summary_lines, sep = "\n"); cat("\n")

# ---- 5. Heatmap of the significant prioritised genes -----------------------
# Columns are the nerve pseudobulks that entered the contrast, grouped by
# category and ordered by patient; values are per-gene z-scores of logcounts,
# so rows are comparable. The plot shows the data the test was run on, in the
# units the test used.
sig_genes <- prio$Gene[prio$significant_InvN_vs_NITB]
if (length(sig_genes) < 2) {
  cat("\nFewer than two significant prioritised genes: heatmap skipped.\n")
} else {
  cols <- cd$grp %in% c(GRP_INV, GRP_NI)
  mat  <- as.matrix(logcounts(sb)[sig_genes, cols, drop = FALSE])
  ann  <- cd[cols, c("grp", "patient", "sample_id"), drop = FALSE]
  ann$pseudobulk <- colnames(sb)[cols]
  ann$category <- factor(ifelse(ann$grp == GRP_INV, "InvN", "NI-TB"),
                         levels = c("NI-TB", "InvN"))
  colnames(mat) <- ann$pseudobulk
  
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  
  col_order <- ann$pseudobulk[order(ann$category, ann$patient, ann$pseudobulk)]
  row_order <- sig_genes[order(prio$logFC_InvN_vs_NITB[match(sig_genes, prio$Gene)])]
  
  hm <- as.data.frame(z) %>%
    mutate(Gene = rownames(z)) %>%
    pivot_longer(-Gene, names_to = "pseudobulk", values_to = "z") %>%
    left_join(ann, by = "pseudobulk") %>%
    mutate(Gene       = factor(Gene, levels = row_order),
           pseudobulk = factor(pseudobulk, levels = col_order))
  
  lim <- as.numeric(quantile(abs(hm$z), 0.99, na.rm = TRUE))
  p <- ggplot(hm, aes(pseudobulk, Gene, fill = pmax(pmin(z, lim), -lim))) +
    geom_tile() +
    facet_grid(. ~ category, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-lim, lim),
                         name = "Expression (z-score)") +
    labs(x = "Nerve pseudobulks", y = NULL) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6, face = "italic", colour = "black"),
          strip.text = element_text(size = 7, face = "bold"),
          legend.position = "bottom",
          legend.key.height = unit(2.5, "mm"), legend.key.width = unit(8, "mm"),
          plot.margin = margin(2, 4, 2, 2, "mm"))
  
  h_mm <- 26 + 2.6 * length(sig_genes)
  ggsave(file.path(FIGURES_DIR, "robustness_InvN_vs_NITB_heatmap.pdf"), p,
         width = 90, height = h_mm, units = "mm", device = pdf_device, limitsize = FALSE)
  cat("\nSaved robustness_InvN_vs_NITB_heatmap.pdf:", length(sig_genes), "genes,",
      sum(cols), "pseudobulks\n")
}

cat("\nDone.\n")
