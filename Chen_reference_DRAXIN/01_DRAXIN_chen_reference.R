# ============================================================
# 03_Draxin_chen_reference.R
# DRAXIN expression across cell types -- Chen et al. (2025)
# sc/snRNA-seq reference cohort (GSE278688 / GSE278689)
# ============================================================
#
# PURPOSE
#   Independent, cell-type-level corroboration that DRAXIN transcript is
#   preferentially found in Schwann cells. Complements, but does not
#   replace, the spatial nerve-level results (Fig 6A-B) and the in situ
#   HCR/IF imaging (Fig S10B-C). Dissociated sc/snRNA-seq has no spatial
#   context, so this analysis says nothing about invaded vs non-invaded
#   nerves; it only addresses cell-type identity.
#
# SCOPE (unchanged)
#   Tumour tissue of PNI-High patients only. Adjacent-normal samples are
#   not paired with tumour in this cohort, and the PNI-Low Schwann cells
#   come almost entirely from one patient (diagnostic tables below), so
#   neither is used as a contrast.
#
# STATISTICS (revised: patients, not cells, are the replicates)
#   Previously: Seurat FindMarkers (Wilcoxon) one-vs-rest, which treats
#   every cell as an independent replicate and gives far too small
#   p-values. Now:
#   1. Counts are summed per PATIENT x CELL TYPE (pseudobulk); a group is
#      used if it has >= MIN_CELLS cells.
#   2. limma-voom on all genes (library-size normalisation needs them),
#      design ~ 0 + cell type + patient (patient = blocking factor).
#   3. For EVERY cell type with >= MIN_PATIENTS patients (not only the top
#      one, so there is no post-selection): DRAXIN in that cell type vs the
#      average of all other cell types, moderated t, BH across cell types.
#      Significant = padj < ALPHA and logFC > 0. (Mean expression level.)
#   3b. PRIMARY for cell identity -- within-patient detection fraction: for
#      every cell type and patient with >= MIN_CELLS_DET cells of that type,
#      % of its cells with DRAXIN detected vs % among all other cells of the
#      same patient; paired one-sided Wilcoxon across patients, BH across
#      cell types. Needed because Schwann cells are rare after dissociation
#      (>= 20 cells in only one PNI-High tumour sample), so the pseudobulk
#      test cannot assess them; the same rule is applied to every cell type.
#   4. Descriptive values are computed per patient and then averaged with
#      equal weight per patient (% expressing; mean expression in positive
#      cells), so a patient with many cells does not dominate. Per-patient
#      values are written out (DRAXIN_celltype_per_patient.csv).
#
# WHAT THIS DOES NOT CLAIM
#   - No spatial / interface localisation of DRAXIN.
#   - No PNI-High vs PNI-Low, and no tumour vs adjacent-normal difference.
#   - A subset of CD8+ T cells also expresses DRAXIN -- reported, not pursued.
#
# INPUTS: input/sc_sn_merged_obj.Rds, input/metadata.csv
# OUTPUTS (output/03_DRAXIN_Chen_reference/):
#   tables/DRAXIN_celltype_summary.csv              main table (per cell type, both tests)
#   tables/DRAXIN_detection_per_patient_vs_other_cells.csv
#   tables/DRAXIN_celltype_per_patient.csv          per patient x cell type
#   figures/DRAXIN_dotplot_tumor_PNIhigh.pdf        % expressing x mean expression
#   figures/DRAXIN_violin_combined_tumor_PNIhigh.pdf violin of cells, % expressing (mean of patients)
#   figures/DRAXIN_per_patient_tumor_PNIhigh.pdf     % expressing per patient (points) and mean
# ============================================================

library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)
library(limma)
library(edgeR)
library(ggplot2)

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
ANALYSIS_NAME <- "03_DRAXIN_Chen_reference"
TABLES_DIR    <- file.path(PROJECT_DIR, "output", ANALYSIS_NAME, "tables")
FIGURES_DIR   <- file.path(PROJECT_DIR, "output", ANALYSIS_NAME, "figures")
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

GENE          <- "DRAXIN"
CELLTYPE_COL  <- "cell_type_refined"
PATIENT_COL   <- "patients"
MIN_CELLS     <- 20      # cells per patient x cell type for a pseudobulk group
MIN_PATIENTS  <- 2       # patients with a group, for a cell type to be tested
MIN_CELLS_DET <- 5       # cells per patient x cell type for the detection-fraction analysis
# (Schwann cells are rare after dissociation: >= 20 cells in only
#  1 PNI-High tumour sample, so the pseudobulk test cannot assess them)
ALPHA         <- 0.05

# Publication PDFs (Helvetica, editable text, no dingbats)
pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")

# ---- 1. Reference and clinical metadata ------------------------------------
if (!exists("merged_obj")) merged_obj <- readRDS(file.path(INPUT_DIR, "sc_sn_merged_obj.Rds"))
DefaultAssay(merged_obj) <- "RNA"
stopifnot(GENE %in% rownames(merged_obj))

clinical <- read.csv(file.path(INPUT_DIR, "metadata.csv")) %>%
  mutate(patient_id_std = str_replace(PatientID, "^PA#", "PA"))
meta <- merged_obj[[]]
cat("Patients in the object but not in metadata.csv (expect none):\n")
print(setdiff(unique(meta[[PATIENT_COL]]), clinical$patient_id_std))
meta$PNI_status <- unname(setNames(clinical$NI.group, clinical$patient_id_std)[meta[[PATIENT_COL]]])
if (all(is.na(meta$PNI_status))) stop("No patient matched metadata.csv -- check the patient ID format.")
cat("\nPNI status of cells (NA = patient not in metadata.csv):\n"); print(table(meta$PNI_status, useNA = "always"))

counts_all <- tryCatch(LayerData(merged_obj, assay = "RNA", layer = "counts"),
                       error = function(e) GetAssayData(merged_obj, assay = "RNA", slot = "counts"))
expr_all <- tryCatch(LayerData(merged_obj, assay = "RNA", layer = "data"),
                     error = function(e) GetAssayData(merged_obj, assay = "RNA", slot = "data"))

# ---- 2. Diagnostics behind the scope (not manuscript figures) ---------------
sch <- data.frame(expr = as.numeric(expr_all[GENE, ]), cell_type = meta[[CELLTYPE_COL]],
                  PNI_status = meta$PNI_status, patient = meta[[PATIENT_COL]]) %>%
  filter(cell_type == "Schwann", !is.na(PNI_status))
cat("\nSchwann cells by PNI status (motivates excluding PNI-Low):\n")
print(sch %>% group_by(PNI_status) %>%
        summarise(n_cells = n(), n_patients = n_distinct(patient), pct_expressing = round(100 * mean(expr > 0), 2)))
cat("\nSchwann cells per patient:\n")
print(sch %>% group_by(patient, PNI_status) %>%
        summarise(n_cells = n(), pct_expressing = round(100 * mean(expr > 0), 2), .groups = "drop") %>%
        arrange(PNI_status, desc(n_cells)))

# ---- 3. Scope: tumour tissue, PNI-High patients ------------------------------
keep <- meta$tissue %in% "Tumor" & meta$PNI_status %in% "High"
m    <- meta[keep, , drop = FALSE]
cnt  <- counts_all[, keep, drop = FALSE]
ex   <- as.numeric(expr_all[GENE, keep])
cat("\nCells in scope:", sum(keep), "| patients:", length(unique(m[[PATIENT_COL]])), "\n")

# ---- 4. Per patient x cell type: descriptive values --------------------------
cells <- data.frame(expr = ex, cell_type = m[[CELLTYPE_COL]], patient = m[[PATIENT_COL]])
per_patient <- cells %>% group_by(cell_type, patient) %>%
  summarise(n_cells = n(), pct_expressing = 100 * mean(expr > 0),
            mean_expr_pos = if (any(expr > 0)) mean(expr[expr > 0]) else NA_real_, .groups = "drop") %>%
  mutate(used_in_pseudobulk = n_cells >= MIN_CELLS, used_in_detection = n_cells >= MIN_CELLS_DET)
write.csv(per_patient, file.path(TABLES_DIR, "DRAXIN_celltype_per_patient.csv"), row.names = FALSE)

desc <- per_patient %>% filter(used_in_detection) %>% group_by(cell_type) %>%
  summarise(n_patients = n(), n_cells = sum(n_cells),
            pct_expressing_patient_mean = mean(pct_expressing),
            pct_expressing_min = min(pct_expressing), pct_expressing_max = max(pct_expressing),
            mean_expr_pos_patient_mean = mean(mean_expr_pos, na.rm = TRUE), .groups = "drop")

# ---- 5. Pseudobulk limma-voom, patient as blocking factor --------------------
grp <- paste(cells$patient, cells$cell_type, sep = "||")
n_g <- table(grp); ok <- names(n_g)[n_g >= MIN_CELLS]
sel <- grp %in% ok
gf  <- factor(grp[sel], levels = ok)
ind <- sparseMatrix(i = seq_along(gf), j = as.integer(gf), x = 1, dims = c(length(gf), nlevels(gf)))
pb  <- as.matrix(cnt[, sel, drop = FALSE] %*% ind); colnames(pb) <- ok
gi  <- do.call(rbind, strsplit(ok, "||", fixed = TRUE))
groups <- data.frame(group = ok, patient = gi[, 1], cell_type = gi[, 2])
cat("\nPseudobulk groups (>= ", MIN_CELLS, " cells), cell type x patient:\n", sep = ""); print(table(groups$cell_type, groups$patient))

dge <- DGEList(pb)
keep_g <- filterByExpr(dge, group = groups$cell_type) | rownames(pb) == GENE
dge <- calcNormFactors(dge[keep_g, , keep.lib.sizes = FALSE])
ct  <- factor(make.names(groups$cell_type)); pat <- factor(groups$patient)
design <- model.matrix(~ 0 + ct + pat); colnames(design) <- make.names(colnames(design))
fit <- lmFit(voom(dge, design), design)

pb_n <- table(groups$cell_type)
testable <- names(pb_n)[pb_n >= MIN_PATIENTS]
lv <- levels(ct); ct_cols <- paste0("ct", lv)
C <- sapply(testable, function(t) {
  v <- setNames(numeric(ncol(design)), colnames(design)); a <- paste0("ct", make.names(t))
  v[a] <- 1; v[setdiff(ct_cols, a)] <- -1 / (length(lv) - 1); v })
f2 <- eBayes(contrasts.fit(fit, C))
tests <- data.frame(cell_type = testable,
                    logFC_vs_other_types = f2$coefficients[GENE, ],
                    t = f2$t[GENE, ], p_value = f2$p.value[GENE, ]) %>%
  mutate(padj = p.adjust(p_value, "BH"),
         significant = padj < ALPHA & logFC_vs_other_types > 0)

# ---- 5b. Within-patient detection fraction (primary for cell identity) -------
# For every cell type and every patient with >= MIN_CELLS_DET cells of that
# type: % of its cells with DRAXIN detected vs % among all OTHER cells of the
# same patient. Paired one-sided Wilcoxon signed-rank test across patients
# (patients = replicates), BH across cell types. Same rule for every type.
tot <- cells %>% group_by(patient) %>% summarise(n_all = n(), pos_all = sum(expr > 0), .groups = "drop")
det <- cells %>% group_by(cell_type, patient) %>% summarise(n = n(), pos = sum(expr > 0), .groups = "drop") %>%
  filter(n >= MIN_CELLS_DET) %>% left_join(tot, by = "patient") %>%
  mutate(pct_type = 100 * pos / n, pct_other = 100 * (pos_all - pos) / (n_all - n), diff = pct_type - pct_other)
det_tests <- det %>% group_by(cell_type) %>%
  summarise(det_n_patients = n(), det_n_higher = sum(diff > 0),
            det_median_diff = median(diff),
            det_p = if (n() >= 3) suppressWarnings(wilcox.test(pct_type, pct_other, paired = TRUE,
                                                               alternative = "greater", exact = FALSE)$p.value) else NA_real_,
            .groups = "drop") %>%
  mutate(det_padj = p.adjust(det_p, "BH"), det_significant = !is.na(det_padj) & det_padj < ALPHA)
write.csv(det, file.path(TABLES_DIR, "DRAXIN_detection_per_patient_vs_other_cells.csv"), row.names = FALSE)

summary_final <- desc %>% left_join(det_tests, by = "cell_type") %>% left_join(tests, by = "cell_type") %>%
  mutate(pseudobulk_tested = cell_type %in% testable) %>% arrange(desc(pct_expressing_patient_mean))
write.csv(summary_final, file.path(TABLES_DIR, "DRAXIN_celltype_summary.csv"), row.names = FALSE)
cat("\nDRAXIN per cell type (patients as replicates; BH across cell types).\n",
    "det_*: % detected vs other cells of the same patient (>= ", MIN_CELLS_DET, " cells);\n",
    "logFC/padj: pseudobulk limma-voom (>= ", MIN_CELLS, " cells; types with < ", MIN_PATIENTS, " patients not tested):\n", sep = "")
print(as.data.frame(summary_final), row.names = FALSE, digits = 3)

# ---- 6. Figures ------------------------------------------------------------
order_lvls <- summary_final$cell_type
n_pat <- length(unique(cells$patient))
cap <- paste0("Tumour tissue, PNI-High patients (Chen et al. 2025); ", n_pat, " patients; patient x cell type groups with >= ",
              MIN_CELLS_DET, " cells. Values averaged over patients (equal weight).\n* = detected in a higher fraction of cells ",
              "than the other cells of the same patient (paired Wilcoxon across patients, BH across cell types, padj < ", ALPHA, ").")
theme_pub <- theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 7))

dot_df <- summary_final %>% mutate(cell_type = factor(cell_type, order_lvls),
                                   star = ifelse(det_significant %in% TRUE, "*", ""))
p_dot <- ggplot(dot_df, aes(cell_type, GENE, size = pct_expressing_patient_mean, colour = mean_expr_pos_patient_mean)) +
  geom_point() +
  geom_text(aes(label = star), colour = "black", size = 6, vjust = -1.1, show.legend = FALSE) +
  scale_size(range = c(0, 10), name = "% expressing\n(mean of patients)") +
  scale_colour_gradient(low = "lightgrey", high = "green4", name = "Mean expression\n(positive cells)", na.value = "grey90") +
  labs(title = paste(GENE, "expression across cell types"), x = NULL, y = NULL, caption = cap) +
  theme_pub + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
ggsave(file.path(FIGURES_DIR, "DRAXIN_dotplot_tumor_PNIhigh.pdf"), p_dot, width = 9, height = 4.2, device = pdf_device)

# Violin + jitter of all cells (descriptive, as in the original figure); the
# % label above each violin is the patient-averaged % expressing, and the
# asterisk comes from the within-patient detection test (patients as replicates).
vio_df <- cells %>% filter(cell_type %in% order_lvls) %>% mutate(cell_type = factor(cell_type, order_lvls))
lab_df <- dot_df %>% mutate(y = max(vio_df$expr) * 1.12,
                            label = paste0(round(pct_expressing_patient_mean, 1), "%", star))
p_vio <- ggplot(vio_df, aes(cell_type, expr)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, fill = "grey80") +
  geom_jitter(height = 0, width = 0.15, size = 0.3, alpha = 0.4, colour = "black") +
  geom_text(data = lab_df, aes(cell_type, y, label = label), inherit.aes = FALSE, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(title = paste(GENE, "expression across cell types"), x = NULL, y = "Expression level",
       caption = paste0("Violins: all cells (descriptive). Labels: % expressing cells, mean of patients. ", cap)) +
  theme_pub + theme(legend.position = "none")
ggsave(file.path(FIGURES_DIR, "DRAXIN_violin_combined_tumor_PNIhigh.pdf"), p_vio, width = 11, height = 6, device = pdf_device)

# Per-patient values: one point per patient (groups with >= MIN_CELLS_DET
# cells), black bar = mean of patients, asterisk from the detection test.
pp_df <- per_patient %>% filter(used_in_detection, cell_type %in% order_lvls) %>%
  mutate(cell_type = factor(cell_type, order_lvls))
star_df <- dot_df %>% filter(star == "*") %>% mutate(y = max(pp_df$pct_expressing) * 1.08)
p_pat <- ggplot(pp_df, aes(cell_type, pct_expressing)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "black") +
  geom_point(shape = 21, size = 2.6, stroke = 0.4, fill = "grey35", colour = "white",
             position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  geom_text(data = star_df, aes(cell_type, y, label = star), size = 6, inherit.aes = FALSE) +
  labs(title = paste(GENE, "expression across cell types"), x = NULL, y = "% of cells expressing",
       caption = paste0("Points: individual patients; bar: mean of patients. ", cap)) +
  theme_pub
ggsave(file.path(FIGURES_DIR, "DRAXIN_per_patient_tumor_PNIhigh.pdf"), p_pat, width = 8, height = 5, device = pdf_device)
# Patient x cell type overview: colour = % of cells with DRAXIN detected,
# dot size = number of cells; groups below MIN_CELLS_DET cells (not used in the
# tests) are drawn as open grey circles, so the data behind every value is visible.
ov <- per_patient %>% filter(cell_type %in% order_lvls) %>%
  mutate(cell_type = factor(cell_type, rev(order_lvls)),
         patient = factor(patient, sort(unique(patient))))
p_ov <- ggplot(ov, aes(patient, cell_type)) +
  geom_point(data = filter(ov, !used_in_detection), aes(size = n_cells), shape = 21, colour = "grey60", fill = NA, stroke = 0.4) +
  geom_point(data = filter(ov, used_in_detection), aes(size = n_cells, fill = pct_expressing), shape = 21, colour = "grey30", stroke = 0.3) +
  scale_y_discrete(limits = levels(ov$cell_type)) + scale_x_discrete(limits = levels(ov$patient)) +
  scale_fill_gradient(low = "white", high = "green4", name = "% cells\nexpressing", limits = c(0, NA)) +
  scale_size_area(max_size = 7, trans = "sqrt", breaks = c(5, 20, 100, 500, 2000), name = "Cells") +
  labs(title = paste(GENE, "-- detection per patient and cell type"), x = "Patient", y = NULL,
       caption = paste0("Tumour tissue, PNI-High patients (Chen et al. 2025). Filled: groups with >= ", MIN_CELLS_DET,
                        " cells (used in the tests); open grey: fewer cells.")) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 7), panel.grid.major = element_line(colour = "grey93", linewidth = 0.3))
ggsave(file.path(FIGURES_DIR, "DRAXIN_patient_by_celltype_overview.pdf"), p_ov,
       width = 3 + 0.35 * nlevels(ov$patient), height = 1.5 + 0.3 * nlevels(ov$cell_type), device = pdf_device)
cat("\nSaved tables and figures in", file.path(PROJECT_DIR, "output", ANALYSIS_NAME), "\n")

