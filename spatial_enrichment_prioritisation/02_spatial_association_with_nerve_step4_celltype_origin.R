# ======================================================================
# 02_spatial_association_with_nerve_step4_celltype_origin.R
#
# Figure 3D: which cell type most likely expresses each signature gene, in
# an independent sc/snRNA-seq reference (Chen et al.), tumour tissue of
# PNI-high patients.
#
# HOW IT IS CALCULATED
#   1. Cells are grouped by PATIENT x CELL TYPE; a group is used if it has
#      >= MIN_CELLS cells.
#   2. Descriptive (heatmap, "top cell type"): % of cells expressing the gene
#      (count > 0) in each group, then averaged over patients with EQUAL
#      weight -- a patient with many cells does not dominate. The top cell
#      type is the one with the highest patient-averaged %; the gap to the
#      second is reported as specificity.
#   3. Test (replaces the earlier cell-level FindMarkers, which treated every
#      cell as an independent replicate): counts are summed per patient x
#      cell type (pseudobulk) and analysed with limma-voom, design
#      ~ 0 + cell type + patient (patient = blocking factor; optionally the
#      technology, sc vs sn). For each gene, its top cell type is tested
#        (a) against the average of all other cell types  -> "enriched", and
#        (b) against the runner-up cell type               -> "specific",
#      one-sided (higher in the top cell type), BH across the signature genes.
#      The replicates are patients, not cells.
#   4. Summary: number of genes per top cell type, split by whether the
#      assignment is significant (a); Schwann cells highlighted.
#   Caveat: the top cell type is chosen from the same data that are then
#   tested, so the p-values (especially "vs runner-up") are somewhat
#   optimistic; the assignment is supporting evidence of the likely source,
#   not proof of origin.
#
# INPUTS: input/sc_sn_merged_obj.Rds, input/metadata.csv,
#         tables/lee_sham_signature_core.csv (step1)
# OUTPUTS: tables/fig3D_celltype_*.csv, figures/fig3D_celltype_*.pdf
# ======================================================================

library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)
library(limma)
library(edgeR)
library(ggplot2)
library(pheatmap)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
summarise <- dplyr::summarise; group_by <- dplyr::group_by; arrange <- dplyr::arrange
count <- dplyr::count; rename <- dplyr::rename

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input"))
  stop("Run this script from the project folder (the one containing input/).")

PROJECT_DIR   <- "."
INPUT_DIR    <- file.path(PROJECT_DIR, "input")
TABLES_DIR   <- file.path(PROJECT_DIR, "output", "02_spatial_association", "tables")
FIGURES_DIR  <- file.path(PROJECT_DIR, "output", "02_spatial_association", "figures")
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

GENE_SET_FILE <- file.path(TABLES_DIR, "lee_sham_signature_core.csv")
CELLTYPE_COL  <- "cell_type_refined"
PATIENT_COL   <- "patients"
TECH_COL      <- NULL          # e.g. "technology" if the object mixes sc and sn per patient; NULL = not used
MIN_CELLS     <- 20            # cells per patient x cell type
ALPHA         <- 0.05
# Publication PDFs: base pdf device, standard PostScript font (Helvetica) kept as
# editable text, no dingbat glyphs for points -- opens cleanly in Illustrator /
# Acrobat (cairo_pdf embeds Linux fonts such as DejaVu that Adobe does not know).
pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")

# Standard cell-type palette of the paper (color_maping used for the earlier figures)
color_mapping <- c(
  "Acinar" = "gold1", "B_cells" = "maroon", "Lymphocyte(CD4)" = "green4",
  "Lymphocytes" = "green4", "Myeloid(Macrophages)" = "#6A3D9A",
  "Endothelial" = "#FF7F00", "Endocrine" = "yellow4", "Stellate" = "khaki2",
  "CAF" = "skyblue2", "Adipocytes" = "steelblue4", "Lymphocyte(NK)" = "green2",
  "Lymphocyte(CD8)" = "palegreen2", "Myeloid" = "#CAB2D6", "Schwann" = "#FDBF6F",
  "Myeloid(C1q_macrophages)" = "#CAB2D6", "Nerve" = "#FDBF6F", "Muscle" = "orchid1",
  "Ductal_malignant" = "gray70", "Ductal_normal" = "black",
  "B_cells_(Plasma cells)" = "brown1"
)
# Reference cell-type names that differ from the palette names are matched by
# keyword (first match wins); the mapping is printed so it can be checked.
# Anything unmatched gets a neutral grey and is reported.
palette_alias <- c("schwann|glia" = "Schwann", "plasma" = "B_cells_(Plasma cells)", "b.?cell" = "B_cells",
                   "cd8" = "Lymphocyte(CD8)", "cd4" = "Lymphocyte(CD4)", "nk" = "Lymphocyte(NK)",
                   "t.?cell|lymph" = "Lymphocytes", "c1q" = "Myeloid(C1q_macrophages)",
                   "macroph" = "Myeloid(Macrophages)", "myeloid|monocyt|dendrit|mast|neutro" = "Myeloid",
                   "endothel" = "Endothelial", "endocrin|islet" = "Endocrine", "stellate|pericyte" = "Stellate",
                   "caf|fibro" = "CAF", "adipo" = "Adipocytes", "muscle|smooth" = "Muscle", "acinar" = "Acinar",
                   "normal.*duct|duct.*normal" = "Ductal_normal",
                   "duct|epithel|malignan|tumou?r|cancer" = "Ductal_malignant")
celltype_colours <- function(types) {
  key <- vapply(types, function(t) {
    if (t %in% names(color_mapping)) return(t)
    hit <- names(palette_alias)[vapply(names(palette_alias), function(a) grepl(a, t, ignore.case = TRUE), logical(1))]
    if (length(hit)) palette_alias[[hit[1]]] else NA_character_
  }, character(1))
  cols <- ifelse(is.na(key), NA_character_, color_mapping[key])
  if (any(is.na(cols))) {
    cat("Cell types without a palette colour (grey used):", paste(types[is.na(cols)], collapse = ", "), "\n")
    cols[is.na(cols)] <- grDevices::gray.colors(sum(is.na(cols)), start = 0.35, end = 0.85)
  }
  cat("Cell-type colours:\n"); print(data.frame(cell_type = types, palette_entry = key, colour = unname(cols)), row.names = FALSE)
  setNames(unname(cols), types)
}

# ---- 1. Reference and scope ------------------------------------------------
if (!exists("merged_obj")) merged_obj <- readRDS(file.path(INPUT_DIR, "sc_sn_merged_obj.Rds"))
meta <- merged_obj[[]]
clinical <- read.csv(file.path(INPUT_DIR, "metadata.csv")) %>%
  mutate(patient_id_std = str_replace(PatientID, "^PA#", "PA"))
meta$PNI_status <- unname(setNames(clinical$NI.group, clinical$patient_id_std)[meta[[PATIENT_COL]]])
unmatched <- unique(meta[[PATIENT_COL]][is.na(meta$PNI_status)])
if (length(unmatched)) cat("Patients without PNI status in metadata.csv (excluded):", paste(unmatched, collapse = ", "), "\n")
if (all(is.na(meta$PNI_status))) stop("No patient matched metadata.csv -- check the patient ID format.")
keep <- meta$tissue %in% "Tumor" & meta$PNI_status %in% "High"        # NA-safe
cat("Cells in scope (tumour, PNI-high):", sum(keep), "| patients:", length(unique(meta[[PATIENT_COL]][keep])), "\n")

genes_sig <- read.csv(GENE_SET_FILE)$Gene
counts <- tryCatch(LayerData(merged_obj, assay = "RNA", layer = "counts"),
                   error = function(e) GetAssayData(merged_obj, assay = "RNA", slot = "counts"))
counts <- counts[, keep, drop = FALSE]
m <- meta[keep, , drop = FALSE]
genes_in <- intersect(genes_sig, rownames(counts))
cat("Signature genes:", length(genes_sig), "| in reference:", length(genes_in), "\n")
if (length(genes_in) < length(genes_sig))
  cat("Not in reference:", paste(setdiff(genes_sig, genes_in), collapse = ", "), "\n")

# ---- 2. Patient x cell-type groups ------------------------------------------
grp <- paste(m[[PATIENT_COL]], m[[CELLTYPE_COL]], sep = "||")
if (!is.null(TECH_COL)) grp <- paste(grp, m[[TECH_COL]], sep = "||")
n_cells <- table(grp)
ok_grp <- names(n_cells)[n_cells >= MIN_CELLS]
cells_ok <- grp %in% ok_grp
grp_f <- factor(grp[cells_ok], levels = ok_grp)
ind <- sparseMatrix(i = seq_along(grp_f), j = as.integer(grp_f), x = 1,
                    dims = c(length(grp_f), nlevels(grp_f)), dimnames = list(NULL, levels(grp_f)))
cnt <- counts[, cells_ok, drop = FALSE]
gi <- do.call(rbind, strsplit(levels(grp_f), "||", fixed = TRUE))
groups <- data.frame(group = levels(grp_f), patient = gi[, 1], cell_type = gi[, 2],
                     tech = if (!is.null(TECH_COL)) gi[, 3] else NA, n_cells = as.integer(n_cells[levels(grp_f)]))
cat("Groups used (>= ", MIN_CELLS, " cells):\n", sep = ""); print(table(groups$cell_type, groups$patient))

# ---- 3. Descriptive: % expressing, averaged over patients -------------------
pct_grp <- as.matrix((cnt[genes_in, , drop = FALSE] > 0) %*% ind) %*% diag(1 / groups$n_cells) * 100
colnames(pct_grp) <- groups$group
pct <- sapply(sort(unique(groups$cell_type)), function(ct)
  rowMeans(pct_grp[, groups$cell_type == ct, drop = FALSE]))      # equal weight per patient
pct <- matrix(pct, nrow = length(genes_in), dimnames = list(genes_in, sort(unique(groups$cell_type))))
n_pat_ct <- table(groups$cell_type)
ord <- t(apply(pct, 1, function(x) order(x, decreasing = TRUE)))
desc <- data.frame(Gene = genes_in,
                   top_cell_type = colnames(pct)[ord[, 1]], top_pct = pct[cbind(seq_along(genes_in), ord[, 1])],
                   second_cell_type = colnames(pct)[ord[, 2]], second_pct = pct[cbind(seq_along(genes_in), ord[, 2])]) %>%
  mutate(specificity_gap = top_pct - second_pct,
         n_patients_top = as.integer(n_pat_ct[top_cell_type]))
write.csv(cbind(Gene = genes_in, as.data.frame(round(pct, 2))),
          file.path(TABLES_DIR, "fig3D_celltype_pct_expressing_patient_mean.csv"), row.names = FALSE)

# ---- 4. Pseudobulk limma-voom, patient as blocking factor -------------------
pb <- as.matrix(cnt %*% ind)
dge <- DGEList(pb)
keep_g <- filterByExpr(dge, group = groups$cell_type) | rownames(pb) %in% genes_in
dge <- calcNormFactors(dge[keep_g, , keep.lib.sizes = FALSE])
ct  <- factor(make.names(groups$cell_type)); pat <- factor(groups$patient)
design <- if (!is.null(TECH_COL)) model.matrix(~ 0 + ct + pat + factor(groups$tech)) else model.matrix(~ 0 + ct + pat)
colnames(design) <- make.names(colnames(design))
v   <- voom(dge, design)
fit <- lmFit(v, design)
lv  <- levels(ct); ct_cols <- paste0("ct", lv)

contrast_vectors <- function(top, second) {
  a <- make.names(top); b <- make.names(second)
  c_rest <- setNames(numeric(ncol(design)), colnames(design))
  c_rest[paste0("ct", a)] <- 1; c_rest[setdiff(ct_cols, paste0("ct", a))] <- -1 / (length(lv) - 1)
  c_sec <- setNames(numeric(ncol(design)), colnames(design)); c_sec[paste0("ct", a)] <- 1; c_sec[paste0("ct", b)] <- -1
  list(c_rest = c_rest, c_sec = c_sec)
}
# One contrast fit per (top, second) combination, then pick each gene's row
combos <- unique(desc[, c("top_cell_type", "second_cell_type")])
res_list <- lapply(seq_len(nrow(combos)), function(k) {
  cc <- contrast_vectors(combos$top_cell_type[k], combos$second_cell_type[k])
  C <- cbind(rest = cc$c_rest, second = cc$c_sec)
  f2 <- eBayes(contrasts.fit(fit, C))
  gs <- desc$Gene[desc$top_cell_type == combos$top_cell_type[k] & desc$second_cell_type == combos$second_cell_type[k]]
  gs <- intersect(gs, rownames(f2$coefficients))
  df <- f2$df.total[match(gs, rownames(f2$coefficients))]
  data.frame(Gene = gs,
             logFC_vs_rest = f2$coefficients[gs, "rest"], t_rest = f2$t[gs, "rest"],
             logFC_vs_second = f2$coefficients[gs, "second"], t_second = f2$t[gs, "second"],
             p_vs_rest = pt(f2$t[gs, "rest"], df, lower.tail = FALSE),        # one-sided: higher in top type
             p_vs_second = pt(f2$t[gs, "second"], df, lower.tail = FALSE))
})
tests <- bind_rows(res_list) %>%
  mutate(padj_vs_rest = p.adjust(p_vs_rest, "BH"), padj_vs_second = p.adjust(p_vs_second, "BH"))

final <- desc %>% left_join(tests, by = "Gene") %>%
  mutate(enriched = !is.na(padj_vs_rest) & padj_vs_rest < ALPHA & logFC_vs_rest > 0,
         specific = enriched & !is.na(padj_vs_second) & padj_vs_second < ALPHA & logFC_vs_second > 0,
         # Assignment categories:
         #   Specific: significantly higher than the average of the other cell types AND
         #             than the next-highest cell type
         #   Enriched: significantly higher than the average of the other cell types only
         #   Highest % (not significant): top cell type by %, no significant difference
         #   Not assignable: no significant difference AND < 5% of cells express it
         #             even in the top cell type -> the reference cannot tell its origin;
         #             these genes are NOT counted under any cell type
         assignment = case_when(specific ~ "Specific",
                                enriched ~ "Enriched",
                                top_pct < 5 ~ "Not assignable",
                                TRUE ~ "Highest % (not significant)")) %>%
  mutate(assigned_cell_type = ifelse(assignment == "Not assignable", NA_character_, top_cell_type)) %>%
  arrange(top_cell_type, desc(specificity_gap))
write.csv(final, file.path(TABLES_DIR, "fig3D_celltype_per_gene.csv"), row.names = FALSE)

NOT_ASSIGNABLE_BAR <- "Not assignable\n(< 5% of cells)"
SHOW_NOT_ASSIGNABLE_BAR <- FALSE   # FALSE: unassignable genes are left out of the bar plot and named in the
# console (for the legend); TRUE: shown as their own bar
not_assignable <- final$Gene[final$assignment == "Not assignable"]
cat("\nNot assignable (< 5% of cells in every cell type, no significant difference):", length(not_assignable),
    if (length(not_assignable)) paste0("-- ", paste(sort(not_assignable), collapse = ", ")) else "", "\n")
summary_tab <- final %>%
  mutate(top_cell_type = ifelse(assignment == "Not assignable", NOT_ASSIGNABLE_BAR, top_cell_type)) %>%
  count(top_cell_type, assignment, name = "n_genes")
plot_tab <- if (SHOW_NOT_ASSIGNABLE_BAR) summary_tab else filter(summary_tab, assignment != "Not assignable")
write.csv(summary_tab, file.path(TABLES_DIR, "fig3D_celltype_summary.csv"), row.names = FALSE)
cat("\nGenes per top cell type (patients as replicates):\n")
print(final %>% count(top_cell_type, assignment) %>% tidyr::pivot_wider(names_from = assignment, values_from = n, values_fill = 0))
cat("\nSchwann: ", sum(final$top_cell_type == "Schwann"), " of ", nrow(final), " genes have Schwann as top cell type; ",
    sum(final$top_cell_type == "Schwann" & final$enriched), " significantly enriched, ",
    sum(final$top_cell_type == "Schwann" & final$specific), " specific.\n", sep = "")

# ---- 5. Figures -----------------------------------------------------------
lev_order <- summary_tab %>% filter(top_cell_type != NOT_ASSIGNABLE_BAR) %>%
  group_by(top_cell_type) %>% summarise(n = sum(n_genes)) %>%
  arrange(desc(n)) %>% pull(top_cell_type)
bar_order <- c(lev_order, intersect(NOT_ASSIGNABLE_BAR, summary_tab$top_cell_type))

asg_col <- c("Specific" = "#08519C", "Enriched" = "#6BAED6",
             "Highest % (not significant)" = "grey80", "Not assignable" = "grey92")

p_sum <- ggplot(plot_tab %>%
                  mutate(top_cell_type = factor(top_cell_type,
                                                levels = rev(intersect(bar_order, top_cell_type))),
                         assignment = factor(assignment, names(asg_col))),
                aes(n_genes, top_cell_type, fill = assignment)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = asg_col[intersect(names(asg_col), plot_tab$assignment)],
                    name = NULL) +
  guides(fill = guide_legend(ncol = 1, keyheight = unit(3, "mm"))) +
  labs(y = NULL,
       x = paste0("Spatially prioritised genes (", sum(plot_tab$n_genes),
                  " of ", nrow(final), ")")) +
  theme_classic(base_size = 7) +
  theme(axis.text  = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.margin = margin(t = 0),
        plot.margin = margin(2, 3, 2, 2, "mm"))

ggsave(file.path(FIGURES_DIR, "fig3D_celltype_summary.pdf"), p_sum,
       width = 50, height = 85, units = "mm", device = pdf_device)

# ---- Heatmaps of % cells expressing (mean over patients) --------------------
# Columns in one fixed order (clustered once on all genes) and one colour scale
# (0-100 %) for every heatmap, so they are directly comparable.
#   fig3D_celltype_heatmap_all_genes.pdf -- all genes, rows clustered without
#     supervision, no gene names; the annotation bars (top cell type,
#     assignment) show whether the clustering reproduces the assignment.
#     Gene-level values: fig3D_celltype_per_gene.csv (supplementary table).
#   fig3D_celltype_heatmap_<cell type>.pdf -- named heatmaps for the cell types
#     forming the invasion site (HEATMAP_CELLTYPES: nerve sheath, vessels,
#     tumour): every gene whose top cell type is that cell type, whatever its
#     assignment (shown in the annotation bar), clustered -- except genes that
#     are "Not assignable" (< 5% of cells in every cell type), which belong to no
#     cell type and appear only in the all-genes heatmap.
#   Per-assignment heatmaps only if HEATMAPS_BY_ASSIGNMENT = TRUE.
HEATMAP_CELLTYPES <- c(Schwann = "schwann", Endothelial = "endothel", Ductal_malignant = "duct.*malig|malig.*duct")
HEATMAPS_BY_ASSIGNMENT <- FALSE
# This section can be run on its own after section 0 (settings, palette,
# pdf_device): if the results of sections 1-5 are not in memory, it reads them
# back from the tables written above, so the reference does not need reloading.
if (!exists("final") || !exists("pct")) {
  final <- read.csv(file.path(TABLES_DIR, "fig3D_celltype_per_gene.csv"), check.names = FALSE)
  pct_tab <- read.csv(file.path(TABLES_DIR, "fig3D_celltype_pct_expressing_patient_mean.csv"), check.names = FALSE)
  pct <- as.matrix(pct_tab[, setdiff(names(pct_tab), "Gene"), drop = FALSE]); rownames(pct) <- pct_tab$Gene
  cat("Heatmaps: read", nrow(final), "genes x", ncol(pct), "cell types from the step4 tables.\n")
}
# Tables written by an earlier version of this script use the original
# category names; translate them to the current ones (no change in meaning).
OLD_ASSIGNMENT_NAMES <- c("Specific (vs runner-up)"        = "Specific",
                          "Enriched (vs other cell types)" = "Enriched",
                          "Highest %, not significant"     = "Highest % (not significant)",
                          "Barely detected (< 5 % cells)"  = "Not assignable",
                          "Barely detected"                = "Not assignable")
is_old <- final$assignment %in% names(OLD_ASSIGNMENT_NAMES)
final$assignment[is_old] <- unname(OLD_ASSIGNMENT_NAMES[final$assignment[is_old]])
if (!exists("lev_order"))
  lev_order <- final %>% filter(assignment != "Not assignable") %>% count(top_cell_type) %>%
  arrange(desc(n)) %>% pull(top_cell_type) %>% union(unique(final$top_cell_type))
if (!exists("ct_col"))
  ct_col <- celltype_colours(sort(unique(c(final$top_cell_type, colnames(pct)))))
if (!exists("asg_col"))
  asg_col <- c("Specific" = "#08519C", "Enriched" = "#6BAED6",
               "Highest % (not significant)" = "grey80", "Not assignable" = "grey92")
stopifnot(all(final$Gene %in% rownames(pct)))
unknown <- setdiff(unique(final$assignment), names(asg_col))
if (length(unknown)) stop("Unknown assignment categories in the table: ", paste(unknown, collapse = ", "))

HM_BREAKS <- seq(0, 100, length.out = 101)
HM_COLS   <- colorRampPalette(c("white", "green4"))(100)
col_order <- colnames(pct)[hclust(dist(t(pct[final$Gene, , drop = FALSE])))$order]
order_rows <- function(genes, by_celltype = TRUE) {
  if (!by_celltype) return(if (length(genes) > 2) genes[hclust(dist(pct[genes, col_order, drop = FALSE]))$order] else genes)
  d <- final %>% filter(Gene %in% genes) %>% mutate(top_cell_type = factor(top_cell_type, lev_order)) %>%
    arrange(top_cell_type)
  unlist(lapply(split(d$Gene, droplevels(d$top_cell_type)), function(g)
    if (length(g) > 2) g[hclust(dist(pct[g, col_order, drop = FALSE]))$order] else g), use.names = FALSE)
}
draw_heatmap <- function(genes, file, title, show_names, by_celltype = TRUE) {
  cat("Heatmap", file, ":", length(genes), "genes\n")
  if (length(genes) < 2) { cat("  skipped (a heatmap needs at least 2 genes)\n"); return(invisible(NULL)) }
  g <- order_rows(genes, by_celltype)
  ct <- as.character(final$top_cell_type[match(g, final$Gene)])
  gaps <- if (by_celltype) which(head(ct, -1) != tail(ct, -1)) else integer(0)
  ann <- data.frame(Top_cell_type = ct, Assignment = final$assignment[match(g, final$Gene)], row.names = g)
  ann_col <- data.frame(Cell_type = col_order, row.names = col_order)
  h <- if (show_names) 2.5 + 0.14 * length(g) else 6
  # Build the heatmap first (silent = TRUE: nothing drawn yet), so that an error
  # never leaves an empty, open PDF behind; then draw it into the PDF.
  ph <- tryCatch(
    pheatmap(pct[g, col_order, drop = FALSE], cluster_rows = FALSE, cluster_cols = FALSE,
             color = HM_COLS, breaks = HM_BREAKS, gaps_row = if (length(gaps)) gaps else NULL,
             annotation_row = ann, annotation_col = ann_col,
             annotation_colors = list(Top_cell_type = ct_col[unique(ct)], Cell_type = ct_col[col_order],
                                      Assignment = asg_col[intersect(names(asg_col), unique(ann$Assignment))]),
             show_rownames = show_names, fontsize = 8, fontsize_row = 7, fontsize_col = 8,
             main = paste0(title, " (n = ", length(g), "): % cells expressing, mean over patients"), silent = TRUE),
    error = function(e) { cat("  ERROR building heatmap:", conditionMessage(e), "\n"); NULL })
  if (is.null(ph)) return(invisible(NULL))
  pdf_device(file.path(FIGURES_DIR, file), width = 8, height = h)
  tryCatch({ grid::grid.newpage(); grid::grid.draw(ph$gtable) },
           error = function(e) cat("  ERROR drawing heatmap:", conditionMessage(e), "\n"),
           finally = dev.off())
  cat("  saved\n")
}
draw_heatmap(final$Gene, "fig3D_celltype_heatmap_all_genes.pdf", "All signature genes",
             show_names = FALSE, by_celltype = FALSE)
for (k in names(HEATMAP_CELLTYPES)) {
  types <- unique(final$top_cell_type[grepl(HEATMAP_CELLTYPES[[k]], final$top_cell_type, ignore.case = TRUE)])
  if (!length(types)) { cat("Heatmap", k, ": no top cell type matches '", HEATMAP_CELLTYPES[[k]], "'\n"); next }
  cat("Heatmap", k, "uses top cell type(s):", paste(types, collapse = ", "), "\n")
  # genes the reference cannot assign (< 5% of cells everywhere) are not shown under a cell type
  draw_heatmap(final$Gene[final$top_cell_type %in% types & final$assignment != "Not assignable"],
               paste0("fig3D_celltype_heatmap_", k, ".pdf"), paste0("Top cell type: ", paste(types, collapse = "/")),
               show_names = TRUE, by_celltype = FALSE)
}
if (HEATMAPS_BY_ASSIGNMENT) for (a in names(asg_col)) {
  draw_heatmap(final$Gene[final$assignment == a],
               paste0("fig3D_celltype_heatmap_", gsub("[^A-Za-z]+", "_", tolower(a)), ".pdf"), a, show_names = TRUE)
}
# Counts for the text: genes per top cell type x assignment
txt_tab <- final %>% count(top_cell_type, assignment) %>% tidyr::pivot_wider(names_from = assignment, values_from = n, values_fill = 0)
cat("\nGenes per top cell type and assignment (for the text):\n"); print(as.data.frame(txt_tab))
cat("\nSaved fig3D tables and figures.\n")
