# ======================================================================
# 02_spatial_association_with_nerve_step2_decay_examples.R
#
# Figure 3A: expression of example genes (default KRT19 and NES) as a
# function of distance from invaded nerves, compared with the same curve
# around geometry-matched relocated shapes of the same nerves (the shams of
# the plasma-cell analysis), plus a true-to-scale zoom of one typical
# invaded nerve and one of its relocated shapes.
#
# Same framework as the plasma-cell decay figure (01_plasma_decay_analysis):
#   * SAME distance tables (bundle from plasma step1), SAME bins and SAME
#     bin function (bin_mid() from 00_decay_shared.R: 6 x 150 um bins,
#     0-900 um, points at bin midpoints 75 ... 825; bin 0 = the nerve or
#     sham spots themselves); check_bundle_settings() stops if the bundle
#     was built with other settings;
#   * same aggregation and same plot style as the plasma curve figure;
#   * per nerve: curve around the nerve, and the mean curve of its own
#     relocated replicates (the control);
#   * per patient: mean over nerves; shown: mean +/- SE across patients.
#   * Illustration only: no per-gene statistic is computed here (the test
#     for all genes is step1); the curves show the shape of the profile.
# Expression: logcounts of the gene per spot.
#
# Example nerve for the zoom (same nerve for all genes): among invaded nerves
# with complete profiles, the one whose profile (observed minus own control,
# per bin) is closest to the median profile, summed over the genes after
# scaling each gene by its spread -- a typical nerve, not a clear one.
#
# Also writes, for the same genes, the per-patient Lee's L of the real masks
# and the mean / 95th percentile of their 200 relocated masks (from step1),
# for quoting in the text; no extra figure.
#
# INPUTS: processed/figure_bundle_for_real_vs_relocated.rds (plasma step1),
#   tables/lee_sham_DE_genes_scores.csv (Lee step1),
#   input/Integrated_NMF.rds (+ DE_SPE_annotaion.csv), um_per_pixel
# OUTPUTS: output/02_spatial_association/figures/fig3A_decay_<genes>.pdf,
#   tables/fig3A_decay_curves_patient_level.csv,
#   tables/fig3A_lee_real_vs_sham_example_genes.csv
# ======================================================================

library(SpatialExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggforce)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
summarise <- dplyr::summarise; group_by <- dplyr::group_by; arrange <- dplyr::arrange
distinct <- dplyr::distinct

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input"))
  stop("Run this script from the project folder (the one containing input/).")

PROJECT_DIR   <- "."
INPUT_DIR     <- file.path(PROJECT_DIR, "input")
PROCESSED_DIR <- file.path(PROJECT_DIR, "processed")
TABLES_DIR    <- file.path(PROJECT_DIR, "output", "02_spatial_association", "tables")
FIGURES_DIR   <- file.path(PROJECT_DIR, "output", "02_spatial_association", "figures")
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

GENES          <- c("KRT19", "NES")        # plotted top to bottom
GROUP          <- "Invaded nerves (Diseased)"
EXPR_ASSAY     <- "logcounts"
CROP_MARGIN_UM <- 300    # zoom: spots up to this distance around the nerve outline
NERVE_OUTLINE_COL <- "#E41A1C"
MAKE_RELOCATED_ZOOM <- TRUE   # also write the supplementary version with the relocated-shape zoom
SPOT_DIAMETER_UM <- 55
PATIENT_LABELS <- c("54703" = "P1", "6249" = "P2", "70865" = "P3")

SHARED_CANDIDATES <- file.path(c(file.path(PROJECT_DIR, "scripts"), file.path(PROJECT_DIR, "script"), PROJECT_DIR),
                               "00_decay_shared.R")
SHARED_PATH <- SHARED_CANDIDATES[file.exists(SHARED_CANDIDATES)][1]
if (is.na(SHARED_PATH)) stop("00_decay_shared.R not found.")
source(SHARED_PATH)
# Publication PDFs: base pdf device, standard PostScript font (Helvetica) kept as
# editable text, no dingbat glyphs for points -- opens cleanly in Illustrator /
# Acrobat (cairo_pdf embeds Linux fonts such as DejaVu that Adobe does not know).
pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")

# ---- 1. Data ---------------------------------------------------------------
bundle <- readRDS(file.path(PROCESSED_DIR, "figure_bundle_for_real_vs_relocated.rds"))
check_bundle_settings(bundle)
meta    <- bundle$spe_meta
coords  <- bundle$coords
shams   <- bundle$sham_shapes_df %>% mutate(source_NerveID = as.character(source_NerveID))
groups  <- bundle$nerve_groups %>% mutate(NerveID = as.character(NerveID), patient_id = as.character(patient_id))
pnd     <- bundle$per_nerve_dist_all %>% mutate(NerveID = as.character(NerveID), bin = bin_mid(dist_um))
psd     <- bundle$per_shape_dist_all %>% mutate(bin = bin_mid(dist_um))
um_per_pixel <- bundle$um_per_pixel

if (!exists("spe") || !inherits(spe, "SpatialExperiment")) spe <- readRDS(file.path(INPUT_DIR, "Integrated_NMF.rds"))
# Attach the spot annotation BY BARCODE, checked against the object's own
# identifiers before anything is overwritten (spatial coordinates belong to
# the object's columns, so a row-order mismatch would silently move labels).
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
# Each nerve must lie in exactly one section (distances, shams and nerve
# spots are all taken within the nerve's own section).
nerve_sections <- tapply(as.character(spe$sample_id[!is.na(spe$NerveID)]), as.character(spe$NerveID[!is.na(spe$NerveID)]),
                         function(x) length(unique(x)))
if (any(nerve_sections > 1)) stop("NerveIDs found in more than one section: ",
                                  paste(names(nerve_sections)[nerve_sections > 1], collapse = ", "))
stopifnot(identical(as.character(spe$barcode), as.character(meta$barcode)), all(GENES %in% rownames(spe)))
expr <- as.matrix(assay(spe, EXPR_ASSAY)[GENES, , drop = FALSE]); colnames(expr) <- spe$barcode

nerves <- groups %>% filter(Group == GROUP)
cat("Nerves in", GROUP, "per patient:\n"); print(table(nerves$patient_id))

# ---- 2. Per-nerve curves (real and own-control), per gene ------------------
nerve_curves <- function(g) {
  sc <- expr[g, ]
  nerve_spots <- data.frame(NerveID = as.character(meta$NerveID), barcode = meta$barcode) %>%
    filter(NerveID %in% nerves$NerveID)
  real <- bind_rows(
    nerve_spots %>% mutate(bin = 0),
    pnd %>% filter(NerveID %in% nerves$NerveID, !is.na(bin)) %>% select(NerveID, barcode, bin)) %>%
    mutate(score = sc[barcode]) %>%
    group_by(NerveID, bin) %>% summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(source = "Real")
  sh <- shams %>% filter(source_NerveID %in% nerves$NerveID)
  sham <- bind_rows(
    sh %>% transmute(ShamNerveID, barcode, bin = 0),
    psd %>% filter(ShamNerveID %in% sh$ShamNerveID, !is.na(bin)) %>% select(ShamNerveID, barcode, bin)) %>%
    mutate(score = sc[barcode]) %>%
    group_by(ShamNerveID, bin) %>% summarise(m = mean(score, na.rm = TRUE), .groups = "drop") %>%
    left_join(distinct(sh, ShamNerveID, NerveID = source_NerveID), by = "ShamNerveID") %>%
    group_by(NerveID, bin) %>% summarise(mean_score = mean(m, na.rm = TRUE), .groups = "drop") %>%
    mutate(source = "Sham")
  paired <- intersect(real$NerveID, sham$NerveID)
  bind_rows(real, sham) %>% filter(NerveID %in% paired) %>% mutate(gene = g) %>%
    left_join(nerves %>% select(NerveID, patient_id), by = "NerveID")
}
curves_nerve <- bind_rows(lapply(GENES, nerve_curves))

curves_patient <- curves_nerve %>% group_by(gene, source, patient_id, bin) %>%
  summarise(mean_score = mean(mean_score, na.rm = TRUE), n_nerves = n(), .groups = "drop")
curves_avg <- curves_patient %>% group_by(gene, source, bin) %>%
  summarise(mean = mean(mean_score, na.rm = TRUE),
            se   = sd(mean_score, na.rm = TRUE) / sqrt(sum(is.finite(mean_score))),
            n_patients = sum(is.finite(mean_score)), n_nerves = sum(n_nerves), .groups = "drop")
stopifnot(setequal(unique(curves_avg$bin), c(0, bin_breaks[-1] - BIN_WIDTH_UM / 2)))   # same bins as plasma figure
write.csv(curves_avg, file.path(TABLES_DIR, "fig3A_decay_curves_patient_level.csv"), row.names = FALSE)

# ---- 3. Lee's L of the example genes (numbers for the text, from step1) ------
# Taken from step1's final results table (no recomputation): per patient and
# test (nerve, interface), the real L and delta L = real - mean of the 200
# relocated masks (what panel C plots), and in how many patients the real L
# exceeds the relocated 95th percentile.
lee_tab <- read.csv(file.path(TABLES_DIR, "lee_sham_DE_genes_scores.csv")) %>% filter(Gene %in% GENES)
missing_genes <- setdiff(GENES, lee_tab$Gene)
if (length(missing_genes)) warning("Not in step1 results table: ", paste(missing_genes, collapse = ", "))
lee_real <- bind_rows(lapply(c("nerve", "interface"), function(tst) bind_rows(lapply(PATIENT_LABELS, function(pl)
  data.frame(gene = lee_tab$Gene, test = tst, patient = pl,
             L_real  = lee_tab[[paste0(tst, "_L_real_", pl)]],
             delta_L = lee_tab[[paste0(tst, "_delta_L_", pl)]],
             n_patients_above_q95 = lee_tab[[paste0(tst, "_n_patients_exceed_q95")]])))))
write.csv(lee_real, file.path(TABLES_DIR, "fig3A_lee_real_vs_sham_example_genes.csv"), row.names = FALSE)
cat("\nLee's L, real vs relocated masks (example genes):\n"); print(lee_real, digits = 3)

# ---- 4. Example nerve: typical profile over all genes -----------------------
prof <- curves_nerve %>% select(gene, NerveID, source, bin, mean_score) %>%
  pivot_wider(names_from = source, values_from = mean_score) %>% mutate(diff = Real - Sham)
n_bins <- length(bin_breaks)                                       # bin 0 + curve bins
complete <- prof %>% filter(!is.na(diff)) %>% count(gene, NerveID) %>% filter(n == n_bins) %>%
  count(NerveID) %>% filter(n == length(GENES)) %>% pull(NerveID)
pr <- prof %>% filter(NerveID %in% complete)
example_rank <- pr %>% group_by(gene, bin) %>% mutate(med = median(diff)) %>%
  group_by(gene) %>% mutate(scale = sd(diff)) %>%
  group_by(NerveID) %>% summarise(distance = sum(((diff - med) / scale)^2), .groups = "drop") %>%
  arrange(distance)
EXAMPLE_NERVE <- example_rank$NerveID[1]
ex_info <- nerves %>% filter(NerveID == EXAMPLE_NERVE)
cat("\nExample nerve:", EXAMPLE_NERVE, "| patient", ex_info$patient_id, "| spots",
    sum(meta$NerveID %in% EXAMPLE_NERVE), "| rank 1 of", nrow(example_rank), "complete invaded nerves\n")

# ---- 5. Plot helpers (true-scale spots, as in the plasma supplement) -------
lab <- function(p) ifelse(p %in% names(PATIENT_LABELS), PATIENT_LABELS[p], p)
hull <- function(xy) { h <- chull(xy[, 1], xy[, 2]); data.frame(x = xy[h, 1], y = xy[h, 2]) }
crop_df <- function(anchor_xy, smp, g) {
  m <- CROP_MARGIN_UM / um_per_pixel
  ix <- which(meta$sample_id == smp &
                coords[, 1] >= min(anchor_xy[, 1]) - m & coords[, 1] <= max(anchor_xy[, 1]) + m &
                coords[, 2] >= min(anchor_xy[, 2]) - m & coords[, 2] <= max(anchor_xy[, 2]) + m)
  data.frame(x = coords[ix, 1], y = coords[ix, 2], score = expr[g, meta$barcode[ix]],
             r = (SPOT_DIAMETER_UM / um_per_pixel) / 2)
}
crop_panel <- function(anchor_xy, smp, g, lims, title) {
  d <- crop_df(anchor_xy, smp, g); xr <- range(d$x)
  bar <- 100 / um_per_pixel
  ggplot(d) + ggforce::geom_circle(aes(x0 = x, y0 = y, r = r, fill = score), colour = NA) +
    scale_fill_viridis_c(option = "inferno", limits = lims, breaks = scales::breaks_pretty(n = 5),
                         name = paste0(g, "\n(", EXPR_ASSAY, ")")) +
    geom_polygon(data = hull(anchor_xy), aes(x, y), fill = NA, colour = NERVE_OUTLINE_COL, linewidth = 0.8) +
    annotate("segment", x = xr[1], xend = xr[1] + bar, y = min(d$y) - 2 * d$r[1], yend = min(d$y) - 2 * d$r[1],
             linewidth = 1.2) +
    annotate("text", x = xr[1] + bar / 2, y = min(d$y) - 4.5 * d$r[1], label = "100 µm", size = 3) +
    coord_fixed(clip = "off") + labs(title = title) +
    guides(fill = guide_colourbar(title.position = "left", title.vjust = 1, barwidth = unit(2.8, "cm"),
                                  barheight = unit(0.3, "cm"))) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = "bottom",
          legend.title = element_text(size = 9), legend.text = element_text(size = 8))
}

ex_smp <- unique(meta$sample_id[meta$NerveID %in% EXAMPLE_NERVE])
ex_xy  <- coords[which(meta$sample_id == ex_smp & meta$NerveID %in% EXAMPLE_NERVE), , drop = FALSE]
ex_sh  <- shams %>% filter(source_NerveID == EXAMPLE_NERVE)
far_rep <- ex_sh %>% group_by(ShamNerveID) %>%
  summarise(d = sqrt(sum((colMeans(coords[match(barcode, meta$barcode), , drop = FALSE]) - colMeans(ex_xy))^2)),
            .groups = "drop") %>% arrange(desc(d)) %>% slice(1) %>% pull(ShamNerveID)
sh_xy <- coords[match(ex_sh$barcode[ex_sh$ShamNerveID == far_rep], meta$barcode), , drop = FALSE]

# ---- 6. Figure: one row per gene (decay curve | zoom on the example nerve) --
pal <- c(Real = "#1B4F8C", Sham = "grey50")            # identical to the plasma-cell curve figure
CURVE_LABELS <- c(Real = "Real", Sham = "Control")
# Same y-axis range and tick spacing for all example genes, so the curves are comparable
Y_LIMS   <- range(c(curves_avg$mean - curves_avg$se, curves_avg$mean + curves_avg$se), na.rm = TRUE)
Y_BREAKS <- scales::breaks_pretty(n = 6)(Y_LIMS)
# Same colour scale (0 to the highest value) for the zooms of all example genes,
# so the spatial panels are comparable between genes
ZOOM_LIMS <- c(0, max(unlist(lapply(GENES, function(g)
  c(crop_df(ex_xy, ex_smp, g)$score, crop_df(sh_xy, ex_smp, g)$score))), na.rm = TRUE))
cat("Zoom colour scale (all genes):", round(ZOOM_LIMS, 2), "\n")

curve_panel <- function(g) {
  d <- curves_avg %>% filter(gene == g)
  ggplot(d, aes(bin, mean, colour = source, fill = source, group = source)) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.4, alpha = 0.9) +
    scale_colour_manual(values = pal, labels = CURVE_LABELS, name = NULL) +
    scale_fill_manual(values = pal, labels = CURVE_LABELS, name = NULL) +
    scale_y_continuous(limits = Y_LIMS, breaks = Y_BREAKS) +
    labs(title = g, x = "Distance from nerve (µm, binned)",
         y = paste0("Mean ", EXPR_ASSAY, "\n(\u00b1 SE across patients)")) +
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          plot.title = element_text(face = "italic", hjust = 0.5))
}
# Shared y-axis title (one label for both curves) when patchwork supports it
collect_axes <- packageVersion("patchwork") >= "1.2.0"
build_figure <- function(with_relocated) {
  rows <- lapply(GENES, function(g) {
    p_c <- curve_panel(g)
    z <- crop_panel(ex_xy, ex_smp, g, ZOOM_LIMS, "InvN")
    if (with_relocated) z <- z | crop_panel(sh_xy, ex_smp, g, ZOOM_LIMS, "Relocated shape")
    list(p_c, wrap_elements(full = z))     # wrapped: the zoom keeps its own colour bar
  })
  w <- if (with_relocated) c(1.5, 2) else c(1.5, 1)
  wrap_plots(unlist(rows, recursive = FALSE), ncol = 2, widths = w) +
    plot_layout(guides = "collect", axis_titles = if (collect_axes) "collect" else "keep") &
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 8))
}
caption_txt <- paste0("Example nerve ", EXAMPLE_NERVE, " (patient ", lab(ex_info$patient_id),
                      "), the invaded nerve with the most typical profile across genes;\nred outline = nerve; ",
                      "zoom shows spots up to ", CROP_MARGIN_UM, " µm from the nerve.")

fig <- build_figure(FALSE) + plot_annotation(caption = caption_txt)
out <- file.path(FIGURES_DIR, paste0("fig3A_decay_", paste(GENES, collapse = "_"), ".pdf"))
ggsave(out, fig, width = 7.5, height = 4 * length(GENES), device = pdf_device)
cat("\nSaved:", out, "\n")

if (MAKE_RELOCATED_ZOOM) {
  fig_s <- build_figure(TRUE) +
    plot_annotation(caption = paste0(caption_txt, "\nRelocated shape = its farthest relocated replicate."))
  out_s <- file.path(FIGURES_DIR, paste0("figS_decay_", paste(GENES, collapse = "_"), "_with_relocated_zoom.pdf"))
  ggsave(out_s, fig_s, width = 10.5, height = 4 * length(GENES), device = pdf_device)
  cat("Saved:", out_s, "\n")
}

