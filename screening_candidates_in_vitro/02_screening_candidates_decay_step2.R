# ======================================================================
# 06_screening_candidates_decay_step2.R
#
# Figure 4: distance-decay panels for the 15 screening candidates, one panel
# per gene, ordered by nerve association (delta L from 02 step1), each with
# the curves of the three nerve groups and a dotted line for the mean of
# tumour background far from any nerve.
#
# Same framework as Figure 3A and the plasma-cell decay figure:
#   * distances and bins from 00_decay_shared.R (150 um bins to 900 um,
#     bin 0 = the nerve spots; distances beyond 900 um are dropped, never
#     folded into the last bin), taken from the plasma bundle, whose
#     settings are checked by check_bundle_settings() -- the script stops
#     if the bundle was built with other settings;
#   * per nerve -> per patient (mean of nerves) -> mean +/- SE across
#     patients (patients are the replicates);
#   * control: geometry-matched relocated shapes of the same nerves within the
#     same section (plasma-analysis shams, made for EVERY nerve of all three
#     groups; replicates averaged per nerve, 1 nerve = 1 control). Version A
#     plots real - control; supplementary versions plot both curves.
# Background: Diseased spots that are neither nerve nor interface, with the
#   interface defined by the SAME contact rule as the spatial test (02 step1:
#   nneighbors = 2, measured 98-205 um; the measured range is printed). This
#   replaces the earlier nneighbors = 8 zone (184 coordinate units, 49-
#   neighbour cap), which was never calibrated and differed from 02 step1.
# No caches: everything is computed from the inputs on each run.
#
# INPUTS: input/Integrated_NMF.rds, input/DE_SPE_annotaion.csv,
#   processed/figure_bundle_for_real_vs_relocated.rds (plasma step1),
#   output/06_candidate_gene_screening/tables/screening_candidates_summary.csv (06 step1)
# OUTPUTS (output/06_candidate_gene_screening/):
#   figures/figure4_screening_candidates_excess_over_relocated.pdf / .png (version A:
#     three nerve groups, real minus own relocated control, line at 0)
#   figures/figure4_screening_candidates_by_deltaL.pdf / .png (version B: three nerve
#     groups, real curves + Diseased background line, as the original figure)
#   if MAKE_RELOCATED_VERSIONS: figS_..._real_vs_relocated_all_groups / _invaded
#     (each real curve with its relocated-shape control, absolute levels)
#   tables/figure4_screening_candidates_curves_patient_level.csv,
#   tables/figure4_screening_candidates_deltaL_and_background.csv
# ======================================================================

library(SpatialExperiment)
library(RANN)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

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
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# Interface = the same contact rule as 02 step1 (Lee's L): a non-nerve spot is
# interface if a nerve spot is among its 12 nearest spots within 46 coordinate
# units (nneighbors = 2; measured 98-205 um in 02 step1). Background =
# Diseased spots that are neither nerve nor interface -- the same background
# population as the spatial test.
IF_NNEIGHBORS  <- 2
IF_RADIUS      <- IF_NNEIGHBORS * 23
IF_K           <- IF_NNEIGHBORS * 6 + 1
N_COLS         <- 5
MAKE_RELOCATED_VERSIONS <- TRUE   # also save supplementary versions with the relocated-shape controls
GROUP_COLORS   <- c("Invaded nerves (Diseased)"     = "#F8766D",
                    "Non-invaded nerves (Diseased)" = "#00BA38",
                    "Non-invaded nerves (Adjacent)" = "#00BFC4")

pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")

# 00_decay_shared.R is the shared helper of spatial_enrichment_prioritisation;
# place a copy in scripts/ or input/ so the zones and bins are identical.
SHARED_CANDIDATES <- file.path(c(file.path(PROJECT_DIR, "scripts"), INPUT_DIR, PROJECT_DIR),
                               "00_decay_shared.R")
SHARED_PATH <- SHARED_CANDIDATES[file.exists(SHARED_CANDIDATES)][1]
if (is.na(SHARED_PATH)) stop("00_decay_shared.R not found.")
source(SHARED_PATH)

# ---- 1. Data -------------------------------------------------------------
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

# Written by 01_plasma_decay_analysis_step1.Rmd
# (spatial_enrichment_prioritisation); place a copy in input/ before running.
bundle <- readRDS(file.path(INPUT_DIR, "figure_bundle_for_real_vs_relocated.rds"))
check_bundle_settings(bundle)                      # stops if bins / distances were built differently
stopifnot(identical(as.character(spe$barcode), as.character(bundle$spe_meta$barcode)))
um_per_pixel <- bundle$um_per_pixel
# Each nerve must lie in exactly one section (distances, shams and nerve
# spots are all taken within the nerve's own section).
nerve_sections <- tapply(as.character(spe$sample_id[!is.na(spe$NerveID)]), as.character(spe$NerveID[!is.na(spe$NerveID)]),
                         function(x) length(unique(x)))
if (any(nerve_sections > 1)) stop("NerveIDs found in more than one section: ",
                                  paste(names(nerve_sections)[nerve_sections > 1], collapse = ", "))
groups <- bundle$nerve_groups %>% mutate(NerveID = as.character(NerveID), patient_id = as.character(patient_id),
                                         Group = as.character(Group))
pnd <- bundle$per_nerve_dist_all %>% mutate(NerveID = as.character(NerveID), bin = bin_mid(dist_um)) %>%
  filter(!is.na(bin))
shams <- bundle$sham_shapes_df %>% mutate(source_NerveID = as.character(source_NerveID))
psd <- bundle$per_shape_dist_all %>% mutate(bin = bin_mid(dist_um)) %>% filter(!is.na(bin))

# ---- 2. Candidates, ordered by nerve association -----------------------------
summ <- read.csv(file.path(TABLES_DIR, "screening_candidates_summary.csv"))
cand <- summ %>% select(Gene, Nominated_as, Nerve_deltaL, Nerve_patients_above_q95) %>%
  filter(Gene %in% rownames(spe)) %>% arrange(desc(Nerve_deltaL))
dropped <- setdiff(summ$Gene, cand$Gene)
if (length(dropped)) cat("Not in the expression data (no panel):", paste(dropped, collapse = ", "), "\n")

# ---- 3. Interface (02 step1 contact rule) and background -------------------
coords <- spatialCoords(spe)
is_nerve <- !is.na(spe$NerveID) & spe$NerveID != ""
is_interface <- rep(FALSE, ncol(spe)); if_dist_um <- numeric(0)
for (s in unique(spe$sample_id)) {
  idx <- which(spe$sample_id == s)
  nn <- nn2(coords[idx, , drop = FALSE], coords[idx, , drop = FALSE],
            k = min(IF_K, length(idx)), searchtype = "radius", radius = IF_RADIUS)
  nb <- nn$nn.idx[, -1, drop = FALSE]                        # neighbours (0 = none within radius)
  touches <- apply(nb, 1, function(r) any(is_nerve[idx[r[r != 0]]]))
  ring <- idx[touches & !is_nerve[idx]]
  is_interface[ring] <- TRUE
  nv <- idx[is_nerve[idx]]
  if (length(ring) && length(nv)) if_dist_um <- c(if_dist_um, nn2(coords[nv, , drop = FALSE], coords[ring, , drop = FALSE], k = 1)$nn.dists[, 1] * um_per_pixel)
}
cat("Interface spots:", sum(is_interface), "| distance to nearest nerve spot:",
    paste(round(range(if_dist_um)), collapse = "-"), "um (02 step1 reported 98-205 um)\n")
bg_idx <- which(spe$context %in% "Diseased" & !is_nerve & !is_interface)
cat("Background spots (Diseased, not nerve, not interface):", length(bg_idx), "\n")

# ---- 4. Curves: per nerve -> per patient -> mean +/- SE across patients ------
# Real: around each nerve. Relocated: around the geometry-matched relocated
# shapes of the same nerve (plasma-analysis shams), replicates averaged per
# nerve, so 1 nerve = 1 control curve. Only nerves with both are used.
expr <- as.matrix(assay(spe, "logcounts")[cand$Gene, , drop = FALSE]); colnames(expr) <- spe$barcode
nerve_spots <- data.frame(NerveID = as.character(spe$NerveID[is_nerve]), barcode = spe$barcode[is_nerve]) %>%
  filter(NerveID %in% groups$NerveID)
real_bins <- bind_rows(nerve_spots %>% mutate(bin = 0),
                       pnd %>% filter(NerveID %in% groups$NerveID) %>% select(NerveID, barcode, bin))
sham_bins <- bind_rows(shams %>% filter(source_NerveID %in% groups$NerveID) %>% transmute(ShamNerveID, barcode, bin = 0),
                       psd %>% filter(ShamNerveID %in% shams$ShamNerveID) %>% select(ShamNerveID, barcode, bin)) %>%
  left_join(distinct(shams, ShamNerveID, NerveID = source_NerveID), by = "ShamNerveID")
paired <- intersect(unique(real_bins$NerveID), unique(sham_bins$NerveID))
curves <- bind_rows(lapply(cand$Gene, function(g) {
  real <- real_bins %>% filter(NerveID %in% paired) %>% mutate(x = expr[g, barcode]) %>%
    group_by(NerveID, bin) %>% summarise(m = mean(x), .groups = "drop") %>% mutate(Source = "Real")
  sham <- sham_bins %>% filter(NerveID %in% paired) %>% mutate(x = expr[g, barcode]) %>%
    group_by(ShamNerveID, NerveID, bin) %>% summarise(m = mean(x), .groups = "drop") %>%
    group_by(NerveID, bin) %>% summarise(m = mean(m), .groups = "drop") %>% mutate(Source = "Relocated")
  # Excess over the nerve's own control: real - relocated, per nerve and bin
  excess <- inner_join(real %>% select(NerveID, bin, r = m), sham %>% select(NerveID, bin, c = m),
                       by = c("NerveID", "bin")) %>% transmute(NerveID, bin, m = r - c, Source = "Excess")
  bind_rows(real, sham, excess) %>%
    left_join(groups %>% select(NerveID, Group, patient_id), by = "NerveID") %>%
    group_by(Source, Group, patient_id, bin) %>% summarise(m = mean(m), n_nerves = n(), .groups = "drop") %>%
    group_by(Source, Group, bin) %>%
    summarise(mean = mean(m), se = sd(m) / sqrt(n()), n_patients = n(), n_nerves = sum(n_nerves), .groups = "drop") %>%
    mutate(Gene = g)
}))
cat("Nerves with real and relocated curves:", length(paired), "\n")
write.csv(curves, file.path(TABLES_DIR, "figure4_screening_candidates_curves_patient_level.csv"), row.names = FALSE)
bg_mean <- setNames(rowMeans(expr[, bg_idx, drop = FALSE]), cand$Gene)
write.csv(cand %>% mutate(Background_mean = bg_mean[Gene]),
          file.path(TABLES_DIR, "figure4_screening_candidates_deltaL_and_background.csv"), row.names = FALSE)

# ---- 5. Panels -------------------------------------------------------------
# Figures are drawn at their final print size (FIG_W_PT x FIG_H_PT points =
# Illustrator units; 72 pt = 1 inch), with text sizes in pt as printed, and
# without in-plot captions (method text goes in the figure legend). Text is
# standard Helvetica (pdf_device), editable in Illustrator.
# This section can be run on its own after section 0 (settings): if the curves
# are not in memory, it reads the two tables written in section 4.
FIG_W_PT  <- 589          # width  (pt) -- all figures of this script
FIG_H_PT  <- 620          # height (pt)
FIG_FONT  <- 8            # pt: axis titles, legend, panel titles
FIG_SMALL <- 7            # pt: tick labels, panel subtitles
FIG_ASPECT <- 1.2         # panel height / width (1 = square; 1.5 fills the full 589 x 620 pt height with 5 x 3)
FIG_N_COLS <- 5           # panels per row (15 genes: 5 x 3);
# if the last row is incomplete, the legend fills the empty cell
# Same y-axis in every panel, so candidates are directly comparable
EXCESS_YLIM   <- c(-0.5, 1.5)                 # real - relocated (logcounts)
EXCESS_BREAKS <- c(-0.5, 0, 0.5, 1.0, 1.5)
if (!exists("curves") || !exists("cand") || !exists("bg_mean")) {
  curves <- read.csv(file.path(TABLES_DIR, "figure4_screening_candidates_curves_patient_level.csv"))
  cand   <- read.csv(file.path(TABLES_DIR, "figure4_screening_candidates_deltaL_and_background.csv"))
  bg_mean <- setNames(cand$Background_mean, cand$Gene)
  cat("Figures: read curves for", length(unique(curves$Gene)), "genes from the section-4 tables.\n")
}
# Real-curve figures (absolute logcounts): one common range for all genes,
# from 0 to the highest curve (+SE) or background line of any gene
REAL_YLIM   <- c(0, max(c(curves$mean[curves$Source != "Excess"] + curves$se[curves$Source != "Excess"], bg_mean),
                        na.rm = TRUE))
REAL_BREAKS <- scales::breaks_pretty(n = 4)(REAL_YLIM)
cat("Common y-axis: excess", paste(EXCESS_YLIM, collapse = " to "), "| real curves 0 to", round(REAL_YLIM[2], 2), "\n")
GROUP_SHORT <- c("Invaded nerves (Diseased)" = "InvN", "Non-invaded nerves (Diseased)" = "NI (tumour-bearing)",
                 "Non-invaded nerves (Adjacent)" = "NI (adjacent)")
theme_fig <- function() theme_classic(base_size = FIG_FONT, base_family = "Helvetica") +
  theme(plot.title = element_text(face = "bold.italic", size = FIG_FONT, margin = margin(0, 0, 1, 0)),
        plot.subtitle = element_text(size = FIG_SMALL, margin = margin(0, 0, 2, 0)),
        axis.text = element_text(size = FIG_SMALL, colour = "black"), axis.title = element_text(size = FIG_FONT),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.text = element_text(size = FIG_FONT), legend.key.width = unit(5, "mm"),
        legend.position = "bottom", plot.margin = margin(3, 3, 3, 3), aspect.ratio = FIG_ASPECT)
panel_title <- function(g) {
  r <- cand[cand$Gene == g, ]
  list(title = if (r$Gene != r$Nominated_as) paste0(g, " (as ", r$Nominated_as, ")") else g,
       sub = if (is.na(r$Nerve_deltaL)) "not in the spatial test" else
         sprintf("dL = %.3f (%d/3 patients)", r$Nerve_deltaL, as.integer(r$Nerve_patients_above_q95)))
}
x_scale <- scale_x_continuous(breaks = c(0, 300, 600, 900), limits = c(0, 900))
panel <- function(g, groups_shown = names(GROUP_COLORS), cols = GROUP_COLORS, sources = c("Real", "Relocated")) {
  tt <- panel_title(g)
  d <- filter(curves, Gene == g, Group %in% groups_shown, Source %in% sources) %>%
    mutate(Source = factor(Source, c("Real", "Relocated")))
  ggplot(d, aes(bin, mean, colour = Group, fill = Group, linetype = Source, group = interaction(Group, Source))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = bg_mean[[g]], linetype = "dotted", colour = "black", linewidth = 0.4) +
    geom_ribbon(data = filter(d, Source == "Real"), aes(ymin = mean - se, ymax = mean + se), alpha = 0.15, colour = NA) +
    geom_line(aes(linewidth = Source)) + geom_point(data = filter(d, Source == "Real"), size = 0.7) +
    scale_colour_manual(values = cols, breaks = groups_shown, labels = GROUP_SHORT[groups_shown], name = NULL) +
    scale_fill_manual(values = cols, breaks = groups_shown, labels = GROUP_SHORT[groups_shown], name = NULL) +
    scale_linetype_manual(values = c(Real = "solid", Relocated = "22"), name = NULL,
                          labels = c(Real = "Real nerves", Relocated = "Relocated shapes"),
                          guide = if (length(sources) > 1) "legend" else "none") +
    scale_linewidth_manual(values = c(Real = 0.6, Relocated = 0.4), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = REAL_BREAKS) + x_scale +
    coord_cartesian(ylim = REAL_YLIM) +
    labs(title = tt$title, subtitle = tt$sub, x = "Distance from nerve (\u00b5m)", y = "Mean logcounts") +
    theme_fig()
}
# Excess over own relocated control (three nerve groups, reference line at 0)
panel_excess <- function(g) {
  tt <- panel_title(g)
  d <- filter(curves, Gene == g, Source == "Excess")
  ggplot(d, aes(bin, mean, colour = Group, fill = Group)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.6) + geom_point(size = 0.7) +
    scale_colour_manual(values = GROUP_COLORS, breaks = names(GROUP_COLORS), labels = GROUP_SHORT, name = NULL) +
    scale_fill_manual(values = GROUP_COLORS, breaks = names(GROUP_COLORS), labels = GROUP_SHORT, name = NULL) +
    scale_y_continuous(breaks = EXCESS_BREAKS) + x_scale +
    coord_cartesian(ylim = EXCESS_YLIM) +
    labs(title = tt$title, subtitle = tt$sub, x = "Distance from nerve (\u00b5m)",
         y = "Real - relocated (logcounts)") +
    theme_fig()
}
# Shared legend; axis titles only on the outer panels
assemble <- function(panels) {
  # axis titles only on the outer panels: y title in the first column, x title
  # on the lowest panel of each column
  n <- length(panels); n_rows <- ceiling(n / FIG_N_COLS)
  for (i in seq_len(n)) {
    col <- (i - 1) %% FIG_N_COLS + 1; row <- (i - 1) %/% FIG_N_COLS + 1
    lowest_in_col <- (i + FIG_N_COLS) > n
    if (col > 1) panels[[i]] <- panels[[i]] + labs(y = NULL)
    if (!lowest_in_col) panels[[i]] <- panels[[i]] + labs(x = NULL)
  }
  n_cells <- FIG_N_COLS * ceiling(length(panels) / FIG_N_COLS)
  legend_in_grid <- n_cells > length(panels)          # an empty cell is left: put the legend there
  if (legend_in_grid) panels <- c(panels, list(guide_area()))
  fig <- wrap_plots(panels, ncol = FIG_N_COLS) + plot_layout(guides = "collect")
  fig & theme(legend.position = if (legend_in_grid) "right" else "bottom",
              legend.direction = if (legend_in_grid) "vertical" else "horizontal",
              legend.box = if (legend_in_grid) "vertical" else "horizontal")
}
save_both <- function(fig, file_stub) {
  out <- file.path(FIGURES_DIR, file_stub)
  ggsave(paste0(out, ".pdf"), fig, width = FIG_W_PT / 72, height = FIG_H_PT / 72, device = pdf_device)
  ggsave(paste0(out, ".png"), fig, width = FIG_W_PT / 72, height = FIG_H_PT / 72, dpi = 300)
  cat("Saved:", file_stub, "(", FIG_W_PT, "x", FIG_H_PT, "pt )\n")
}
# Legend text (not drawn in the figures):
#  excess: real - relocated per nerve and bin (geometry-matched relocated shapes
#    of the same nerve within the same section), mean per patient, mean +/- SE
#    across patients; 0 = no different from the same shape elsewhere in the
#    tissue. Bin 0 = nerve spots; 150 um bins to 900 um. Panels ordered by nerve
#    delta L (02 step1).
#  real curves: dotted line = mean of Diseased spots that are neither nerve nor
#    in direct contact with a nerve (same background as the spatial test).
save_both(assemble(lapply(cand$Gene, panel_excess)), "figure4_screening_candidates_excess_over_relocated")
# Figure 4 (main): the three nerve groups, real curves (as in the original figure)
save_both(assemble(lapply(cand$Gene, panel, sources = "Real")), "figure4_screening_candidates_by_deltaL")
# Optional supplementary versions with the relocated-shape controls
if (MAKE_RELOCATED_VERSIONS) {
  save_both(assemble(lapply(cand$Gene, panel, sources = c("Real", "Relocated"))),
            "figS_screening_candidates_real_vs_relocated_all_groups")
  save_both(assemble(lapply(cand$Gene, panel, groups_shown = "Invaded nerves (Diseased)", sources = c("Real", "Relocated"))),
            "figS_screening_candidates_real_vs_relocated_invaded")
}

cat("Saved Figure 4 and tables in", file.path(PROJECT_DIR, "output", ANALYSIS_NAME), "\n")

