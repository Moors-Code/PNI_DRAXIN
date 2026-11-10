# ======================================================================
# Zoom comparison: observed nerve vs. its own relocated geometry (same
# shape and size, placed in the same tissue section, so same patient), one
# example per nerve category, with true-to-scale spot rendering and scale
# bars, using the CALIBRATED um_per_pixel. The three examples may come from
# different patients (see EXAMPLE_SELECTION below).
#
# STANDALONE: loads its required objects from the bundle saved by
# 01_plasma_decay_analysis_step1.Rmd (run that first).
#
# REVISION (review):
#   - Uses 00_decay_shared.R, so the near/far zones and curve bins are the
#     same as in step1/step3 (near <= 150 um, far 650-1000 um; 150 um bins
#     to 900 um). Previously candidates were scored on BINNED curves, where
#     "far" meant the 675 and 825 bins (= 600-900 um, with 750-1000 um
#     folded into 825) -- a different contrast from the analysis.
#   - Candidates are scored on the same paired quantities used in the
#     statistics (each nerve vs the mean of its own relocated replicates),
#     taken from bundle$paired_nerve; only nerves with a control qualify.
#   - Distance profile: relocated line = mean of the nerve's replicates
#     (the control used in the statistics); the map still shows one
#     replicate. Bins beyond 900 um are excluded.
#   - Purpose: show that relocation preserves the nerve geometry and how a
#     TYPICAL nerve compares with its own relocated controls.
#   - EXAMPLE_SELECTION switch (the default modes never look at the controls
#     or at the ordering of the categories):
#       "median_profile_size_matched" (default): as "median_profile", but
#         among the TYPICAL_FRACTION of nerves with the most typical profile
#         in each category (at least 3), the three examples are chosen to be
#         as similar in size (spot count) as possible, so the true-to-scale
#         maps are comparable. Size is unrelated to the outcome, so this
#         keeps the selection neutral.
#       "median_profile": per category, the nerve whose whole
#         distance profile (observed minus mean of its own relocated
#         replicates, per 150 um bin, 0-900 um) is closest to the category's
#         median profile (sum of squared differences across bins). Typical in
#         shape, not just in one number -- e.g. a nerve whose curve is
#         dominated by a neighbouring nerve's halo is not chosen.
#       "median_decay": per category, the nerve whose own decay
#         (near <= 150 um minus far 650-1000 um) is closest to the category
#         median. Rows may come from different patients.
#       "median_decay_same_patient": same rule, but all three rows from one
#         patient (the patient whose nerves are closest to the three medians).
#       "clearest_ordered": the original selection (clearest examples that
#         preserve NI-TB > InvN > NI-Adj). Comparison only; if used, the
#         legend must say the examples were selected for clarity.
#   - Every mode writes a representativeness table: each example's decay
#     and its percentile within its category.
#   - PDFs saved with cairo_pdf when available, so the "\u00b5" in labels prints.
#   - Refuses to run on a bundle built with different decay settings.
# ======================================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(ggforce)
library(SpatialExperiment)

count     <- dplyr::count
rename    <- dplyr::rename
summarise <- dplyr::summarise
summarize <- dplyr::summarize
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
filter    <- dplyr::filter
select    <- dplyr::select
slice     <- dplyr::slice
pull      <- dplyr::pull
distinct  <- dplyr::distinct
group_by  <- dplyr::group_by

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input"))
  stop("Run this script from the project folder (the one containing input/).")

PROJECT_DIR   <- "."
INPUT_DIR     <- file.path(PROJECT_DIR, "input")
PROCESSED_DIR <- file.path(PROJECT_DIR, "processed")

ANALYSIS_NAME <- "01_plasma_cells_supplement"
FIGURES_DIR   <- file.path(PROJECT_DIR, "output", ANALYSIS_NAME, "figures")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 0a. Shared code ----------------------------------------------------
# 00_decay_shared.R lives next to the scripts; the scripts/ folder is
# checked first, then the project folder itself.
SHARED_CANDIDATES <- file.path(c(file.path(PROJECT_DIR, "scripts"), file.path(PROJECT_DIR, "script"), PROJECT_DIR),
                               "00_decay_shared.R")
SHARED_PATH <- SHARED_CANDIDATES[file.exists(SHARED_CANDIDATES)][1]
if (is.na(SHARED_PATH)) stop("00_decay_shared.R not found. Looked in:\n  ",
                             paste(SHARED_CANDIDATES, collapse = "\n  "))
source(SHARED_PATH)
# cairo_pdf prints the micro sign; defined here too in case the shared file
# is an older version without pdf_device
if (!exists("pdf_device")) pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

# "median_decay" (recommended), "median_decay_same_patient", or
# "clearest_ordered" (original behaviour, for comparison only) -- see header
# "median_profile_size_matched" (recommended), "median_profile", "median_decay", "median_decay_same_patient",
# or "clearest_ordered" (original behaviour, comparison only) -- see header
EXAMPLE_SELECTION <- "median_profile_size_matched"
stopifnot(EXAMPLE_SELECTION %in% c("median_profile_size_matched", "median_profile", "median_decay",
                                   "median_decay_same_patient", "clearest_ordered"))
TYPICAL_FRACTION <- 0.25   # size-matched mode: pool = most typical quarter of each category (min 3 nerves)

# ---- 0b. Load the saved object bundle -----------------------------------
BUNDLE_PATH <- file.path(PROCESSED_DIR, "figure_bundle_for_real_vs_relocated.rds")
if (!file.exists(BUNDLE_PATH)) {
  stop("Bundle not found at: ", BUNDLE_PATH,
       "\nRun 01_plasma_decay_analysis_step1.Rmd first to generate it.")
}
cat("Loading object bundle from:", BUNDLE_PATH, "\n")
bundle <- readRDS(BUNDLE_PATH)
check_bundle_settings(bundle)

spe                 <- bundle$spe_meta
coords              <- bundle$coords
plasma_score        <- bundle$plasma_score
sham_shapes_df      <- bundle$sham_shapes_df
nerve_groups        <- bundle$nerve_groups
paired_nerve        <- bundle$paired_nerve
per_nerve_dist_all  <- bundle$per_nerve_dist_all
per_shape_dist_all  <- bundle$per_shape_dist_all
cat("Bundle loaded successfully.\n\n")

# ---- 0c. Local repair: rebuild spe metadata if the bundled copy is broken -
if (is.null(spe$NerveID) || !("NerveID" %in% colnames(spe))) {
  cat("Bundled spe metadata is missing NerveID -- repairing locally.\n")
  spe_full <- readRDS(file.path(INPUT_DIR, "Integrated_NMF.rds"))
  colData(spe_full) <- DataFrame(read.csv(file.path(INPUT_DIR, "DE_SPE_annotaion.csv")))
  colnames(spe_full) <- spe_full$barcode
  spe <- data.frame(barcode = spe_full$barcode, sample_id = spe_full$sample_id,
                    NerveID = spe_full$NerveID, stringsAsFactors = FALSE)
  rm(spe_full)
  cat("Repaired spe metadata:", nrow(spe), "spots,", sum(!is.na(spe$NerveID)), "nerve spots.\n\n")
}

# ---- 0d. Load CALIBRATED um_per_pixel -----------------------------------
UM_PER_PIXEL_PATH <- file.path(INPUT_DIR, "um_per_pixel_CALIBRATED.rds")
if (!file.exists(UM_PER_PIXEL_PATH)) {
  stop("Calibrated um_per_pixel not found at: ", UM_PER_PIXEL_PATH,
       "\nRun um_per_pixel_calibration.R first.")
}
um_per_pixel <- readRDS(UM_PER_PIXEL_PATH)
cat("Loaded calibrated um_per_pixel:", round(um_per_pixel, 6),
    "(equivalent PIXELS_PER_100UM:", round(100/um_per_pixel, 3), ")\n\n")

CROP_MARGIN_UM <- 130
N_RELOCATED_OVERLAY <- 5
SPOT_DIAMETER_UM <- 55

base_theme <- theme_void(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 11))

# ---- 1. Geometry helpers ----------------------------------------------------
get_bbox <- function(pt_coords, margin_px) {
  c(xmin = min(pt_coords[, 1]) - margin_px, xmax = max(pt_coords[, 1]) + margin_px,
    ymin = min(pt_coords[, 2]) - margin_px, ymax = max(pt_coords[, 2]) + margin_px)
}

crop_spots <- function(bbox, spe, coords, sample_id_filter) {
  which(spe$sample_id == sample_id_filter &
          coords[, 1] >= bbox["xmin"] & coords[, 1] <= bbox["xmax"] &
          coords[, 2] >= bbox["ymin"] & coords[, 2] <= bbox["ymax"])
}

convex_hull_df <- function(pt_coords) {
  h <- chull(pt_coords[, 1], pt_coords[, 2])
  data.frame(x = pt_coords[h, 1], y = pt_coords[h, 2])
}

# ---- Scale bar as a SEPARATE panel, stacked below the spatial plot -------
build_scale_bar_panel <- function(um_per_pixel, x_range, bar_um = 50, bar_color = "black", label_size = 3.2) {
  bar_length_px <- bar_um / um_per_pixel
  x_start <- x_range[1] + 0.05 * diff(x_range)
  x_end <- x_start + bar_length_px
  df <- data.frame(x = c(x_start, x_end), y = c(0, 0))
  ggplot(df, aes(x = x, y = y)) +
    geom_line(color = bar_color, linewidth = 1.4) +
    annotate("text", x = (x_start + x_end) / 2, y = -0.15, label = paste0(bar_um, " µm"),
             color = bar_color, size = label_size, vjust = 1) +
    coord_cartesian(xlim = x_range, clip = "off") +
    ylim(-0.5, 0.3) +
    theme_void()
}

# ---- Crop panel builder, TRUE-SCALE spots via geom_circle ----------------
build_crop_panel <- function(anchor_coords, spe, coords, plasma_score, sample_id_filter,
                             um_per_pixel, crop_margin_um, color_limits, title,
                             extra_outlines = NULL, spot_diameter_um = SPOT_DIAMETER_UM) {
  margin_px <- crop_margin_um / um_per_pixel
  bbox <- get_bbox(anchor_coords, margin_px)
  crop_idx <- crop_spots(bbox, spe, coords, sample_id_filter)
  crop_df <- data.frame(x = coords[crop_idx, 1], y = coords[crop_idx, 2],
                        score = plasma_score[spe$barcode[crop_idx]])
  crop_df$r <- (spot_diameter_um / um_per_pixel) / 2
  outline <- convex_hull_df(anchor_coords)
  
  p <- ggplot(crop_df) +
    ggforce::geom_circle(aes(x0 = x, y0 = y, r = r, fill = score), color = NA) +
    scale_fill_viridis_c(option = "inferno", limits = color_limits)
  
  if (!is.null(extra_outlines) && nrow(extra_outlines) > 0) {
    p <- p + geom_polygon(data = extra_outlines, aes(x = x, y = y, group = shape_id), inherit.aes = FALSE,
                          fill = NA, color = "white", alpha = 0.35, linewidth = 0.5, linetype = "dashed")
  }
  
  p + geom_polygon(data = outline, aes(x = x, y = y), inherit.aes = FALSE,
                   fill = NA, color = "white", linewidth = 1) +
    coord_fixed(xlim = range(crop_df$x), ylim = range(crop_df$y)) +
    labs(title = title) +
    base_theme
}

build_crop_panel_with_scalebar <- function(anchor_coords, spe, coords, plasma_score, sample_id_filter,
                                           um_per_pixel, crop_margin_um, color_limits, title,
                                           extra_outlines = NULL, spot_diameter_um = SPOT_DIAMETER_UM,
                                           bar_um = 50) {
  p_spatial <- build_crop_panel(anchor_coords, spe, coords, plasma_score, sample_id_filter,
                                um_per_pixel, crop_margin_um, color_limits, title,
                                extra_outlines, spot_diameter_um)
  # Scale-bar panel shares the spatial panel's x-range, so patchwork cannot
  # stretch the bar to a different physical length.
  margin_px <- crop_margin_um / um_per_pixel
  bbox <- get_bbox(anchor_coords, margin_px)
  crop_idx <- crop_spots(bbox, spe, coords, sample_id_filter)
  x_range <- range(coords[crop_idx, 1])
  
  p_scale <- build_scale_bar_panel(um_per_pixel, x_range = x_range, bar_um = bar_um)
  p_spatial / p_scale + plot_layout(heights = c(20, 1))
}

get_crop_score_range <- function(anchor_coords, spe, coords, plasma_score, sample_id_filter, um_per_pixel, crop_margin_um) {
  margin_px <- crop_margin_um / um_per_pixel
  bbox <- get_bbox(anchor_coords, margin_px)
  crop_idx <- crop_spots(bbox, spe, coords, sample_id_filter)
  range(plasma_score[spe$barcode[crop_idx]], na.rm = TRUE)
}

select_farthest_replicate <- function(nerve_id, sample_id, nerve_coords, sham_shapes_df, coords, spe, um_per_pixel) {
  this_nerve_shams <- sham_shapes_df %>% filter(source_NerveID == nerve_id, sample_id == !!sample_id)
  if (nrow(this_nerve_shams) == 0) stop("No sham shapes for NerveID=", nerve_id, ", sample_id=", sample_id)
  nerve_centroid <- colMeans(nerve_coords)
  replicate_ids <- unique(this_nerve_shams$ShamNerveID)
  replicate_dists <- sapply(replicate_ids, function(sid) {
    sc <- coords[match(this_nerve_shams$barcode[this_nerve_shams$ShamNerveID == sid], spe$barcode), , drop = FALSE]
    sqrt(sum((colMeans(sc) - nerve_centroid)^2)) * um_per_pixel
  })
  list(shape_id = replicate_ids[which.max(replicate_dists)], separation_um = max(replicate_dists),
       this_nerve_shams = this_nerve_shams)
}

# ---- 2. Example selection ------------------------------------------------
# Selection metric (default modes): each nerve's own decay effect size,
#   decay_magnitude = mean proportion near (<= 150 um) - far (650-1000 um),
# from bundle$paired_nerve (shared decay function; only nerves that have a
# relocated control qualify). The example per category is the nerve whose
# decay is closest to that category's MEDIAN -- a typical nerve, robust to
# the few strongly enriched nerves that pull the mean up. The selection never
# looks at the relocated controls or at the ordering of categories.
cat_labels <- c(adj = "Non-invaded nerves (Adjacent)", ni = "Non-invaded nerves (Diseased)",
                invn = "Invaded nerves (Diseased)")
cat_short  <- c(adj = "NI-Adj", ni = "NI-TB", invn = "InvN")

candidate_table <- function(group_label) {
  paired_nerve %>%
    filter(Group == group_label) %>%
    mutate(NerveID = as.character(NerveID)) %>%
    left_join(nerve_groups %>% select(NerveID, sample_id), by = "NerveID") %>%
    transmute(NerveID, sample_id, patient_id = as.character(patient_id),
              decay = decay_magnitude,                 # selection metric
              delta_decay,                             # sham-corrected (reported only)
              near_sep = delta_near, far_div = delta_far)
}
cats <- lapply(cat_labels, candidate_table)
if (any(sapply(cats, nrow) == 0)) stop("A nerve category has no paired nerves.")

closest_to_median <- function(df) df[which.min(abs(df$decay - median(df$decay))), ]

# Per-nerve profiles: observed minus mean of own relocated replicates, per bin
nerve_profiles <- function() {
  real_bins <- per_nerve_dist_all %>% filter(!is.na(bin)) %>%
    group_by(NerveID = as.character(NerveID), bin) %>%
    summarise(real = mean(score, na.rm = TRUE), .groups = "drop")
  sham_bins <- per_shape_dist_all %>% filter(!is.na(bin)) %>%
    group_by(ShamNerveID, bin) %>% summarise(m = mean(score, na.rm = TRUE), .groups = "drop") %>%
    left_join(distinct(sham_shapes_df, ShamNerveID, source_NerveID), by = "ShamNerveID") %>%
    group_by(NerveID = as.character(source_NerveID), bin) %>%
    summarise(sham = mean(m, na.rm = TRUE), .groups = "drop")
  inner_join(real_bins, sham_bins, by = c("NerveID", "bin")) %>% mutate(diff = real - sham)
}

# Nerve size = number of spots in the nerve (in its own section)
nerve_size <- function(nid, smp) sum(spe$sample_id == smp & !is.na(spe$NerveID) & spe$NerveID == nid)

# Profile distance of every nerve (with all bins) to its category's median profile
profile_ranked <- function(df, prof) {
  n_bins <- length(bin_breaks) - 1
  pr <- prof %>% filter(NerveID %in% df$NerveID)
  complete_ids <- pr %>% count(NerveID) %>% filter(n == n_bins) %>% pull(NerveID)   # all bins present
  pr <- pr %>% filter(NerveID %in% complete_ids)
  med <- pr %>% group_by(bin) %>% summarise(med = median(diff), .groups = "drop")
  dist <- pr %>% left_join(med, by = "bin") %>% group_by(NerveID) %>%
    summarise(profile_distance = sum((diff - med)^2), .groups = "drop")
  df %>% inner_join(dist, by = "NerveID") %>% arrange(profile_distance) %>%
    mutate(profile_rank = row_number(),
           n_spots = mapply(nerve_size, NerveID, sample_id))
}

if (EXAMPLE_SELECTION %in% c("median_profile", "median_profile_size_matched")) {
  prof   <- nerve_profiles()
  ranked <- lapply(cats, profile_ranked, prof = prof)
  
  if (EXAMPLE_SELECTION == "median_profile") {
    best <- lapply(ranked, function(df) df[1, ])
    cat("\n=== Examples: nerve whose profile (observed - own control) is closest to the category median profile ===\n")
  } else {
    # Pool: most typical TYPICAL_FRACTION of each category (at least 3 nerves)
    pools <- lapply(ranked, function(df) df[seq_len(min(nrow(df), max(3, ceiling(TYPICAL_FRACTION * nrow(df))))), ])
    cat("\n=== Size-matching pools (most typical nerves per category) ===\n")
    for (k in names(pools)) cat(cat_short[[k]], ": ", nrow(pools[[k]]), " of ", nrow(ranked[[k]]),
                                " nerves; spot counts ", paste(sort(pools[[k]]$n_spots), collapse = ", "), "\n", sep = "")
    # Combination with the most similar sizes (smallest log-size range);
    # ties broken by the sum of profile ranks (more typical first)
    combos <- expand.grid(a = seq_len(nrow(pools$adj)), n = seq_len(nrow(pools$ni)), v = seq_len(nrow(pools$invn)))
    sz <- cbind(pools$adj$n_spots[combos$a], pools$ni$n_spots[combos$n], pools$invn$n_spots[combos$v])
    combos$size_spread <- apply(log(sz), 1, function(x) diff(range(x)))
    combos$rank_sum    <- pools$adj$profile_rank[combos$a] + pools$ni$profile_rank[combos$n] + pools$invn$profile_rank[combos$v]
    top <- combos[order(round(combos$size_spread, 6), combos$rank_sum), ][1, ]
    best <- list(adj = pools$adj[top$a, ], ni = pools$ni[top$n, ], invn = pools$invn[top$v, ])
    cat("\n=== Examples: typical profiles, matched for size (spot counts ",
        paste(sz[as.integer(rownames(top)), ], collapse = " / "), ") ===\n", sep = "")
  }
  
} else if (EXAMPLE_SELECTION == "median_decay") {
  # Per category independently (rows may come from different patients).
  best <- lapply(cats, closest_to_median)
  cat("\n=== Examples: nerve closest to each category's median decay (near - far) ===\n")
  
} else if (EXAMPLE_SELECTION == "median_decay_same_patient") {
  # One patient for all three rows: the patient whose best-matching nerve in
  # each category is closest to that category's median, distances scaled by
  # the category IQR so no category dominates. Then that patient's nerve
  # closest to each category median.
  pats <- Reduce(intersect, lapply(cats, function(df) unique(df$patient_id)))
  if (length(pats) == 0) stop("No patient has paired nerves in all three categories.")
  dev_tab <- sapply(names(cats), function(k) {
    df <- cats[[k]]; med <- median(df$decay)
    sc <- IQR(df$decay); if (!is.finite(sc) || sc == 0) sc <- 1
    sapply(pats, function(pp) min(abs(df$decay[df$patient_id == pp] - med)) / sc)
  })
  dev_tab <- matrix(dev_tab, nrow = length(pats), dimnames = list(pats, names(cats)))
  patient_distance <- rowSums(dev_tab)
  chosen_patient <- names(which.min(patient_distance))
  cat("\n=== Patient choice: scaled distance of each patient's closest nerve from the category medians ===\n")
  print(round(cbind(dev_tab, total = patient_distance), 3))
  cat("Chosen patient:", chosen_patient, "\n")
  best <- lapply(cats, function(df) {
    sub <- df[df$patient_id == chosen_patient, ]
    sub[which.min(abs(sub$decay - median(df$decay))), ]
  })
  cat("\n=== Examples from patient", chosen_patient, ": nerve closest to each category's median decay ===\n")
  
} else {
  # "clearest_ordered": the ORIGINAL selection, kept for comparison only.
  # Picks visually clear examples that preserve NI-TB > InvN > NI-Adj; the
  # legend must then say examples were selected for clarity.
  ENRICHMENT_GROUPS <- c("ni", "invn")
  scored <- lapply(names(cats), function(k) {
    df <- cats[[k]]
    df$score <- if (k %in% ENRICHMENT_GROUPS) df$near_sep - abs(df$far_div) else -(abs(df$near_sep) + abs(df$far_div))
    df
  })
  names(scored) <- names(cats)
  clean <- lapply(scored, function(df) df %>% filter(score >= median(score, na.rm = TRUE)))
  combos <- expand.grid(ai = seq_len(nrow(clean$adj)), ni = seq_len(nrow(clean$ni)), vi = seq_len(nrow(clean$invn)))
  combos$min_score <- pmin(clean$adj$score[combos$ai], clean$ni$score[combos$ni], clean$invn$score[combos$vi])
  valid <- combos[clean$ni$near_sep[combos$ni] > clean$invn$near_sep[combos$vi] &
                    clean$invn$near_sep[combos$vi] > clean$adj$near_sep[combos$ai], ]
  valid <- valid[order(-valid$min_score), ]
  if (nrow(valid) == 0) {
    warning("No combination among 'reasonably clean' candidates preserves the expected ordering.")
    best <- lapply(scored, function(df) df[which.max(df$score), ])
  } else {
    best <- list(adj = clean$adj[valid$ai[1], ], ni = clean$ni[valid$ni[1], ], invn = clean$invn[valid$vi[1], ])
  }
  best <- lapply(best, function(d) d[, setdiff(names(d), "score")])
  cat("\n=== Clearest examples preserving the expected ordering (NOT representative) ===\n")
}
print(bind_rows(best, .id = "category"))

# ---- Representativeness check (all modes) ---------------------------------
# For each example: its decay (near - far) and percentile within its
# category, the category median and range, and -- for reference -- its
# sham-corrected enrichment percentile. Report the percentiles in the legend.
# Distance (um) from a nerve to the nearest spot of any other nerve in its section
nearest_other_nerve_um <- function(nid, smp) {
  this   <- coords[which(spe$sample_id == smp & !is.na(spe$NerveID) & spe$NerveID == nid), , drop = FALSE]
  others <- which(spe$sample_id == smp & !is.na(spe$NerveID) & spe$NerveID != nid)
  if (nrow(this) == 0 || length(others) == 0) return(NA_real_)
  min(RANN::nn2(this, coords[others, , drop = FALSE], k = 1)$nn.dists[, 1]) * um_per_pixel
}
representativeness <- bind_rows(lapply(names(cats), function(k) {
  df <- cats[[k]]; ex <- best[[k]]
  data.frame(category = cat_short[[k]], NerveID = ex$NerveID, patient_id = ex$patient_id,
             example_decay = ex$decay,
             decay_percentile = round(100 * mean(df$decay <= ex$decay)),
             category_median_decay = median(df$decay),
             category_min_decay = min(df$decay), category_max_decay = max(df$decay),
             example_sham_corrected = ex$delta_decay,
             sham_corrected_percentile = round(100 * mean(df$delta_decay <= ex$delta_decay)),
             nearest_other_nerve_um = nearest_other_nerve_um(ex$NerveID, ex$sample_id),
             n_spots = nerve_size(ex$NerveID, ex$sample_id),
             category_median_n_spots = median(mapply(nerve_size, df$NerveID, df$sample_id)),
             profile_rank = if (!is.null(ex$profile_rank)) ex$profile_rank else NA_integer_,
             n_nerves_in_category = nrow(df))
}))
cat("\n=== Representativeness check (selection mode: ", EXAMPLE_SELECTION, ") ===\n", sep = "")
print(representativeness, digits = 3, row.names = FALSE)
SUPP_TABLES_DIR <- file.path(PROJECT_DIR, "output", ANALYSIS_NAME, "tables")
dir.create(SUPP_TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(representativeness,
          file.path(SUPP_TABLES_DIR, paste0("zoom_examples_representativeness_", EXAMPLE_SELECTION, ".csv")),
          row.names = FALSE)

EXAMPLE_ADJ_NERVE     <- best$adj$NerveID;   EXAMPLE_ADJ_SAMPLE   <- best$adj$sample_id;   EXAMPLE_ADJ_PATIENT   <- best$adj$patient_id
EXAMPLE_NI_TB_NERVE   <- best$ni$NerveID;    EXAMPLE_NI_TB_SAMPLE <- best$ni$sample_id;    EXAMPLE_NI_TB_PATIENT <- best$ni$patient_id
EXAMPLE_INVN_NERVE    <- best$invn$NerveID;  EXAMPLE_INVN_SAMPLE  <- best$invn$sample_id;  EXAMPLE_INVN_PATIENT  <- best$invn$patient_id

# ---- 3. Distance-curve data ------------------------------------------------
# Real: per-bin mean around the observed nerve. Relocated: per-bin mean for
# each replicate, then averaged over replicates (= the control used in the
# statistics). Bins beyond CURVE_MAX_UM are NA in the bundle and dropped.
compute_nerve_curve_df <- function(nerve_id, sample_id, per_nerve_dist_all, sham_shapes_df,
                                   per_shape_dist_all, plasma_score, spe) {
  real_candidates <- per_nerve_dist_all %>% filter(NerveID == nerve_id, !is.na(bin))
  barcode_samples <- spe$sample_id[match(real_candidates$barcode, spe$barcode)]
  real_d <- real_candidates[barcode_samples == sample_id, ] %>%
    mutate(score = plasma_score[barcode]) %>%
    group_by(bin) %>% summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    mutate(Source = "Real")
  
  shape_ids <- unique(sham_shapes_df$ShamNerveID[sham_shapes_df$source_NerveID == nerve_id &
                                                   sham_shapes_df$sample_id == sample_id])
  reloc_d <- per_shape_dist_all %>% filter(ShamNerveID %in% shape_ids, !is.na(bin)) %>%
    group_by(ShamNerveID, bin) %>% summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    group_by(bin) %>% summarise(mean_score = mean(mean_score, na.rm = TRUE), .groups = "drop") %>%
    mutate(Source = "Relocated")
  
  curve_df <- bind_rows(real_d, reloc_d)
  curve_df$Source <- factor(curve_df$Source, levels = c("Real", "Relocated"))
  curve_df
}

# ---- 4. Zoomed observed-vs-relocated builder (per category) --------------
build_zoom_row <- function(nerve_id, sample_id, spe, coords, plasma_score, sham_shapes_df,
                           um_per_pixel, crop_margin_um, n_overlay, global_limits, row_label,
                           per_nerve_dist_all, per_shape_dist_all, curve_y_limits) {
  nerve_idx <- which(spe$sample_id == sample_id & !is.na(spe$NerveID) & spe$NerveID == nerve_id)
  if (length(nerve_idx) == 0) stop("No spots for NerveID=", nerve_id, ", sample_id=", sample_id)
  nerve_coords <- coords[nerve_idx, , drop = FALSE]
  
  this_nerve_shams <- sham_shapes_df %>% filter(source_NerveID == nerve_id, sample_id == !!sample_id)
  if (nrow(this_nerve_shams) == 0) stop("No sham shapes for NerveID=", nerve_id, ", sample_id=", sample_id)
  
  selection <- select_farthest_replicate(nerve_id, sample_id, nerve_coords, sham_shapes_df, coords, spe, um_per_pixel)
  featured_shape_id <- selection$shape_id
  cat("  [", row_label, "] observed-relocated centroid separation:", round(selection$separation_um), "um",
      "(crop diameter ~", round(2 * crop_margin_um), "um)\n", sep = "")
  
  featured_coords <- coords[match(this_nerve_shams$barcode[this_nerve_shams$ShamNerveID == featured_shape_id], spe$barcode), , drop = FALSE]
  featured_sample <- this_nerve_shams$sample_id[this_nerve_shams$ShamNerveID == featured_shape_id][1]
  
  other_shape_ids <- setdiff(unique(this_nerve_shams$ShamNerveID), featured_shape_id)
  overlay_outlines <- lapply(head(other_shape_ids, n_overlay), function(sid) {
    sc <- coords[match(this_nerve_shams$barcode[this_nerve_shams$ShamNerveID == sid], spe$barcode), , drop = FALSE]
    cbind(convex_hull_df(sc), shape_id = sid)
  })
  overlay_df <- if (length(overlay_outlines) > 0) bind_rows(overlay_outlines) else NULL
  
  p_real <- build_crop_panel_with_scalebar(nerve_coords, spe, coords, plasma_score, sample_id,
                                           um_per_pixel, crop_margin_um, global_limits, "Observed nerve")
  p_relocated <- build_crop_panel_with_scalebar(featured_coords, spe, coords, plasma_score, featured_sample,
                                                um_per_pixel, crop_margin_um, global_limits,
                                                "Relocated nerve geometry", extra_outlines = overlay_df)
  # Wrap each composite (spatial + scale bar) so patchwork treats it as ONE
  # panel in the widths = c(1, 1, 1.1) layout below.
  p_real <- wrap_elements(full = p_real)
  p_relocated <- wrap_elements(full = p_relocated)
  
  curve_df <- compute_nerve_curve_df(nerve_id, sample_id, per_nerve_dist_all, sham_shapes_df,
                                     per_shape_dist_all, plasma_score, spe)
  n_reps <- n_distinct(this_nerve_shams$ShamNerveID)
  
  p_curve <- ggplot(curve_df, aes(x = bin, y = mean_score, linetype = Source, group = Source)) +
    geom_line(color = "black", linewidth = 0.8) +
    geom_point(color = "black", size = 1.3) +
    scale_linetype_manual(values = c(Real = "solid", Relocated = "dashed"),
                          labels = c(Real = "Observed", Relocated = paste0("Relocated (mean of ", n_reps, ")"))) +
    scale_x_continuous(limits = c(0, CURVE_MAX_UM), breaks = bin_breaks) +
    coord_cartesian(ylim = curve_y_limits) +
    labs(title = "Distance profile", x = "Distance (µm, bin midpoints)", y = "Plasma-cell RCTD proportion",
         linetype = NULL) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 11))
  
  row_raw <- p_real + p_relocated + p_curve +
    plot_layout(widths = c(1, 1, 1.1)) +
    plot_annotation(title = row_label, theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.02)))
  wrap_elements(full = row_raw)
}

# ---- 5. Shared scales across all three zoom rows ---------------------------
nerve_coords_adj  <- coords[which(spe$sample_id == EXAMPLE_ADJ_SAMPLE   & spe$NerveID == EXAMPLE_ADJ_NERVE), , drop = FALSE]
nerve_coords_ni   <- coords[which(spe$sample_id == EXAMPLE_NI_TB_SAMPLE & spe$NerveID == EXAMPLE_NI_TB_NERVE), , drop = FALSE]
nerve_coords_invn <- coords[which(spe$sample_id == EXAMPLE_INVN_SAMPLE  & spe$NerveID == EXAMPLE_INVN_NERVE), , drop = FALSE]

sel_adj  <- select_farthest_replicate(EXAMPLE_ADJ_NERVE,   EXAMPLE_ADJ_SAMPLE,   nerve_coords_adj,  sham_shapes_df, coords, spe, um_per_pixel)
sel_ni   <- select_farthest_replicate(EXAMPLE_NI_TB_NERVE, EXAMPLE_NI_TB_SAMPLE, nerve_coords_ni,   sham_shapes_df, coords, spe, um_per_pixel)
sel_invn <- select_farthest_replicate(EXAMPLE_INVN_NERVE,  EXAMPLE_INVN_SAMPLE,  nerve_coords_invn, sham_shapes_df, coords, spe, um_per_pixel)

get_featured_coords_sample <- function(sel) {
  shams <- sel$this_nerve_shams
  fc <- coords[match(shams$barcode[shams$ShamNerveID == sel$shape_id], spe$barcode), , drop = FALSE]
  fs <- shams$sample_id[shams$ShamNerveID == sel$shape_id][1]
  list(coords = fc, sample = fs)
}
feat_adj  <- get_featured_coords_sample(sel_adj)
feat_ni   <- get_featured_coords_sample(sel_ni)
feat_invn <- get_featured_coords_sample(sel_invn)

ZOOM_COLOR_LIMITS <- range(c(
  get_crop_score_range(nerve_coords_adj,  spe, coords, plasma_score, EXAMPLE_ADJ_SAMPLE,   um_per_pixel, CROP_MARGIN_UM),
  get_crop_score_range(feat_adj$coords,   spe, coords, plasma_score, feat_adj$sample,      um_per_pixel, CROP_MARGIN_UM),
  get_crop_score_range(nerve_coords_ni,   spe, coords, plasma_score, EXAMPLE_NI_TB_SAMPLE, um_per_pixel, CROP_MARGIN_UM),
  get_crop_score_range(feat_ni$coords,    spe, coords, plasma_score, feat_ni$sample,       um_per_pixel, CROP_MARGIN_UM),
  get_crop_score_range(nerve_coords_invn, spe, coords, plasma_score, EXAMPLE_INVN_SAMPLE,  um_per_pixel, CROP_MARGIN_UM),
  get_crop_score_range(feat_invn$coords,  spe, coords, plasma_score, feat_invn$sample,     um_per_pixel, CROP_MARGIN_UM)))

curve_df_adj  <- compute_nerve_curve_df(EXAMPLE_ADJ_NERVE,   EXAMPLE_ADJ_SAMPLE,   per_nerve_dist_all, sham_shapes_df, per_shape_dist_all, plasma_score, spe)
curve_df_ni   <- compute_nerve_curve_df(EXAMPLE_NI_TB_NERVE, EXAMPLE_NI_TB_SAMPLE, per_nerve_dist_all, sham_shapes_df, per_shape_dist_all, plasma_score, spe)
curve_df_invn <- compute_nerve_curve_df(EXAMPLE_INVN_NERVE,  EXAMPLE_INVN_SAMPLE,  per_nerve_dist_all, sham_shapes_df, per_shape_dist_all, plasma_score, spe)
CURVE_Y_LIMITS <- range(c(curve_df_adj$mean_score, curve_df_ni$mean_score, curve_df_invn$mean_score), na.rm = TRUE)

# ---- 6. Build the three zoom rows ------------------------------------------
row_adj  <- build_zoom_row(EXAMPLE_ADJ_NERVE, EXAMPLE_ADJ_SAMPLE, spe, coords, plasma_score, sham_shapes_df,
                           um_per_pixel, CROP_MARGIN_UM, N_RELOCATED_OVERLAY, ZOOM_COLOR_LIMITS,
                           paste0("NI -- adjacent tissue (Patient ", EXAMPLE_ADJ_PATIENT, ")"),
                           per_nerve_dist_all, per_shape_dist_all, CURVE_Y_LIMITS)
row_ni   <- build_zoom_row(EXAMPLE_NI_TB_NERVE, EXAMPLE_NI_TB_SAMPLE, spe, coords, plasma_score, sham_shapes_df,
                           um_per_pixel, CROP_MARGIN_UM, N_RELOCATED_OVERLAY, ZOOM_COLOR_LIMITS,
                           paste0("NI -- tumour-bearing tissue (Patient ", EXAMPLE_NI_TB_PATIENT, ")"),
                           per_nerve_dist_all, per_shape_dist_all, CURVE_Y_LIMITS)
row_invn <- build_zoom_row(EXAMPLE_INVN_NERVE, EXAMPLE_INVN_SAMPLE, spe, coords, plasma_score, sham_shapes_df,
                           um_per_pixel, CROP_MARGIN_UM, N_RELOCATED_OVERLAY, ZOOM_COLOR_LIMITS,
                           paste0("InvN (Patient ", EXAMPLE_INVN_PATIENT, ")"),
                           per_nerve_dist_all, per_shape_dist_all, CURVE_Y_LIMITS)

zoom_legend_plot <- ggplot(data.frame(x = 0, y = 0, score = ZOOM_COLOR_LIMITS), aes(x = x, y = y, color = score)) +
  geom_point() + scale_color_viridis_c(option = "inferno", limits = ZOOM_COLOR_LIMITS) +
  labs(color = "Plasma cell\nweight (RCTD)") + theme(legend.position = "right")
zoom_legend <- cowplot::get_legend(zoom_legend_plot)

# ---- 7. Assemble figure ------------------------------------------------------
zoom_block <- wrap_elements(full = (row_adj / row_ni / row_invn)) | wrap_elements(full = zoom_legend)
full_fig <- zoom_block + plot_layout(widths = c(20, 2.2))

print(full_fig)
out_name <- paste0("real_vs_relocated_zoom_only_CALIBRATED_", EXAMPLE_SELECTION, ".pdf")
ggsave(file.path(FIGURES_DIR, out_name), full_fig, width = 13, height = 11, dpi = 600, device = pdf_device)

cat("\nSaved ", out_name, " (example selection: ", EXAMPLE_SELECTION, ").\n", sep = "")
