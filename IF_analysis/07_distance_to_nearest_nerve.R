###############################################################################
# 07_distance_to_nearest_nerve.R
#
# Distance from each annotated tumour ROI to the nearest annotated nerve, using
# the two cleaned exports together:
#   CANCER_FILE  cells of the perineural / distal tumour ROIs (01_clean_cancer.R)
#   NERVE_FILE   cells of the invaded / uninvaded nerve ROIs (02_clean_nerve.R)
#
# WHY
#   Distal tumour ROIs were selected visually as lying away from any viable
#   nerve. This measures how far they actually are, so the criterion can be
#   stated as a number rather than an impression, and so the imaging definition
#   of "distal" can be compared with the >= 1.3 mm threshold used for Tum-DistN
#   in the Visium analysis. Perineural ROIs are measured the same way and serve
#   as the positive control: their distances should be small.
#
# HOW
#   Both exports carry whole-image cell coordinates in um in the same frame, so
#   for every tumour cell the nearest nerve cell within the same image is found
#   (RANN::nn2, Euclidean). Distances are then summarised PER ROI, so each
#   patch - C_1, C_FAR_1 and so on - contributes one value rather than one per
#   cell. Distance to the nearest invaded and nearest uninvaded nerve are
#   reported separately alongside the distance to the nearest nerve of any kind.
#
# WHAT IT IS NOT
#   * Cell-centroid to cell-centroid, so edge-to-edge separation is shorter by
#     roughly one cell radius (~5 um) - negligible at these distances.
#   * Only nerves annotated in NERVE_FILE are counted. A nerve present in the
#     tissue but never annotated cannot contribute, so these are distances to
#     the nearest ANNOTATED nerve.
#   * An image present in only one export gets no distances and is reported.
#
# OUTPUTS (OUT_DIR)
#   Distance_per_ROI.csv        one row per tumour ROI
#   Distance_summary_by_type.csv   Distal vs Perineural, overall and per image
#   Distance_per_cell.csv.gz    optional, one row per tumour cell
#   Fig_Distance_per_ROI.pdf/.png  per-ROI distances by ROI type
#   Distance_report.txt         the console report
###############################################################################

library(RANN)

# ---- Settings ---------------------------------------------------------------
source(file.path(if (basename(getwd()) == "scripts") "." else "scripts", "00_settings.R"))

# Which cleaned exports to use is set once in 00_settings.R (CANCER_INPUT /
# NERVE_INPUT). The two exports must come from the same pair of folders, because
# Image_Name is the key that matches tumour ROIs to nerve ROIs.
CANCER_FILE <- cleaned_file(CANCER_INPUT)
NERVE_FILE  <- cleaned_file(NERVE_INPUT)
OUT_DIR     <- file.path("output", "Distance_to_nerve")

XCOL <- "Cells/X Position in um"
YCOL <- "Cells/Y Position in um"

CANCER_TYPES <- c("Distal", "Perineural")     # ROI types in the cancer export
NERVE_TYPES  <- c("Invaded", "Uninvaded")     # ROI types in the nerve export

# The distal-selection criterion to check against, in um. Reported, not applied.
CRITERION_UM   <- 500
# Visium Tum-DistN threshold, for comparison in the report.
VISIUM_UM      <- 1300
WRITE_PER_CELL <- FALSE   # TRUE also writes the per-cell distances (large)
MAKE_FIGURE    <- TRUE

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
log_lines <- character(0)
say <- function(...) {
  s <- paste0(...)
  log_lines <<- c(log_lines, s)
  cat(s, "\n", sep = "")
}

# ---- 1. Read both exports ---------------------------------------------------
need <- function(d, cols, what) {
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(what, " is missing column(s): ", paste(miss, collapse = "; "))
}

read_cells <- function(path, types, what) {
  if (!file.exists(path)) stop(what, " not found: ", path)
  d <- read.csv(path, check.names = FALSE)
  need(d, c("Image_Name", "ROI", "Type", XCOL, YCOL), what)
  keep <- c("Image_Name", "ROI", "Type", XCOL, YCOL,
            intersect(c("Patient_ID", "Block", "ROI_number"), names(d)))
  d <- d[, keep, drop = FALSE]
  d$x <- suppressWarnings(as.numeric(d[[XCOL]]))
  d$y <- suppressWarnings(as.numeric(d[[YCOL]]))
  bad <- is.na(d$x) | is.na(d$y)
  if (any(bad)) say("  ", what, ": ", sum(bad), " of ", nrow(d),
                    " cells dropped for missing coordinates")
  d <- d[!bad, ]
  other <- setdiff(unique(d$Type), types)
  if (length(other)) say("  ", what, ": ROI types not used: ", paste(other, collapse = ", "))
  d <- d[d$Type %in% types, ]
  say("  ", what, ": ", nrow(d), " cells, ", length(unique(d$Image_Name)), " images, types ",
      paste(sort(unique(d$Type)), collapse = " / "))
  d
}

say("Reading exports")
can <- read_cells(CANCER_FILE, CANCER_TYPES, "cancer export")
ner <- read_cells(NERVE_FILE,  NERVE_TYPES,  "nerve export")

# ---- 2. Which images can be matched ----------------------------------------
ci <- sort(unique(can$Image_Name)); ni <- sort(unique(ner$Image_Name))
both <- intersect(ci, ni)
say("\nImages in both exports (", length(both), "): ", paste(both, collapse = ", "))
if (length(setdiff(ci, ni)))
  say("Cancer export only, no distances computed: ", paste(setdiff(ci, ni), collapse = ", "))
if (length(setdiff(ni, ci)))
  say("Nerve export only, not used: ", paste(setdiff(ni, ci), collapse = ", "))
if (!length(both)) stop("No image name appears in both exports - check Image_Name spelling.")

# Coordinate frames must agree. Overlapping ranges are necessary, not sufficient,
# but a disjoint range means the two exports are not registered to each other.
say("\nCoordinate ranges per image (um), cancer vs nerve:")
for (im in both) {
  a <- can[can$Image_Name == im, ]; b <- ner[ner$Image_Name == im, ]
  ov <- function(p, q) max(0, min(max(p), max(q)) - max(min(p), min(q)))
  frac <- ov(a$x, b$x) / max(diff(range(a$x)), 1) * ov(a$y, b$y) / max(diff(range(a$y)), 1)
  say(sprintf("  %-14s cancer x %5.0f-%5.0f y %5.0f-%5.0f | nerve x %5.0f-%5.0f y %5.0f-%5.0f%s",
              im, min(a$x), max(a$x), min(a$y), max(a$y),
              min(b$x), max(b$x), min(b$y), max(b$y),
              if (frac < 0.25) "   << LITTLE OVERLAP, CHECK REGISTRATION" else ""))
}

# ---- 3. Nearest nerve cell, per tumour cell --------------------------------
nearest <- function(target, query) {            # returns distance and row index
  if (!nrow(target)) return(list(d = rep(NA_real_, nrow(query)), i = rep(NA_integer_, nrow(query))))
  nn <- nn2(as.matrix(target[, c("x", "y")]), as.matrix(query[, c("x", "y")]), k = 1)
  list(d = nn$nn.dists[, 1], i = nn$nn.idx[, 1])
}

per_cell <- do.call(rbind, lapply(both, function(im) {
  q  <- can[can$Image_Name == im, ]
  nv <- ner[ner$Image_Name == im, ]
  any_ <- nearest(nv, q)
  inv  <- nearest(nv[nv$Type == "Invaded", ],   q)
  uni  <- nearest(nv[nv$Type == "Uninvaded", ], q)
  nearest_roi  <- if (all(is.na(any_$i))) NA_character_ else nv$ROI[any_$i]
  nearest_type <- if (all(is.na(any_$i))) NA_character_ else nv$Type[any_$i]
  data.frame(
    Image_Name = im,
    Patient_ID = if ("Patient_ID" %in% names(q)) q$Patient_ID else NA,
    ROI = q$ROI, ROI_type = q$Type,
    d_any = any_$d, d_invaded = inv$d, d_uninvaded = uni$d,
    Nearest_nerve_ROI = nearest_roi, Nearest_nerve_type = nearest_type,
    stringsAsFactors = FALSE)
}))

if (WRITE_PER_CELL) {
  gz <- gzfile(file.path(OUT_DIR, "Distance_per_cell.csv.gz"), "w")
  write.csv(per_cell, gz, row.names = FALSE); close(gz)
}

# ---- 4. One value per ROI ---------------------------------------------------
# Each patch (C_1, C_FAR_1, ...) contributes one row. The mean and median
# describe the patch as a whole; the minimum is its nearest edge, which is the
# value the visual selection criterion is really about.
mode_chr <- function(v) { v <- v[!is.na(v)]; if (!length(v)) NA_character_ else names(sort(table(v), decreasing = TRUE))[1] }

per_roi <- do.call(rbind, lapply(
  split(per_cell, list(per_cell$Image_Name, per_cell$ROI), drop = TRUE),
  function(s) data.frame(
    Image_Name = s$Image_Name[1],
    Patient_ID = s$Patient_ID[1],
    ROI = s$ROI[1], ROI_type = s$ROI_type[1],
    N_cells = nrow(s),
    Min_distance_um    = round(min(s$d_any)),
    Median_distance_um = round(median(s$d_any)),
    Mean_distance_um   = round(mean(s$d_any)),
    Max_distance_um    = round(max(s$d_any)),
    Min_to_invaded_um   = if (all(is.na(s$d_invaded)))   NA else round(min(s$d_invaded,   na.rm = TRUE)),
    Min_to_uninvaded_um = if (all(is.na(s$d_uninvaded))) NA else round(min(s$d_uninvaded, na.rm = TRUE)),
    Nearest_nerve_ROI  = mode_chr(s$Nearest_nerve_ROI),
    Nearest_nerve_type = mode_chr(s$Nearest_nerve_type),
    stringsAsFactors = FALSE)))
per_roi <- per_roi[order(per_roi$ROI_type, per_roi$Image_Name, per_roi$Min_distance_um), ]
write.csv(per_roi, file.path(OUT_DIR, "Distance_per_ROI.csv"), row.names = FALSE)

# ---- 5. Report --------------------------------------------------------------
say("\n================ Distance to the nearest annotated nerve ================")
for (tp in intersect(CANCER_TYPES, unique(per_roi$ROI_type))) {
  s <- per_roi[per_roi$ROI_type == tp, ]
  say("\n", tp, " ROIs: ", nrow(s), " in ", length(unique(s$Image_Name)), " images")
  q <- function(v) sprintf("median %.0f, IQR %.0f-%.0f, range %.0f-%.0f um",
                           median(v), quantile(v, .25), quantile(v, .75), min(v), max(v))
  say("  nearest cell of the ROI : ", q(s$Min_distance_um))
  say("  ROI as a whole (median) : ", q(s$Median_distance_um))
  if (tp == "Distal") {
    say("  ROIs whose nearest cell is < ", CRITERION_UM, " um: ",
        sum(s$Min_distance_um < CRITERION_UM), " of ", nrow(s))
    say("  ROIs whose median is    < ", CRITERION_UM, " um: ",
        sum(s$Median_distance_um < CRITERION_UM), " of ", nrow(s))
    say("  ROIs whose median is    < ", VISIUM_UM, " um (Visium Tum-DistN threshold): ",
        sum(s$Median_distance_um < VISIUM_UM), " of ", nrow(s))
    close_rois <- s[s$Min_distance_um < CRITERION_UM, ]
    if (nrow(close_rois)) {
      say("  ROIs below the criterion at their nearest edge:")
      for (k in seq_len(nrow(close_rois)))
        say(sprintf("    %-14s %-10s min %4.0f, median %4.0f um, nearest nerve %s (%s)",
                    close_rois$Image_Name[k], close_rois$ROI[k],
                    close_rois$Min_distance_um[k], close_rois$Median_distance_um[k],
                    close_rois$Nearest_nerve_ROI[k], close_rois$Nearest_nerve_type[k]))
    }
  }
}

say("\nPer image, by ROI type (median over ROIs of the per-ROI median distance, um):")
tab <- tapply(per_roi$Median_distance_um,
              list(per_roi$Image_Name, per_roi$ROI_type), median)
print(round(tab)); log_lines <- c(log_lines, capture.output(print(round(tab))))

summary_by_type <- do.call(rbind, lapply(
  split(per_roi, list(per_roi$Image_Name, per_roi$ROI_type), drop = TRUE),
  function(s) data.frame(
    Image_Name = s$Image_Name[1], ROI_type = s$ROI_type[1], N_ROIs = nrow(s),
    Min_distance_median_um    = round(median(s$Min_distance_um)),
    Median_distance_median_um = round(median(s$Median_distance_um)),
    Median_distance_min_um    = round(min(s$Median_distance_um)),
    Median_distance_max_um    = round(max(s$Median_distance_um)),
    stringsAsFactors = FALSE)))
write.csv(summary_by_type, file.path(OUT_DIR, "Distance_summary_by_type.csv"), row.names = FALSE)

# Sanity check: perineural ROIs must be closer to nerves than distal ones.
if (all(c("Distal", "Perineural") %in% per_roi$ROI_type)) {
  pm <- median(per_roi$Median_distance_um[per_roi$ROI_type == "Perineural"])
  dm <- median(per_roi$Median_distance_um[per_roi$ROI_type == "Distal"])
  say(sprintf("\nControl check: perineural ROIs sit at a median of %.0f um from a nerve, distal ROIs at %.0f um (%.0f-fold).",
              pm, dm, dm / max(pm, 1)))
  if (pm >= dm) say("  ! Perineural ROIs are not closer than distal ones - check the registration and ROI types.")
}

# ---- 6. Figure --------------------------------------------------------------
if (MAKE_FIGURE && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  g <- ggplot(per_roi, aes(ROI_type, Median_distance_um)) +
    geom_hline(yintercept = c(CRITERION_UM, VISIUM_UM), linetype = c("dotted", "dashed"),
               colour = "grey50", linewidth = 0.3) +
    geom_point(aes(colour = Image_Name), position = position_jitter(width = 0.15, height = 0, seed = 1),
               size = 1.6, alpha = 0.9) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "black") +
    scale_y_log10() +
    labs(x = NULL, y = expression("Distance to nearest annotated nerve (" * mu * "m, log scale)"),
         colour = "Image",
         caption = paste0("One point per ROI (median over its cells). Dotted line ", CRITERION_UM,
                          " um, dashed line ", VISIUM_UM, " um.")) +
    theme_classic(base_size = 8) +
    theme(legend.position = "right", plot.caption = element_text(hjust = 0, size = 6))
  ggsave(file.path(OUT_DIR, "Fig_Distance_per_ROI.pdf"), g, width = 120, height = 70, units = "mm")
  ggsave(file.path(OUT_DIR, "Fig_Distance_per_ROI.png"), g, width = 120, height = 70, units = "mm",
         dpi = 600, bg = "white")
}

writeLines(log_lines, file.path(OUT_DIR, "Distance_report.txt"))
say("\nWritten to ", normalizePath(OUT_DIR, mustWork = FALSE))
