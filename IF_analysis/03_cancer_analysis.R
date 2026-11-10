###############################################################################
# Cancer analysis: PanCK+ tumour cells, PERINEURAL vs DISTAL ROIs
#
# Inputs (paths and input folder set in 00_settings.R)
#   output/<CANCER_INPUT>/cleaned/combined/All_images_combined.csv  (01_clean_cancer.R)
#   input/thresholds.csv - columns Image_Name and PanCK_threshold (other columns and
#                     extra blocks are ignored). Image_Name must be written EXACTLY
#                     as in the data, e.g. "23885 bl5" (the script prints the names)
# Older runs exported some measures under other names (col_aliases) and may lack
# markers (70865 bl15: only PanCK and GPX2). A block without data for a marker is
# left out of that marker only; the script prints which blocks each marker uses.
#
# Tumour cells : PanCK mean intensity > PanCK_threshold of that image
# Per ROI      : RNA     = mean dots per PanCK+ cell
#                protein = mean intensity of PanCK+ cells
#                background = mean CD3 (TRITC) / CD45 (Cy5) of the same cells
#                corrected protein = protein / same-channel background
#                (ROIs with < min_cells PanCK+ cells are excluded and listed)
# log2         : applied to the ROI values (after averaging cells, never to cells);
#                ROI values of 0 (GPX2) get half the smallest non-zero value first
# Statistics   : ROI = observation, patient = blocking factor
#                lmer(log2 ROI value ~ Region + (1 | Patient)), equal weight per ROI
#                -> difference perineural vs distal WITHIN the analysed tumours
#                fold change = 2^coefficient (perineural / distal), 95% CI
#                Holm correction within each figure; leave-one-patient-out
#                Supplementary columns: patient as replicate
#                lmer(log2 ROI value ~ Region + (1 | Patient) + (0 + Region | Patient))
#                (one block per patient - the script checks this)
# Cohort vs discovery patient: blocks in discovery_blocks (the Visium patient) are
#                NOT in the cohort statistics; they are shown in a separate figure
#                (each ROI one dot, perineural vs distal, Wilcoxon rank-sum test on
#                ROIs, which is descriptive for that single patient)
# Figure       : circles = per patient, median of log2 perineural ROI values minus
#                median of log2 distal ROI values (= log2 fold change), shown as ratio
# CD3/CD45 background: only with the "_with_background" input; without these columns
#   the background and corrected analyses are skipped (reported when the script runs)
# Patient colours: 00_settings.R (same colours as the QuPath/COMET figures)
# Outputs (folder output/Cancer_final)
#   Fig_Cohort_forest            cohort: RNA + protein (10 markers)
#   Fig_Cohort_<marker>_superplot  separate SuperPlot for the markers in separate_markers
#   Fig_Discovery_<marker>       discovery patient: per-ROI dot plot
#   Cohort_statistics.csv, Cohort_per_patient.csv, Discovery_<marker>_test.csv
#   Background_statistics.csv, Corrected_statistics.csv (CD3/CD45 controls and
#   background-corrected protein; figures only if plot_background_figures = TRUE)
#   Per_ROI_values.csv, PanCK_threshold_QC.png, Figure_legends.txt
# Blocks without both perineural and distal ROIs (e.g. 80830 bl3, distal only) are
# left out of the cohort analysis automatically (reported when the script runs).
###############################################################################
library(ggplot2)
library(lmerTest)   # lme4 + Satterthwaite degrees of freedom

# ---- Settings ---------------------------------------------------------------
source(file.path(if (basename(getwd()) == "scripts") "." else "scripts", "00_settings.R"))
data_file      <- cleaned_file(CANCER_INPUT)
threshold_file <- THRESHOLD_FILE
out_dir        <- file.path("output", "Cancer_final")
min_cells      <- 10
panck_col      <- "Cells/Mean Intensity (Pan CK - Cy5 (1))"

# Other names of the same measure in older runs (standard name = other name(s));
# values are combined into the standard column
col_aliases <- list(
  "Cells/Mean Intensity (Pan CK - Cy5 (1))" = "Cells/Mean Intensity (Pan CK - TRITC_AF_PROT)",
  "Cells/Dot count (Dots T2 GPX2 - TRITC)"  = "Cells/Dot count (Dots T10 GPX2 - TRITC_RNA_AF)")

# Discovery (Visium) patient: shown separately, not part of the cohort statistics
discovery_blocks <- c("70865 bl15")
discovery_marker <- "GPX2"               # marker plotted for the discovery patient

# Blocks to leave out completely (e.g. for a sensitivity analysis)
exclude_blocks <- c()

# Background-control figures (CD3/CD45, corrected protein); the tables are always written
plot_background_figures <- FALSE

# Markers also shown as a separate SuperPlot (ROIs, perineural vs distal), in the same
# format as the QuPath tenascin / p75NTR figures; they stay in the forest plot as well
separate_markers <- c("ITGB4")
separate_p       <- "Holm"   # "Holm" (adjusted over the markers of the forest plot) or "raw"

# Patient colours: 00_settings.R (block_colours / block_labels)

# ---- Figure settings (edit here) ---------------------------------------------
# Sizes in mm; 1 column ~ 85-89 mm, 1.5 columns ~ 120-140 mm, 2 columns ~ 180-183 mm
fig_font         <- "Arial"        # "Arial" or "Helvetica" (most journals); any installed font
fig_text_pt      <- 7              # all text (axis labels, ticks, legend, p values), in pt
fig_width_mm     <- 180            # width of the supplementary forest plots
fig_height_mm    <- c(main = 120, background = 50, corrected = 75)
# Main forest plot (Fig_Cohort_forest): own size and text size
fig_forest_width_mm  <- 55
fig_forest_height_mm <- 72
fig_forest_text_pt   <- 8          # all text in the main forest plot, in pt
fig_forest_p         <- "holm"     # right-hand p values: "holm" (Holm only), "both" (raw + Holm), "none"
fig_forest_suffix    <- FALSE      # TRUE: "(RNA)" / "(protein)" after marker names
# FALSE: names only; RNA rows are shaded and RNA names in italics
fig_dpi          <- 600            # PNG / TIFF resolution
fig_formats      <- c("pdf", "png", "tiff")   # pdf = vector, fonts embedded
fig_show_title   <- FALSE          # titles/captions usually go in the figure legend text;
fig_show_caption <- FALSE          #   they are always written to Figure_legends.txt
fig_show_p       <- TRUE           # p values on the right-hand side
fig_point_mm     <- 1.4            # patient circles
fig_diamond_mm   <- 2.4            # overall fold change
fig_line_mm      <- 0.4            # CI bars and axis line
fig_x_limits     <- c(0.25, 4)     # fold-change axis; points outside are drawn at the edge
fig_x_breaks     <- c(0.25, 0.5, 1, 2, 4)
fig_disc_width_mm  <- 50           # discovery-patient figure
fig_disc_height_mm <- 60

markers <- read.table(header = TRUE, sep = "|", strip.white = TRUE, text = "
name   | type    | background | column
PFKFB3 | RNA     |            | Cells/Dot count (Dots T12 PFKFB3 - Cy7)
ENO1   | RNA     |            | Cells/Dot count (Dots T7 ENO1 - Cy5)
TFAM   | RNA     |            | Cells/Dot count (Dots T11 TFAM - Cy5)
MT-ND4 | RNA     |            | Cells/Dot count (Dots T8 MT_ND4 - Cy7)
GPX2   | RNA     |            | Cells/Dot count (Dots T2 GPX2 - TRITC)
TOMM20 | Protein | CD3        | Cells/Mean Intensity (TOMM20 - TRITC (1))
NDUFB8 | Protein | CD45       | Cells/Mean Intensity (NDUFB8 - Cy5 (1))
MTCO1  | Protein | CD45       | Cells/Mean Intensity (MTOC1 - Cy5 (1))
ITGB4  | Protein | CD3        | Cells/Mean Intensity (ITGB4 - TRITC (1))
GLUT1  | Protein | CD3        | Cells/Mean Intensity (GLUT1 - TRITC (1))")
background_cols <- c(CD3  = "Cells/Mean Intensity (CD3 - TRITC (1))",
                     CD45 = "Cells/Mean Intensity (CD45 - Cy5 (1))")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Read data and thresholds, check that everything matches --------------
d   <- read.csv(data_file, check.names = FALSE)
thr <- read.csv(threshold_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")

missing_cols <- setdiff(c("Image_Name", "Patient_ID", "ROI", "Type"), names(d))
if (length(missing_cols)) stop("Columns not in data: ", paste(missing_cols, collapse = "; "))
# one column per measure: standard name, filled from its older-run names where empty
has_background <- all(sapply(background_cols, function(b)
  length(intersect(c(b, col_aliases[[b]]), names(d))) > 0))
if (!has_background) {
  message("No CD3/CD45 background columns in ", data_file,
          " -> background and corrected analyses skipped")
  background_cols <- background_cols[0]
}
for (std in c(panck_col, markers$column, background_cols)) {
  present <- intersect(c(std, col_aliases[[std]]), names(d))
  if (!length(present)) stop("Column not in data (under any name): ", std)
  v <- rep(NA_real_, nrow(d))
  for (cn in present) v[is.na(v)] <- d[[cn]][is.na(v)]
  d[[std]] <- v
}
if (length(exclude_blocks)) {
  message("Blocks excluded by setting: ", paste(exclude_blocks, collapse = ", "))
  d <- d[!d$Image_Name %in% exclude_blocks, ]
}
if (!all(c("Image_Name", "PanCK_threshold") %in% names(thr)))
  stop(threshold_file, " needs the columns Image_Name and PanCK_threshold")

dup <- thr$Image_Name[duplicated(thr$Image_Name)]
if (length(dup)) stop("Image_Name listed more than once in ", threshold_file, ": ",
                      paste(dup, collapse = ", "))
data_images <- sort(unique(d$Image_Name))
thr_images  <- thr$Image_Name[!is.na(thr$PanCK_threshold)]
no_thr <- setdiff(data_images, thr_images)
if (length(no_thr)) stop(
  "\nNo PanCK_threshold for: ", paste0('"', no_thr, '"', collapse = ", "),
  "\nImage names in the data:      ", paste0('"', data_images, '"', collapse = ", "),
  "\nImage names in ", threshold_file, ": ", paste0('"', thr$Image_Name, '"', collapse = ", "),
  "\nWrite the names in ", threshold_file, " exactly as in the data.")
blk <- unique(d[, c("Image_Name", "Patient_ID")])
block_cols <- block_colours(blk$Image_Name, blk$Patient_ID)   # 00_settings.R
block_labs <- block_labels(blk$Image_Name, blk$Patient_ID)
unused <- setdiff(thr_images, data_images)
if (length(unused)) message("Thresholds not used (image not in data): ",
                            paste(unused, collapse = ", "))

d <- d[d$Type %in% c("Distal", "Perineural"), ]
d$PanCK_threshold <- thr$PanCK_threshold[match(d$Image_Name, thr$Image_Name)]
d$PanCK_pos <- d[[panck_col]] > d$PanCK_threshold

counts <- aggregate(cbind(Cells = 1, PanCK_pos = PanCK_pos) ~ Image_Name + PanCK_threshold,
                    data = d, FUN = sum)
counts$Pct_PanCK_pos <- round(100 * counts$PanCK_pos / counts$Cells, 1)
write.csv(counts, file.path(out_dir, "PanCK_cells_per_image.csv"), row.names = FALSE)
cat("\nPanCK+ cells per image:\n"); print(counts, row.names = FALSE)

qc <- ggplot(d, aes(.data[[panck_col]])) +
  geom_histogram(bins = 80, fill = "grey60") +
  geom_vline(aes(xintercept = PanCK_threshold), colour = "#eb6834", linewidth = 0.8) +
  facet_wrap(~ Image_Name, scales = "free") +
  labs(x = "PanCK mean intensity per cell", y = "Cells",
       title = "PanCK threshold per image (orange)") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "PanCK_threshold_QC.png"), qc, width = 11, height = 6,
       dpi = 200, bg = "white")

tumour <- d[d$PanCK_pos, ]

# ---- 2. Per-ROI values --------------------------------------------------------
summarise_min <- function(v, f) { v <- v[!is.na(v)]; if (length(v) >= min_cells) f(v) else NA }

roi <- do.call(rbind, lapply(split(tumour, list(tumour$Image_Name, tumour$ROI), drop = TRUE),
                             function(x) {
                               out <- data.frame(Patient_ID = x$Patient_ID[1], Image_Name = x$Image_Name[1],
                                                 ROI = x$ROI[1], Region = x$Type[1], N_PanCK_cells = nrow(x))
                               for (b in names(background_cols))
                                 out[[paste0(b, "_background")]] <- summarise_min(x[[background_cols[b]]], mean)
                               for (i in seq_len(nrow(markers))) {
                                 m <- markers[i, ]
                                 out[[m$name]] <- summarise_min(x[[m$column]], mean)
                                 if (m$type == "Protein" && has_background)
                                   out[[paste0(m$name, "_corrected")]] <- out[[m$name]] /
                                   out[[paste0(m$background, "_background")]]
                               }
                               out
                             }))
roi <- roi[order(roi$Image_Name, roi$Region, roi$ROI), ]
write.csv(roi, file.path(out_dir, "Per_ROI_values.csv"), row.names = FALSE)

small <- roi[roi$N_PanCK_cells < min_cells, c("Image_Name", "ROI", "Region", "N_PanCK_cells")]
write.csv(small, file.path(out_dir, "Excluded_ROIs_few_PanCK_cells.csv"), row.names = FALSE)
cat("\nROIs excluded (<", min_cells, "PanCK+ cells):", nrow(small), "\n")
if (nrow(small)) print(small, row.names = FALSE)
roi_n <- aggregate(ROI ~ Image_Name + Region, data = roi[roi$N_PanCK_cells >= min_cells, ],
                   FUN = length)
cat("\nROIs used per image:\n"); print(reshape(roi_n, idvar = "Image_Name",
                                               timevar = "Region", direction = "wide"),
                                       row.names = FALSE)

no_bg <- if (!has_background) roi[0, ] else roi[is.na(roi$CD3_background) & roi$N_PanCK_cells >= min_cells,
                                                c("Image_Name", "ROI", "Region", "N_PanCK_cells")]
if (nrow(no_bg)) {
  cat("\nROIs without CD3/CD45 values (not in the corrected analysis):\n")
  print(no_bg, row.names = FALSE)
}
roi$Region <- factor(roi$Region, levels = c("Distal", "Perineural"))
missing_disc <- setdiff(discovery_blocks, roi$Image_Name)
if (length(missing_disc)) message("Discovery block(s) not in data: ", paste(missing_disc, collapse = ", "))
roi_disc <- roi[roi$Image_Name %in% discovery_blocks & roi$N_PanCK_cells >= min_cells, ]
roi      <- roi[!roi$Image_Name %in% discovery_blocks, ]          # cohort from here on
# a block needs perineural AND distal ROIs to contribute to the comparison
regions_per_block <- tapply(as.character(roi$Region[roi$N_PanCK_cells >= min_cells]),
                            roi$Image_Name[roi$N_PanCK_cells >= min_cells],
                            function(r) all(c("Distal", "Perineural") %in% r))
one_region <- names(regions_per_block)[!regions_per_block]
one_region <- union(one_region, setdiff(unique(roi$Image_Name), names(regions_per_block)))
if (length(one_region)) {
  message("Not in the cohort analysis (no perineural or no distal ROIs): ",
          paste(one_region, collapse = ", "))
  roi <- roi[!roi$Image_Name %in% one_region, ]
}
cat("\nCohort blocks:", paste(sort(unique(roi$Image_Name)), collapse = ", "),
    "\nDiscovery block(s):", paste(sort(unique(roi_disc$Image_Name)), collapse = ", "), "\n")
blocks_per_patient <- tapply(roi$Image_Name, roi$Patient_ID, function(x) length(unique(x)))
if (any(blocks_per_patient > 1)) stop("More than one block for patient(s): ",
                                      paste(names(which(blocks_per_patient > 1)), collapse = ", "),
                                      " - the model assumes one block per patient")
images <- sort(unique(roi$Image_Name))
cols   <- block_cols[images]
avail <- sapply(c(markers$name, if (length(background_cols)) paste0(names(background_cols), "_background")),
                function(m) tapply(!is.na(roi[[m]]), roi$Image_Name, any))
cat("\nCohort - markers available per block (FALSE = block left out of that marker):\n")
print(t(avail))

# ---- 3. Statistics -------------------------------------------------------------
log_values <- function(df, v) {           # log2 ROI values (small constant if zeros, GPX2)
  df <- df[!is.na(df[[v]]), ]
  pc <- if (any(df[[v]] == 0)) min(df[[v]][df[[v]] > 0]) / 2 else 0
  df$y <- log2(df[[v]] + pc)
  df$peri <- as.numeric(df$Region == "Perineural")
  df
}
fit <- function(df, v, patient_slope = FALSE) {
  df <- log_values(df, v)
  form <- if (patient_slope) y ~ peri + (1 | Patient_ID) + (0 + peri | Patient_ID)
  else y ~ peri + (1 | Patient_ID)
  m  <- suppressMessages(lmer(form, data = df, REML = TRUE))
  cf <- summary(m)$coefficients["peri", ]
  ci <- cf[["Estimate"]] + c(-1, 1) * qt(0.975, cf[["df"]]) * cf[["Std. Error"]]
  c(fc = 2^cf[["Estimate"]], lo = 2^ci[1], hi = 2^ci[2], p = cf[["Pr(>|t|)"]],
    df = cf[["df"]], n_peri = sum(df$peri), n_dist = sum(1 - df$peri))
}
ratio_per_image <- function(v) {          # per patient: median log2 peri - median log2 distal
  df <- log_values(roi, v)
  sapply(images, function(im) {
    s <- df[df$Image_Name == im, ]
    2^(median(s$y[s$peri == 1]) - median(s$y[s$peri == 0]))
  })
}

analyse <- function(set, tag) {
  res <- lapply(seq_len(nrow(set)), function(i) {
    v <- set$var[i]; f <- fit(roi, v); r <- ratio_per_image(v)
    g <- fit(roi, v, patient_slope = TRUE)
    loo <- sapply(images, function(im) fit(roi[roi$Image_Name != im, ], v)[c("fc", "p")])
    list(stats = data.frame(
      Marker = set$label[i], Measure = set$measure[i],
      ROIs_perineural = f[["n_peri"]], ROIs_distal = f[["n_dist"]],
      Fold_change = round(f[["fc"]], 2),
      CI95 = paste(round(f[["lo"]], 2), "to", round(f[["hi"]], 2)),
      df = round(f[["df"]], 1), p = signif(f[["p"]], 2),
      Patients_same_direction = paste0(sum(sign(log(r)) == sign(log(f[["fc"]])), na.rm = TRUE),
                                       "/", sum(is.finite(r))),
      LOO_significant = paste0(sum(loo["p", ] < 0.05), "/", length(images)),
      LOO_fold_change_range = paste(round(min(loo["fc", ]), 2), "to",
                                    round(max(loo["fc", ]), 2)),
      Patient_as_replicate_fold_change = round(g[["fc"]], 2),
      Patient_as_replicate_CI95 = paste(round(g[["lo"]], 2), "to", round(g[["hi"]], 2)),
      Patient_as_replicate_p = signif(g[["p"]], 2),
      lo = f[["lo"]], hi = f[["hi"]]),
      per = data.frame(Marker = set$label[i], Image_Name = images,
                       Ratio_perineural_distal = round(r, 3)))
  })
  stats <- do.call(rbind, lapply(res, `[[`, "stats"))
  stats$p_Holm <- signif(p.adjust(stats$p, "holm"), 2)
  stats$Patient_as_replicate_p_Holm <- signif(p.adjust(stats$Patient_as_replicate_p, "holm"), 2)
  per <- do.call(rbind, lapply(res, `[[`, "per"))
  write.csv(stats[, setdiff(names(stats), c("lo", "hi"))],
            file.path(out_dir, paste0(tag, "_statistics.csv")), row.names = FALSE)
  write.csv(per, file.path(out_dir, paste0(tag, "_per_patient.csv")), row.names = FALSE)
  cat("\n==", tag, "==\n")
  print(stats[, c("Marker", "Fold_change", "CI95", "p", "p_Holm",
                  "Patients_same_direction", "LOO_significant")], row.names = FALSE)
  list(stats = stats, per = per)
}

# ---- 4. Forest plot ---------------------------------------------------------------
legend_common <- paste0(
  "Diamonds and bars: fold change (perineural / distal) across patients with 95% CI, from a ",
  "linear mixed model of log2 ROI values with patient as random intercept (filled: ",
  "Holm-adjusted p < 0.05). Circles: individual patients (median perineural ROI / median ",
  "distal ROI). p values: unadjusted (Holm-adjusted in brackets). Axis limited to ",
  fig_x_limits[1], "-", fig_x_limits[2], "; values outside are drawn at the axis edge.")
legends <- character(0)

save_fig <- function(p, file, height, width = fig_width_mm) {
  for (f in fig_formats) {
    out <- file.path(out_dir, paste0(file, ".", f))
    if (f == "pdf")  ggsave(out, p, device = cairo_pdf, width = width,
                            height = height, units = "mm")
    if (f == "png")  ggsave(out, p, width = width, height = height, units = "mm",
                            dpi = fig_dpi, bg = "white")
    if (f == "tiff") ggsave(out, p, device = "tiff", width = width, height = height,
                            units = "mm", dpi = fig_dpi, bg = "white", compression = "lzw")
  }
}

forest <- function(res, set, title, legend_text, file, height, width = fig_width_mm,
                   text_pt = fig_text_pt, p_mode = "both", compact = FALSE) {
  fig_text_pt <- text_pt                       # local text size for this figure
  legends[[file]] <<- paste0(title, ". ", legend_common, " ", legend_text)
  lab_order <- rev(set$label)
  fs <- res$stats; fs$yy <- match(fs$Marker, lab_order)
  fs$Sig <- factor(ifelse(fs$p_Holm < 0.05, "Holm p < 0.05", "n.s."),
                   levels = c("Holm p < 0.05", "n.s."))
  fmt <- function(x) ifelse(x < 0.001, "<0.001", trimws(formatC(x, format = "fg", digits = 2)))
  fs$p_txt <- if (p_mode == "holm") fmt(fs$p_Holm) else
    paste0(ifelse(fs$p < 0.001, "p < 0.001", paste0("p = ", fmt(fs$p))), " (", fmt(fs$p_Holm), ")")
  show_p <- fig_show_p && p_mode != "none"
  # patient legend: "Patient N" (block added only if a patient has more than one block here)
  leg_labs <- block_labs[names(cols)]
  if (compact) {
    short <- sub(" \\(.*\\)$", "", leg_labs)
    used  <- names(cols) %in% unique(res$per$Image_Name)
    dup   <- short %in% short[used][duplicated(short[used])]
    leg_labs <- ifelse(dup, leg_labs, short)
  }
  names(leg_labs) <- names(cols)
  leg_breaks <- names(cols)[order(block_labs[names(cols)])]   # ordered by patient number
  ital <- set$type[match(lab_order, set$label)] == "RNA"     # RNA names in italics (compact plot)
  ylab_txt <- if (compact) parse(text = ifelse(ital, paste0("italic('", lab_order, "')"),
                                               paste0("'", lab_order, "'"))) else lab_order
  pi <- res$per[is.finite(res$per$Ratio_perineural_distal) &
                  res$per$Ratio_perineural_distal > 0, ]
  pi$yy <- match(pi$Marker, lab_order)
  rna_rows <- fs[set$type == "RNA", ]          # grey shading behind RNA rows
  txt <- element_text(size = fig_text_pt, colour = "black", family = fig_font)
  ln  <- fig_line_mm / 0.75                    # ggplot linewidth unit ~0.75 mm
  
  p <- ggplot() +
    geom_rect(data = rna_rows, aes(ymin = yy - 0.5, ymax = yy + 0.5),
              xmin = -Inf, xmax = Inf, fill = "grey94") +
    geom_vline(xintercept = 1, colour = "grey50", linewidth = ln * 0.8) +
    geom_point(data = pi, aes(Ratio_perineural_distal, yy, colour = Image_Name),
               shape = 16, size = fig_point_mm, alpha = 0.9,
               position = position_jitter(height = 0.15, width = 0, seed = 3)) +
    geom_errorbarh(data = fs, aes(xmin = lo, xmax = hi, y = yy),
                   height = 0.3, linewidth = ln, colour = "black") +
    geom_point(data = fs, aes(Fold_change, yy, fill = Sig), shape = 23,
               size = fig_diamond_mm, stroke = ln, colour = "black") +
    scale_x_log10(breaks = fig_x_breaks, labels = format(fig_x_breaks, drop0trailing = TRUE),
                  limits = fig_x_limits, oob = scales::squish, expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = seq_along(lab_order), labels = ylab_txt,
                       expand = expansion(add = 0.6),
                       sec.axis = if (show_p) dup_axis(name = if (p_mode == "holm") "Holm p" else NULL,
                                                       labels = fs$p_txt[match(lab_order, fs$Marker)])
                       else waiver()) +
    scale_colour_manual(values = cols, name = NULL, breaks = leg_breaks,
                        labels = unname(leg_labs[leg_breaks])) +
    scale_fill_manual(values = c(`Holm p < 0.05` = "black", n.s. = "white"),
                      name = NULL, drop = FALSE) +
    labs(x = "Fold change (perineural / distal)", y = NULL,
         title = if (fig_show_title) title,
         caption = if (fig_show_caption) paste(legend_common, legend_text)) +
    guides(fill = guide_legend(order = 1, override.aes = list(size = fig_diamond_mm)),
           colour = guide_legend(order = 2, nrow = if (compact) NULL else 1,
                                 ncol = if (compact) 3 else NULL,
                                 override.aes = list(size = fig_point_mm * 1.3))) +
    theme_void(base_size = fig_text_pt, base_family = fig_font) +
    theme(text = txt, axis.text = txt, legend.text = txt, legend.title = txt,
          plot.title = element_text(size = fig_text_pt + 1, family = fig_font, hjust = 0),
          plot.caption = element_text(size = fig_text_pt - 1, family = fig_font,
                                      colour = "grey30", hjust = 0),
          axis.title.x = element_text(size = fig_text_pt, family = fig_font,
                                      margin = margin(t = 4)),
          axis.text.y.left  = element_text(size = fig_text_pt, family = fig_font, hjust = 1,
                                           margin = margin(r = 2)),
          axis.title.y.right = element_text(size = fig_text_pt, family = fig_font, angle = 0,
                                            vjust = 1, hjust = 0, margin = margin(l = 2)),
          axis.text.y.right = element_text(size = fig_text_pt, family = fig_font, hjust = 0,
                                           margin = margin(l = 3)),
          axis.text.x = element_text(size = fig_text_pt, family = fig_font,
                                     margin = margin(t = 2)),
          axis.line.x = element_line(colour = "black", linewidth = ln),
          axis.ticks.x = element_line(colour = "black", linewidth = ln),
          axis.ticks.length.x = unit(1.2, "mm"),
          legend.position = "bottom", legend.box = "vertical",
          legend.spacing.y = unit(0.5, "mm"), legend.margin = margin(0, 0, 0, 0),
          legend.key.size = unit(3, "mm"), legend.spacing.x = unit(1, "mm"),
          plot.margin = if (compact) margin(1, 1, 1, 1, "mm") else margin(4, 4, 4, 4, "mm"),
          plot.background = element_rect(fill = "white", colour = NA))
  save_fig(p, file, height, width)
}

# ---- 5. The three figures -------------------------------------------------------
prot <- markers[markers$type == "Protein", ]

set_main <- data.frame(
  var = markers$name, type = markers$type,
  label = if (fig_forest_suffix) paste0(markers$name, ifelse(markers$type == "RNA", " (RNA)", " (protein)"))
  else markers$name,
  measure = ifelse(markers$type == "RNA", "mean dots per PanCK+ cell",
                   "mean intensity of PanCK+ cells"))
res_main <- analyse(set_main, "Cohort")
forest(res_main, set_main,
       title = "PanCK+ tumour cells: perineural vs distal (cohort)",
       legend_text = paste0("Shaded: RNAscope probes, mean dots per PanCK+ cell. ",
                            "White: protein (immunofluorescence), mean intensity of PanCK+ cells."),
       file = "Fig_Cohort_forest", height = fig_forest_height_mm, width = fig_forest_width_mm,
       text_pt = fig_forest_text_pt, p_mode = fig_forest_p, compact = TRUE)

# Separate SuperPlots (e.g. ITGB4 in tumour cells): one dot per ROI, perineural vs distal
for (m in separate_markers) {
  i  <- match(m, set_main$var)
  if (is.na(i)) { message("separate_markers: ", m, " is not a marker of the cohort figure"); next }
  ok <- roi$N_PanCK_cells >= min_cells & !is.na(roi[[m]])
  panel <- m
  df <- data.frame(panel = panel, group = as.character(roi$Region[ok]), value = roi[[m]][ok],
                   block = roi$Image_Name[ok])
  st <- res_main$stats[res_main$stats$Marker == set_main$label[i], ]
  pv <- if (separate_p == "Holm") st$p_Holm else st$p
  p_lab <- paste0(if (pv < 0.001) "p < 0.001" else paste0("p = ", trimws(formatC(pv, format = "fg", digits = 2))),
                  if (separate_p == "Holm") " (Holm)" else "")
  ylab <- if (set_main$type[i] == "RNA") "Dots per PanCK+ tumour cell\n(log scale)" else
    "Mean intensity of PanCK+ tumour cells\n(log scale)"
  g <- superplot_2group(df, c("Distal", "Perineural"), c(Distal = "Distal", Perineural = "Perineural"),
                        block_cols, block_labs, data.frame(panel = panel, p_label = p_lab), ylab, log_y = TRUE)
  save_superplot(g, file.path(out_dir, paste0("Fig_Cohort_", m, "_superplot")), 1)
  legends[[paste0("Fig_Cohort_", m, "_superplot")]] <- paste0(
    m, " in PanCK+ tumour cells, perineural vs distal ROIs. Small dots: ROIs (mean ", m,
    if (set_main$type[i] == "RNA") " dots per PanCK+ cell" else " intensity of PanCK+ cells",
    "), coloured by patient; large dots: median per patient, connected; black bars: median of patient ",
    "medians. Fold change ", st$Fold_change, " (95% CI ", st$CI95, "), linear mixed model of log2 ROI values ",
    "with patient as random intercept; p value ", if (separate_p == "Holm") "Holm-adjusted over the markers of the cohort figure" else "unadjusted", ".")
}

if (has_background) {
  set_bg <- data.frame(
    var = c("CD3_background", "CD45_background"), type = "Protein",
    label = c("CD3 (TRITC channel)", "CD45 (Cy5 channel)"),
    measure = "mean intensity of PanCK+ cells (not expressed by tumour cells)")
  res_bg <- analyse(set_bg, "Background")
  if (plot_background_figures) forest(res_bg, set_bg,
                                      title = "Channel background in PanCK+ tumour cells: perineural vs distal",
                                      legend_text = paste0("CD3 and CD45 are not expressed by tumour cells; their mean ",
                                                           "intensity in PanCK+ cells measures TRITC and Cy5 background."),
                                      file = "Fig_Supp_background_controls", height = fig_height_mm[["background"]])
  
  set_corr <- data.frame(
    var = paste0(prot$name, "_corrected"), type = "Protein",
    label = paste0(prot$name, " (protein / ", prot$background, ")"),
    measure = paste0("mean intensity / mean ", prot$background, " in the same PanCK+ cells"))
  res_corr <- analyse(set_corr, "Corrected")
  if (plot_background_figures) forest(res_corr, set_corr,
                                      title = "PanCK+ tumour cells: protein corrected for channel background",
                                      legend_text = paste0("Mean protein intensity divided by the mean of the same-channel ",
                                                           "control in the same PanCK+ cells of each ROI (TRITC / CD3, Cy5 / CD45)."),
                                      file = "Fig_Supp_forest_corrected", height = fig_height_mm[["corrected"]])
}

# ---- 6. Discovery (Visium) patient: per-ROI dot plot ------------------------------
for (b in sort(unique(roi_disc$Image_Name))) {
  v  <- discovery_marker
  s  <- roi_disc[roi_disc$Image_Name == b & !is.na(roi_disc[[v]]), ]
  pe <- s[[v]][s$Region == "Perineural"]; di <- s[[v]][s$Region == "Distal"]
  wt <- wilcox.test(pe, di, exact = FALSE)
  test <- data.frame(Block = b, Marker = v, ROIs_perineural = length(pe), ROIs_distal = length(di),
                     Median_perineural = round(median(pe), 2), Median_distal = round(median(di), 2),
                     Ratio_of_medians = round(median(pe) / median(di), 2),
                     Wilcoxon_p = signif(wt$p.value, 2))
  write.csv(test, file.path(out_dir, paste0("Discovery_", v, "_test.csv")), row.names = FALSE)
  write.csv(s[, c("Image_Name", "ROI", "Region", "N_PanCK_cells", v)],
            file.path(out_dir, paste0("Discovery_", v, "_per_ROI.csv")), row.names = FALSE)
  cat("\n== Discovery patient", b, "-", v, "(ROIs as units, Wilcoxon) ==\n"); print(test, row.names = FALSE)
  
  txt <- element_text(size = fig_text_pt, colour = "black", family = fig_font)
  ln  <- fig_line_mm / 0.75
  med <- aggregate(s[[v]], list(Region = s$Region), median); names(med)[2] <- "m"
  p_txt <- if (wt$p.value < 0.001) "p < 0.001" else
    paste0("p = ", formatC(wt$p.value, format = "fg", digits = 2))
  y_top <- max(s[[v]]) * 1.12
  unit_lab <- if (markers$type[markers$name == v] == "RNA")
    bquote(.(v) ~ "dots per PanCK"^"+" ~ "cell") else bquote(.(v) ~ "intensity, PanCK"^"+" ~ "cells")
  g <- ggplot(s, aes(Region, .data[[v]])) +
    geom_point(shape = 16, size = fig_point_mm, colour = block_cols[[b]], alpha = 0.85,
               position = position_jitter(width = 0.15, height = 0, seed = 1)) +
    geom_errorbar(data = med, aes(x = Region, ymin = m, ymax = m), inherit.aes = FALSE,
                  width = 0.45, linewidth = ln * 1.5, colour = "black") +
    annotate("segment", x = 1, xend = 2, y = y_top, yend = y_top, linewidth = ln) +
    annotate("text", x = 1.5, y = y_top, label = p_txt, vjust = -0.5,
             size = fig_text_pt / .pt, family = fig_font) +
    scale_x_discrete(labels = paste0(levels(s$Region), "\nn = ", table(s$Region))) +
    scale_y_continuous(limits = c(0, y_top * 1.1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = unit_lab, title = if (fig_show_title) paste0("Discovery patient (", b, ")")) +
    theme_classic(base_size = fig_text_pt, base_family = fig_font) +
    theme(text = txt, axis.text = txt, axis.title = txt,
          plot.title = element_text(size = fig_text_pt, family = fig_font, hjust = 0.5),
          axis.line = element_line(linewidth = ln, colour = "black"),
          axis.ticks = element_line(linewidth = ln, colour = "black"),
          axis.ticks.length = unit(1, "mm"), plot.margin = margin(3, 3, 3, 3, "mm"))
  for (f in fig_formats) {
    out <- file.path(out_dir, paste0("Fig_Discovery_", v, ".", f))
    if (f == "pdf")  ggsave(out, g, device = cairo_pdf, width = fig_disc_width_mm,
                            height = fig_disc_height_mm, units = "mm")
    if (f == "png")  ggsave(out, g, width = fig_disc_width_mm, height = fig_disc_height_mm,
                            units = "mm", dpi = fig_dpi, bg = "white")
    if (f == "tiff") ggsave(out, g, device = "tiff", width = fig_disc_width_mm,
                            height = fig_disc_height_mm, units = "mm", dpi = fig_dpi,
                            bg = "white", compression = "lzw")
  }
  legends[[paste0("Fig_Discovery_", v)]] <- paste0(
    v, " in PanCK+ tumour cells of the discovery (Visium) patient ", b, ". Each dot is one ROI ",
    "(mean ", v, " dots per PanCK+ cell); bars: medians. Perineural vs distal: Wilcoxon rank-sum ",
    "test with ROIs as units (", length(pe), " perineural, ", length(di), " distal ROIs). ",
    "This patient is not included in the cohort analysis.")
}

writeLines(paste0(names(legends), ":\n", legends, "\n"), file.path(out_dir, "Figure_legends.txt"))
message("Done. Results in: ", out_dir)

