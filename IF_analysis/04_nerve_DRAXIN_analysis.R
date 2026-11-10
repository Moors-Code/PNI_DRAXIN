###############################################################################
# Nerve analysis: DRAXIN vs markers of the TGFBI+ (repair-like) Schwann cell state
#
# Inputs
#   data_file      - output/<NERVE_INPUT>/cleaned/combined/All_images_combined.csv (02_clean_nerve.R)
#   threshold_file - input/thresholds.csv, one row per block: Image_Name plus any of P75NGF_threshold,
#                    S100B_threshold, NESTIN_threshold (other columns, e.g. PanCK,
#                    and extra blocks are ignored). Image_Name must be written
#                    EXACTLY as in the data, e.g. "80830 bl4".
#
# Missing data are handled per marker:
#   - a block in the data without ANY row in the threshold file -> script stops
#     (almost always a naming mismatch)
#   - a block with a row but an empty threshold for one marker, or without data
#     for a marker (e.g. 2777_block6: no GAP43 dots, no NESTIN threshold) -> that
#     block is left out of that marker only, and this is reported
#
# Per nerve (ROI, >= min_cells cells)
#   DRAXIN, GAP43 : total dots of all cells / total cell area x 100
#                   (dots per 100 um2)
#   % P75NGF+, % S100B+, % NESTIN+ : % cells with mean intensity > the block's threshold
#
# Statistics - question: within the analysed cohort, does DRAXIN go with GAP43,
#   % P75NGF+ (repair-state markers), % NESTIN+ (supplementary) or % S100B+
#   (pan-Schwann reference) across nerves?
#   lmer(log2 DRAXIN ~ marker + (1 | Patient_ID / Image_Name))
#   - patient and block as random intercepts: nerves from the same patient /
#     block are not independent (patients 80830 and 27777 have two blocks;
#     blocks whose name does not give the patient are assigned in patient_of_block)
#   - invasion status is NOT in the primary model: invasion is one of the
#     sources of variation in Schwann cell state that the question is about
#   - marker on the scale shown in the figure: log2 GAP43 (fold change of DRAXIN
#     per 2-fold GAP43), % positive cells (fold change per +10 percentage points)
#   - each marker uses all nerves that have a value for it; Satterthwaite df;
#     Holm correction over the 4 markers
#   Robustness (statistics table): adjusted for invasion status; uninvaded nerves
#   only; per-block Spearman rho and the number of blocks with the same direction
#   Zero dot densities: half the smallest non-zero value is added before log2
#
# Paths, input folder and patient colours: 00_settings.R
# Output folder (output/Nerve_final): Per_nerve_values.csv, Statistics.csv, Per_block_correlations.csv,
#   Fig_Nerve_DRAXIN (main: GAP43, P75NGF, S100B),
#   Fig_Supp_Nerve_DRAXIN_NESTIN (pdf/png/tiff), Figure_legend.txt
###############################################################################
library(ggplot2)
library(patchwork)
library(lmerTest)   # lme4 + Satterthwaite degrees of freedom

# ---- Settings ---------------------------------------------------------------
source(file.path(if (basename(getwd()) == "scripts") "." else "scripts", "00_settings.R"))
data_file      <- cleaned_file(NERVE_INPUT)
threshold_file <- THRESHOLD_FILE
out_dir        <- file.path("output", "Nerve_final")
min_cells      <- 10

col_area   <- "Cells/Area in um2"
col_draxin <- "Cells/Dot count (Dots T3 DRAXIN - Cy5)"
col_gap43  <- "Cells/Dot count (Dots T4 GAP43 - Cy7)"
col_p75    <- "Cells/Mean Intensity (P75NGF - TRITC (1))"
col_s100   <- "Cells/Mean Intensity (S100b - Cy5 (1))"
col_nestin <- "Cells/Mean Intensity (NESTIN - Cy5 (1))"

# Patient of blocks such as "2777_block6": PATIENT_OF_BLOCK in 00_settings.R
patient_of_block <- PATIENT_OF_BLOCK
# Block colours: 00_settings.R (patient colour; a second block of the same patient
# is a darker shade, e.g. 80830 bl3, 2777_block6)

# ---- Figure settings (edit here) ---------------------------------------------
# Sizes in mm; 1 column ~ 85-89 mm, 1.5 columns ~ 120-140 mm, 2 columns ~ 180-183 mm
fig_font       <- "Arial"
fig_text_pt    <- 12              # all text, in pt
fig_width_mm   <- 300
fig_height_mm  <- 150
fig_dpi        <- 600
fig_formats    <- c("pdf", "png", "tiff")
fig_point_mm   <- 1.2            # nerve circles
fig_alpha      <- 0.8
fig_line_mm    <- 0.4            # model line and axis lines
fig_show_stats <- TRUE           # fold change and p value in each panel
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Read data and thresholds, check that everything matches --------------
d   <- read.csv(data_file, check.names = FALSE)
thr <- read.csv(threshold_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")

missing_cols <- setdiff(c("Image_Name", "Patient_ID", "ROI", "Type", col_area, col_draxin,
                          col_gap43, col_p75, col_s100, col_nestin), names(d))
if (length(missing_cols)) stop("Columns not in data: ", paste(missing_cols, collapse = "; "))
thr_cols <- c(P75NGF = "P75NGF_threshold", S100B = "S100B_threshold", NESTIN = "NESTIN_threshold")
missing_thr_cols <- setdiff(c("Image_Name", thr_cols), names(thr))
if (length(missing_thr_cols)) stop(threshold_file, " has no column(s): ",
                                   paste(missing_thr_cols, collapse = ", "))
dup <- thr$Image_Name[duplicated(thr$Image_Name)]
if (length(dup)) stop("Image_Name listed more than once in ", threshold_file, ": ",
                      paste(dup, collapse = ", "))

data_images <- sort(unique(d$Image_Name))
no_row <- setdiff(data_images, thr$Image_Name)
if (length(no_row)) stop(
  "\nNo row in ", threshold_file, " for: ", paste0('"', no_row, '"', collapse = ", "),
  "\nImage names in the data:      ", paste0('"', data_images, '"', collapse = ", "),
  "\nImage names in ", threshold_file, ": ", paste0('"', thr$Image_Name, '"', collapse = ", "),
  "\nWrite the names in ", threshold_file, " exactly as in the data.")
unused <- setdiff(thr$Image_Name, data_images)
if (length(unused)) message("Blocks in ", threshold_file, " not in this data (ignored): ",
                            paste(unused, collapse = ", "))

fix <- d$Image_Name %in% names(patient_of_block)
d$Patient_ID[fix] <- patient_of_block[d$Image_Name[fix]]
blk <- unique(d[, c("Image_Name", "Patient_ID")])
block_cols <- block_colours(blk$Image_Name, blk$Patient_ID)   # 00_settings.R
block_labs <- block_labels(blk$Image_Name, blk$Patient_ID)
cat("\nBlocks per patient:\n")
print(unique(d[order(d$Patient_ID, d$Image_Name), c("Patient_ID", "Image_Name")]), row.names = FALSE)

thr <- thr[thr$Image_Name %in% data_images, c("Image_Name", thr_cols)]
cat("\nThresholds used (NA = marker not analysed for that block):\n")
print(thr, row.names = FALSE)

d <- d[d$Type %in% c("Invaded", "Uninvaded"), ]
i <- match(d$Image_Name, thr$Image_Name)
d$P75_pos  <- d[[col_p75]]    > thr$P75NGF_threshold[i]   # NA if no threshold
d$S100_pos <- d[[col_s100]]   > thr$S100B_threshold[i]
d$NES_pos  <- d[[col_nestin]] > thr$NESTIN_threshold[i]

# ---- 2. Per-nerve values -----------------------------------------------------
dens <- function(x, col) {                     # dots per 100 um2; NA if no dot data
  ok <- !is.na(x[[col]])
  if (any(ok)) 100 * sum(x[[col]][ok]) / sum(x[[col_area]][ok]) else NA
}
pct <- function(v) if (any(!is.na(v))) 100 * mean(v, na.rm = TRUE) else NA

pn <- do.call(rbind, lapply(split(d, list(d$Image_Name, d$ROI), drop = TRUE), function(x) {
  data.frame(Patient_ID = x$Patient_ID[1], Image_Name = x$Image_Name[1], ROI = x$ROI[1],
             Invasion = x$Type[1], N_cells = nrow(x),
             Total_cell_area_um2 = sum(x[[col_area]]),
             DRAXIN_per_100um2 = dens(x, col_draxin),
             GAP43_per_100um2  = dens(x, col_gap43),
             Pct_P75NGF_pos    = pct(x$P75_pos),
             Pct_S100B_pos     = pct(x$S100_pos),
             Pct_NESTIN_pos    = pct(x$NES_pos),
             Mean_NESTIN_intensity = mean(x[[col_nestin]], na.rm = TRUE))
}))
pn <- pn[order(pn$Image_Name, pn$ROI), ]
pn$Mean_NESTIN_intensity[is.nan(pn$Mean_NESTIN_intensity)] <- NA
write.csv(pn, file.path(out_dir, "Per_nerve_values.csv"), row.names = FALSE)

cat("\nNerves excluded (<", min_cells, "cells):", sum(pn$N_cells < min_cells), "\n")
pn <- pn[pn$N_cells >= min_cells & !is.na(pn$DRAXIN_per_100um2), ]
cat("Nerves used per block:\n"); print(table(pn$Image_Name, pn$Invasion))

# log2 with a small constant for zero densities
pc_dr <- min(pn$DRAXIN_per_100um2[pn$DRAXIN_per_100um2 > 0]) / 2
pc_ga <- min(pn$GAP43_per_100um2[pn$GAP43_per_100um2 > 0], na.rm = TRUE) / 2
pn$DRAXIN_plot <- pn$DRAXIN_per_100um2 + pc_dr
pn$GAP43_plot  <- pn$GAP43_per_100um2 + pc_ga
pn$y         <- log2(pn$DRAXIN_plot)
pn$x_GAP43   <- log2(pn$GAP43_plot)
pn$x_P75NGF  <- pn$Pct_P75NGF_pos / 10
pn$x_S100B   <- pn$Pct_S100B_pos / 10
pn$x_NESTIN  <- pn$Pct_NESTIN_pos / 10
pn$Invasion  <- factor(pn$Invasion, levels = c("Uninvaded", "Invaded"))

# ---- 3. Statistics -------------------------------------------------------------
markers <- data.frame(
  x     = c("x_GAP43", "x_P75NGF", "x_S100B", "x_NESTIN"),
  name  = c("GAP43", "P75NGF", "S100B", "NESTIN"),
  raw   = c("GAP43_per_100um2", "Pct_P75NGF_pos", "Pct_S100B_pos", "Pct_NESTIN_pos"),
  unit  = c("per 2-fold GAP43", "per +10% P75NGF+", "per +10% S100B+", "per +10% NESTIN+"))

cat("\nBlocks per marker (blocks without a value for a marker are left out of it):\n")
for (k in seq_len(nrow(markers))) {
  has <- tapply(!is.na(pn[[markers$raw[k]]]), pn$Image_Name, any)
  cat(sprintf("  %-7s used: %s%s\n", markers$name[k], paste(names(has)[has], collapse = ", "),
              if (any(!has)) paste0("  | left out: ", paste(names(has)[!has], collapse = ", ")) else ""))
}

fit <- function(df, x, invasion = FALSE) {
  df <- df[!is.na(df[[x]]), ]
  form <- as.formula(paste("y ~", x, if (invasion) "+ Invasion", "+ (1 | Patient_ID / Image_Name)"))
  m  <- suppressMessages(lmer(form, data = df, REML = TRUE))
  cf <- summary(m)$coefficients
  ci <- cf[x, "Estimate"] + c(-1, 1) * qt(0.975, cf[x, "df"]) * cf[x, "Std. Error"]
  list(model = m, data = df, b = cf[x, "Estimate"], lo = ci[1], hi = ci[2],
       p = cf[x, "Pr(>|t|)"], df = cf[x, "df"])
}

per_block <- do.call(rbind, lapply(seq_len(nrow(markers)), function(k) {
  s0 <- pn[!is.na(pn[[markers$raw[k]]]), ]
  do.call(rbind, lapply(sort(unique(s0$Image_Name)), function(b) {
    s <- s0[s0$Image_Name == b, ]
    ct <- suppressWarnings(cor.test(s$DRAXIN_per_100um2, s[[markers$raw[k]]],
                                    method = "spearman", exact = FALSE))
    data.frame(Marker = markers$name[k], Image_Name = b, N_nerves = nrow(s),
               Spearman_rho = round(ct$estimate, 2), p = signif(ct$p.value, 2))
  }))
}))
write.csv(per_block, file.path(out_dir, "Per_block_correlations.csv"), row.names = FALSE)

fits <- lapply(markers$x, function(x) fit(pn, x))           # primary: no invasion term
stats <- do.call(rbind, lapply(seq_len(nrow(markers)), function(k) {
  f  <- fits[[k]]
  u  <- fit(pn[pn$Invasion == "Uninvaded", ], markers$x[k])
  ai <- fit(pn, markers$x[k], invasion = TRUE)
  rho <- per_block$Spearman_rho[per_block$Marker == markers$name[k]]
  data.frame(
    Marker = markers$name[k], Unit = markers$unit[k],
    Nerves = nrow(f$data), Blocks = length(unique(f$data$Image_Name)),
    Patients = length(unique(f$data$Patient_ID)),
    DRAXIN_fold_change = round(2^f$b, 3),
    CI95 = paste(round(2^f$lo, 3), "to", round(2^f$hi, 3)),
    df = round(f$df, 1), p = signif(f$p, 2),
    Blocks_same_direction = paste0(sum(sign(rho) == sign(f$b), na.rm = TRUE), "/", length(rho)),
    Adjusted_for_invasion_fold_change = round(2^ai$b, 3), Adjusted_for_invasion_p = signif(ai$p, 2),
    Uninvaded_only_fold_change = round(2^u$b, 3), Uninvaded_only_p = signif(u$p, 2))
}))
stats$p_Holm <- signif(p.adjust(stats$p, "holm"), 2)
stats <- stats[, c(1:9, ncol(stats), 10:(ncol(stats) - 1))]
write.csv(stats, file.path(out_dir, "Statistics.csv"), row.names = FALSE)
cat("\nDRAXIN vs markers (mixed model, patient and block as random intercepts):\n")
print(stats[, c("Marker", "Nerves", "Blocks", "Patients", "DRAXIN_fold_change", "CI95", "p",
                "p_Holm", "Blocks_same_direction", "Adjusted_for_invasion_p", "Uninvaded_only_p")],
      row.names = FALSE)

# ---- 4. Figure ----------------------------------------------------------------------
fmt_p <- function(p) if (p < 0.001) "p < 0.001" else
  paste0("p = ", trimws(formatC(p, format = "fg", digits = 2)))
txt <- element_text(size = fig_text_pt, colour = "black", family = fig_font)
ln  <- fig_line_mm / 0.75
all_blocks <- unique(pn$Image_Name)
all_blocks <- all_blocks[order(block_labs[all_blocks])]       # legend ordered by patient

panel <- function(k, xlab, log_x) {
  f   <- fits[[k]]; xv <- markers$x[k]; df <- f$data
  a   <- fixef(f$model)[["(Intercept)"]]
  xs  <- seq(min(df[[xv]]), max(df[[xv]]), length.out = 100)
  line <- data.frame(xx = if (log_x) 2^xs else xs * 10, yy = 2^(a + f$b * xs))
  pdat <- data.frame(xx = if (log_x) df$GAP43_plot else df[[markers$raw[k]]],
                     yy = df$DRAXIN_plot, Block = df$Image_Name)
  lab <- paste0(formatC(2^f$b, format = "f", digits = 2), "-fold DRAXIN ", markers$unit[k],
                "\n", fmt_p(stats$p_Holm[k]), " (Holm-adjusted)")
  g <- ggplot(pdat, aes(xx, yy)) +
    geom_point(aes(colour = Block), shape = 16, size = fig_point_mm, alpha = fig_alpha) +
    geom_line(data = line, linewidth = ln * 1.3, colour = "black") +
    scale_y_log10(labels = function(v) format(v, drop0trailing = TRUE, scientific = FALSE)) +
    scale_colour_manual(values = block_cols[all_blocks], limits = all_blocks,
                        labels = block_labs[all_blocks], name = "Patient (block)") +
    labs(x = xlab, y = if (k == 1) expression("DRAXIN dots per 100" ~ mu * m^2) else NULL,
         title = if (fig_show_stats) lab) +
    theme_classic(base_size = fig_text_pt, base_family = fig_font) +
    theme(text = txt, axis.text = txt, axis.title = txt, legend.text = txt,
          legend.title = txt, axis.line = element_line(linewidth = ln, colour = "black"),
          axis.ticks = element_line(linewidth = ln, colour = "black"),
          axis.ticks.length = unit(1, "mm"), legend.key.size = unit(3, "mm"),
          plot.title = element_text(size = fig_text_pt, family = fig_font, hjust = 0,
                                    lineheight = 1, margin = margin(b = 2)))
  g <- g + if (log_x) scale_x_log10(labels = function(v) format(v, drop0trailing = TRUE,
                                                                scientific = FALSE))
  else scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25))
  g
}

fig <- (panel(1, expression("GAP43 dots per 100" ~ mu * m^2), TRUE) +
          panel(2, expression("% P75NGF"^"+" ~ "cells"), FALSE) +
          panel(3, expression("% S100B"^"+" ~ "cells"), FALSE)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.margin = margin(2, 3, 2, 3, "mm")) &
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = fig_point_mm * 1.4,
                                                             alpha = 1)))

save_fig <- function(p, file, width) for (f in fig_formats) {
  out <- file.path(out_dir, paste0(file, ".", f))
  if (f == "pdf")  ggsave(out, p, device = cairo_pdf, width = width,
                          height = fig_height_mm, units = "mm")
  if (f == "png")  ggsave(out, p, width = width, height = fig_height_mm,
                          units = "mm", dpi = fig_dpi, bg = "white")
  if (f == "tiff") ggsave(out, p, device = "tiff", width = width,
                          height = fig_height_mm, units = "mm", dpi = fig_dpi,
                          bg = "white", compression = "lzw")
}
save_fig(fig, "Fig_Nerve_DRAXIN", fig_width_mm)

supp <- panel(4, expression("% NESTIN"^"+" ~ "cells"), FALSE) +
  labs(y = expression("DRAXIN dots per 100" ~ mu * m^2)) +
  theme(legend.position = "bottom", plot.margin = margin(2, 3, 2, 3, "mm")) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = fig_point_mm * 1.4, alpha = 1)))
save_fig(supp, "Fig_Supp_Nerve_DRAXIN_NESTIN", fig_width_mm / 2)

n_txt <- paste(sprintf("%s: %d nerves, %d blocks, %d patients", stats$Marker, stats$Nerves,
                       stats$Blocks, stats$Patients), collapse = "; ")
writeLines(paste0(
  "DRAXIN transcript density versus markers of the TGFBI+ Schwann cell state. ",
  "Each circle is one nerve (>= ", min_cells, " cells), coloured by patient (block) (", n_txt, "). ",
  "DRAXIN and GAP43: total dots divided by total cell area of the nerve. ",
  "% P75NGF+, % S100B+ (pan-Schwann reference) and % NESTIN+ (supplementary figure): share of ",
  "cells above a per-block intensity threshold. ",
  "Lines and fold changes: linear mixed model of log2 DRAXIN with the marker as fixed effect ",
  "and patient and block as random intercepts; p values Holm-adjusted over the four markers. ",
  "Results were similar after adjustment for invasion status and in uninvaded nerves only."),
  file.path(out_dir, "Figure_legend.txt"))

# ---- 5. DRAXIN in invaded vs uninvaded nerves (within tumour-bearing blocks) ----------
# Question: within the blocks that contain both nerve types, and allowing each patient /
# block its own DRAXIN level, is DRAXIN density different in invaded nerves (InvN) than in
# uninvaded nerves of the same sections (NI-TB)?
#   lmer(log2 DRAXIN ~ Invasion + (1 | Patient_ID / Image_Name))  - random intercepts,
#   same model as the NESTIN / ITGB4 intensity analysis (05); if every patient contributes
#   one block, patient and block are the same grouping -> (1 | Image_Name)
#   Random slope (robustness, patient-level question): (1 + Invasion | block) - the effect may
#   differ between patients; tests the average effect against that variation
#   Consistency: median per block and nerve type, number of blocks higher in InvN,
#   paired t-test over block medians (log2), leave-one-block-out estimates
both_dr <- tapply(pn$Invasion, pn$Image_Name, function(v) all(c("Uninvaded", "Invaded") %in% v))
inv <- pn[pn$Image_Name %in% names(both_dr)[both_dr], ]
if (any(!both_dr)) message("DRAXIN invaded vs uninvaded: not used (only one nerve type): ",
                           paste(names(both_dr)[!both_dr], collapse = ", "))
cat("\nDRAXIN invaded vs uninvaded - nerves per block:\n"); print(table(inv$Image_Name, inv$Invasion))

re_inv <- if (max(tapply(inv$Image_Name, inv$Patient_ID, function(x) length(unique(x)))) > 1)
  "(1 | Patient_ID / Image_Name)" else "(1 | Image_Name)"
fit_inv <- function(s) {
  m  <- suppressMessages(lmer(as.formula(paste("y ~ Invasion +", re_inv)), data = s, REML = TRUE))
  cf <- summary(m)$coefficients["InvasionInvaded", ]
  ci <- cf[["Estimate"]] + c(-1, 1) * qt(0.975, cf[["df"]]) * cf[["Std. Error"]]
  list(model = m, fc = 2^cf[["Estimate"]], lo = 2^ci[1], hi = 2^ci[2], p = cf[["Pr(>|t|)"]],
       df = cf[["df"]], singular = isSingular(m))
}
fi <- fit_inv(inv)
# Random slope (robustness): invasion effect allowed to differ between patients / blocks
rs_inv <- if (re_inv == "(1 | Image_Name)") "(1 + Invasion | Image_Name)" else
  "(1 + Invasion | Patient_ID) + (1 | Patient_ID:Image_Name)"
rs_warn <- character(0)
m_rs <- withCallingHandlers(suppressMessages(lmer(as.formula(paste("y ~ Invasion +", rs_inv)), data = inv, REML = TRUE)),
                            warning = function(w) { rs_warn <<- c(rs_warn, conditionMessage(w)); invokeRestart("muffleWarning") })
cr <- summary(m_rs)$coefficients["InvasionInvaded", ]
cir <- cr[["Estimate"]] + c(-1, 1) * qt(0.975, cr[["df"]]) * cr[["Std. Error"]]
rs_note <- paste(c(if (isSingular(m_rs)) "singular fit", rs_warn), collapse = " | ")
med_inv <- tapply(inv$DRAXIN_per_100um2, list(inv$Image_Name, inv$Invasion), median)
blk_lfc <- log2(med_inv[, "Invaded"] + pc_dr) - log2(med_inv[, "Uninvaded"] + pc_dr)
loo <- sapply(names(blk_lfc), function(b) fit_inv(inv[inv$Image_Name != b, ])$fc)

inv_block <- data.frame(Image_Name = rownames(med_inv), Patient = block_labs[rownames(med_inv)],
                        Nerves_uninvaded = as.vector(table(inv$Image_Name, inv$Invasion)[rownames(med_inv), "Uninvaded"]),
                        Nerves_invaded   = as.vector(table(inv$Image_Name, inv$Invasion)[rownames(med_inv), "Invaded"]),
                        Median_uninvaded = signif(med_inv[, "Uninvaded"], 3),
                        Median_invaded   = signif(med_inv[, "Invaded"], 3),
                        Ratio_invaded_uninvaded = round(2^blk_lfc, 2),
                        Fold_change_without_this_block = round(loo, 2), row.names = NULL)
inv_stats <- data.frame(
  Comparison = "DRAXIN dots per 100 um2, invaded (InvN) vs uninvaded (NI-TB)",
  Fold_change = round(fi$fc, 2), CI95 = paste(round(fi$lo, 2), "to", round(fi$hi, 2)),
  df = round(fi$df, 1), p = signif(fi$p, 2),
  Nerves_uninvaded = sum(inv$Invasion == "Uninvaded"), Nerves_invaded = sum(inv$Invasion == "Invaded"),
  Blocks = nrow(inv_block), Patients = length(unique(inv$Patient_ID)),
  Blocks_higher_in_invaded = paste0(sum(blk_lfc > 0), "/", length(blk_lfc)),
  Block_level_paired_t_p = if (length(blk_lfc) >= 2) signif(t.test(blk_lfc)$p.value, 2) else NA,
  Leave_one_block_out_range = paste(round(min(loo), 2), "to", round(max(loo), 2)),
  Random_effects = re_inv, Singular_fit = fi$singular,
  Random_slope_fold_change = round(2^cr[["Estimate"]], 2),
  Random_slope_CI95 = paste(round(2^cir[1], 2), "to", round(2^cir[2], 2)),
  Random_slope_p = signif(cr[["Pr(>|t|)"]], 2),
  Random_slope_note = if (nzchar(rs_note)) rs_note else "none")
write.csv(inv_stats, file.path(out_dir, "Statistics_DRAXIN_invaded_vs_uninvaded.csv"), row.names = FALSE)
write.csv(inv_block, file.path(out_dir, "Per_block_DRAXIN_invaded_vs_uninvaded.csv"), row.names = FALSE)
cat("\nDRAXIN invaded vs uninvaded nerves (within blocks, random intercept):\n")
print(inv_stats[, c("Fold_change", "CI95", "p", "Blocks_higher_in_invaded", "Block_level_paired_t_p",
                    "Leave_one_block_out_range", "Random_slope_fold_change", "Random_slope_CI95",
                    "Random_slope_p")], row.names = FALSE)
print(inv_block, row.names = FALSE)

# SuperPlot, same format as the NESTIN / ITGB4 figures (superplot_2group in 00_settings.R)
sp_df <- data.frame(panel = "DRAXIN", group = as.character(inv$Invasion),
                    value = inv$DRAXIN_plot, block = inv$Image_Name)
sp <- superplot_2group(sp_df, c("Uninvaded", "Invaded"), c(Uninvaded = "NI-TB", Invaded = "InvN"),
                       block_cols, block_labs, data.frame(panel = "DRAXIN", p_label = fmt_p(fi$p)),
                       expression("DRAXIN dots per 100" ~ mu * m^2 ~ "(log scale)"), log_y = TRUE)
save_superplot(sp, file.path(out_dir, "Fig_Nerve_DRAXIN_invaded_vs_uninvaded"), 1)
cat(paste0(
  "\nLegend (Fig_Nerve_DRAXIN_invaded_vs_uninvaded): DRAXIN transcript density (dots per 100 um2) in ",
  "uninvaded (NI-TB) and invaded (InvN) nerves of the same tumour-bearing sections. Each small dot is one ",
  "nerve (>= ", min_cells, " cells), coloured by patient; large dots: median per patient, connected between ",
  "nerve types; black bars: median of patient medians. p value: linear mixed model of log2 DRAXIN density ",
  "with invasion status as fixed effect and patient as random intercept (", inv_stats$Nerves_uninvaded,
  " NI-TB and ", inv_stats$Nerves_invaded, " InvN nerves, ", inv_stats$Patients, " patients). ",
  "DRAXIN was higher in InvN in ", inv_stats$Blocks_higher_in_invaded, " patients.\n"),
  file = file.path(out_dir, "Figure_legend.txt"), append = TRUE)

message("Done. Results in: ", out_dir)

