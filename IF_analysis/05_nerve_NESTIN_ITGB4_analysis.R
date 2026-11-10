###############################################################################
# Nerve analysis: NESTIN and ITGB4 in invaded vs uninvaded nerves
#
# Inputs (paths and input folder set in 00_settings.R)
#   output/<NERVE_INPUT>/cleaned/combined/All_images_combined.csv  (02_clean_nerve.R)
#   input/thresholds.csv - one row per block: Image_Name plus NESTIN_threshold and/or
#                          ITGB4_threshold (other columns and extra blocks are ignored).
#                          A block with an empty threshold is left out of that marker's
#                          % positive analysis only (its intensity analysis is kept).
#
# Horizon gives no pixel-level intensity per nerve, so two cell-based measures are used:
#   % positive cells (primary): cells with mean intensity > the block's threshold
#   median cell intensity (threshold-free check): median of the cells' mean intensity
# Only nerves with >= min_cells cells are used.
#
# Comparison: invaded vs uninvaded nerves WITHIN the same blocks (same section, staining
# and imaging); blocks without both nerve types are left out of the statistics.
#
# Statistics
#   % positive : binomial mixed model on positive / negative cell counts per nerve,
#                glmer(cbind(pos, neg) ~ Invasion + (1 | Patient_ID / Image_Name) + (1 | Nerve))
#                (if every patient has one block in the comparison: (1 | Image_Name) instead,
#                 because patient and block are then the same grouping)
#                (nerve-level random effect for overdispersion; nerves with more cells weigh more)
#                -> odds ratio invaded / uninvaded with 95% CI, and model-estimated % positive
#   intensity  : lmer(log2 median intensity ~ Invasion + (1 | Patient_ID / Image_Name))
#                -> fold change invaded / uninvaded with 95% CI (Satterthwaite df)
#   Holm correction over the markers within each measure.
#   Fit check: every model is refitted with a second optimizer (Estimate_other_optimizer);
#   lme4 convergence notes are written to Statistics.csv (Fit_note) instead of the console.
#   Random slope (robustness, patient-level question): same models with the invasion effect
#                allowed to differ between patients (1 + Invasion | block or patient);
#                columns Random_slope_* in Statistics.csv
#   Block level (robustness): median per block and nerve type, number of blocks in the
#   same direction, and a paired t-test over blocks (one value per block).
#
# Figure (SuperPlot, same format as the QuPath tenascin / p75NTR figures): small dots = nerves,
#   coloured by patient (block); large dots = block medians connected per block;
#   black bar = median of block medians; p value (Holm).
# Output folder (output/Nerve_NESTIN_ITGB4):
#   Per_nerve_values.csv, Statistics.csv, Per_block_values.csv,
#   Fig_Nerve_NESTIN_ITGB4_pct_positive(_<marker>), Fig_Supp_Nerve_NESTIN_ITGB4_intensity(_<marker>)
#   (pdf/png; one file per marker plus both markers side by side),
#   Figure_legends.txt
###############################################################################
library(ggplot2)
library(lmerTest)   # lme4 + Satterthwaite degrees of freedom

# ---- Settings ---------------------------------------------------------------
source(file.path(if (basename(getwd()) == "scripts") "." else "scripts", "00_settings.R"))
data_file      <- cleaned_file(NERVE_INPUT)
threshold_file <- THRESHOLD_FILE
out_dir        <- file.path("output", "Nerve_NESTIN_ITGB4")
min_cells      <- 10

markers <- data.frame(
  name      = c("NESTIN",                                  "ITGB4"),
  column    = c("Cells/Mean Intensity (NESTIN - Cy5 (1))", "Cells/Mean Intensity (ITGB4 - TRITC (1))"),
  threshold = c("NESTIN_threshold",                        "ITGB4_threshold"),
  stringsAsFactors = FALSE)

# Nerve groups as labelled in the figure (all nerves lie in tumour-bearing sections)
group_labels <- c(Uninvaded = "NI-TB", Invaded = "InvN")

# Figure format: superplot_2group() / save_superplot() in 00_settings.R
# (same size and style as the QuPath tenascin / p75NTR figures: 4 x 4.3 inch per panel)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Read data and thresholds ----------------------------------------------
d   <- read.csv(data_file, check.names = FALSE)
thr <- read.csv(threshold_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")

missing_cols <- setdiff(c("Image_Name", "Patient_ID", "ROI", "Type", markers$column), names(d))
if (length(missing_cols)) stop("Columns not in data: ", paste(missing_cols, collapse = "; "))
if (!"Image_Name" %in% names(thr)) stop(threshold_file, " has no column Image_Name")
for (tc in setdiff(markers$threshold, names(thr))) {
  message("No column ", tc, " in ", threshold_file, " -> % positive not analysed for that marker")
  thr[[tc]] <- NA_real_
}
dup <- thr$Image_Name[duplicated(thr$Image_Name)]
if (length(dup)) stop("Image_Name listed more than once in ", threshold_file, ": ", paste(dup, collapse = ", "))

fix <- d$Image_Name %in% names(PATIENT_OF_BLOCK)
d$Patient_ID[fix] <- PATIENT_OF_BLOCK[d$Image_Name[fix]]
d <- d[d$Type %in% c("Invaded", "Uninvaded"), ]
data_images <- sort(unique(d$Image_Name))
no_row <- setdiff(data_images, thr$Image_Name)
if (length(no_row)) message("No row in ", threshold_file, " for: ", paste(no_row, collapse = ", "),
                            " -> % positive not analysed for these blocks")

blk <- unique(d[, c("Image_Name", "Patient_ID")])
block_cols <- block_colours(blk$Image_Name, blk$Patient_ID)   # 00_settings.R
block_labs <- block_labels(blk$Image_Name, blk$Patient_ID)

# ---- 2. Per-nerve values -------------------------------------------------------
i <- match(d$Image_Name, thr$Image_Name)
pn <- do.call(rbind, lapply(split(d, list(d$Image_Name, d$ROI), drop = TRUE), function(x) {
  out <- data.frame(Patient_ID = x$Patient_ID[1], Image_Name = x$Image_Name[1], ROI = x$ROI[1],
                    Invasion = x$Type[1], N_cells = nrow(x))
  for (k in seq_len(nrow(markers))) {
    m <- markers$name[k]; v <- x[[markers$column[k]]]
    t0 <- thr[[markers$threshold[k]]][match(x$Image_Name[1], thr$Image_Name)]
    ok <- !is.na(v)
    out[[paste0(m, "_cells_measured")]] <- sum(ok)
    out[[paste0(m, "_pos_cells")]]      <- if (length(t0) && !is.na(t0)) sum(v[ok] > t0) else NA
    out[[paste0("Pct_", m, "_pos")]]    <- if (length(t0) && !is.na(t0) && sum(ok))
      100 * mean(v[ok] > t0) else NA
    out[[paste0("Median_", m, "_intensity")]] <- if (sum(ok)) median(v[ok]) else NA
  }
  out
}))
pn <- pn[order(pn$Image_Name, pn$Invasion, pn$ROI), ]
write.csv(pn, file.path(out_dir, "Per_nerve_values.csv"), row.names = FALSE)

cat("\nNerves excluded (<", min_cells, "cells):", sum(pn$N_cells < min_cells), "\n")
pn <- pn[pn$N_cells >= min_cells, ]
pn$Invasion <- factor(pn$Invasion, levels = c("Uninvaded", "Invaded"))
pn$Nerve    <- paste(pn$Image_Name, pn$ROI)

both <- tapply(pn$Invasion, pn$Image_Name, function(v) all(c("Uninvaded", "Invaded") %in% v))
one_type <- names(both)[!both]
if (length(one_type)) message("Not in the statistics (only one nerve type): ", paste(one_type, collapse = ", "))
st <- pn[pn$Image_Name %in% names(both)[both], ]
cat("Nerves used per block:\n"); print(table(st$Image_Name, st$Invasion))

# ---- 3. Statistics -------------------------------------------------------------
# Random effects: patient and block nested; if every patient contributes only one block,
# patient and block are the same grouping and only the block intercept is used
# (otherwise the model is unidentifiable and lme4 warns about convergence)
re_term <- function(s) {
  if (max(tapply(s$Image_Name, s$Patient_ID, function(x) length(unique(x)))) > 1)
    "(1 | Patient_ID / Image_Name)" else "(1 | Image_Name)"
}
# Random-slope version (robustness, patient-level question): the invasion effect may differ
# between patients/blocks, so the test asks whether the AVERAGE effect differs from 1
# relative to that between-patient variation (closer to "does it hold for patients in general")
rs_term <- function(s) {
  if (max(tapply(s$Image_Name, s$Patient_ID, function(x) length(unique(x)))) > 1)
    "(1 + Invasion | Patient_ID) + (1 | Patient_ID:Image_Name)" else "(1 + Invasion | Image_Name)"
}
# Fit and keep lme4's warnings as a note (reported in Statistics.csv instead of printed)
fit_quiet <- function(expr) {
  notes <- character(0)
  fit <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") })
  conv <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
  list(fit = fit, note = paste(unique(c(notes, conv)), collapse = " | "))
}

fit_pct <- function(m) {
  s <- st[!is.na(st[[paste0(m, "_pos_cells")]]), ]
  if (length(unique(s$Image_Name)) < 2) return(NULL)
  s$pos <- s[[paste0(m, "_pos_cells")]]; s$neg <- s[[paste0(m, "_cells_measured")]] - s$pos
  form <- as.formula(paste("cbind(pos, neg) ~ Invasion +", re_term(s), "+ (1 | Nerve)"))
  r  <- fit_quiet(glmer(form, data = s, family = binomial, control = glmerControl(optimizer = "bobyqa")))
  r2 <- fit_quiet(glmer(form, data = s, family = binomial, control = glmerControl(optimizer = "Nelder_Mead")))
  g  <- r$fit
  cf <- summary(g)$coefficients
  b  <- cf["InvasionInvaded", "Estimate"]; se <- cf["InvasionInvaded", "Std. Error"]
  rs <- fit_quiet(glmer(as.formula(paste("cbind(pos, neg) ~ Invasion +", rs_term(s), "+ (1 | Nerve)")),
                        data = s, family = binomial, control = glmerControl(optimizer = "bobyqa")))
  cr <- summary(rs$fit)$coefficients["InvasionInvaded", ]
  list(data = s, est = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
       p = cf["InvasionInvaded", "Pr(>|z|)"],
       rs_est = exp(cr[["Estimate"]]), rs_lo = exp(cr[["Estimate"]] - 1.96 * cr[["Std. Error"]]),
       rs_hi = exp(cr[["Estimate"]] + 1.96 * cr[["Std. Error"]]), rs_p = cr[["Pr(>|z|)"]],
       rs_note = paste(c(if (isSingular(rs$fit)) "singular fit", rs$note)[nzchar(c(if (isSingular(rs$fit)) "singular fit", rs$note))], collapse = " | "),
       pct_un = 100 * plogis(cf["(Intercept)", "Estimate"]),
       pct_inv = 100 * plogis(cf["(Intercept)", "Estimate"] + b),
       singular = isSingular(g), re = re_term(s), note = r$note,
       est_alt = exp(fixef(r2$fit)[["InvasionInvaded"]]))
}
fit_int <- function(m) {
  s <- st[!is.na(st[[paste0("Median_", m, "_intensity")]]) & st[[paste0("Median_", m, "_intensity")]] > 0, ]
  if (length(unique(s$Image_Name)) < 2) return(NULL)
  s$y <- log2(s[[paste0("Median_", m, "_intensity")]])
  form <- as.formula(paste("y ~ Invasion +", re_term(s)))
  r  <- fit_quiet(lmer(form, data = s, REML = TRUE))
  r2 <- fit_quiet(lmer(form, data = s, REML = TRUE, control = lmerControl(optimizer = "Nelder_Mead")))
  f  <- r$fit
  cf <- summary(f)$coefficients["InvasionInvaded", ]
  ci <- cf[["Estimate"]] + c(-1, 1) * qt(0.975, cf[["df"]]) * cf[["Std. Error"]]
  rs <- fit_quiet(lmer(as.formula(paste("y ~ Invasion +", rs_term(s))), data = s, REML = TRUE))
  cr <- summary(rs$fit)$coefficients["InvasionInvaded", ]
  cir <- cr[["Estimate"]] + c(-1, 1) * qt(0.975, cr[["df"]]) * cr[["Std. Error"]]
  list(data = s, est = 2^cf[["Estimate"]], lo = 2^ci[1], hi = 2^ci[2], p = cf[["Pr(>|t|)"]],
       rs_est = 2^cr[["Estimate"]], rs_lo = 2^cir[1], rs_hi = 2^cir[2], rs_p = cr[["Pr(>|t|)"]],
       rs_note = paste(c(if (isSingular(rs$fit)) "singular fit", rs$note)[nzchar(c(if (isSingular(rs$fit)) "singular fit", rs$note))], collapse = " | "),
       df = cf[["df"]], singular = isSingular(f), re = re_term(s), note = r$note,
       est_alt = 2^fixef(r2$fit)[["InvasionInvaded"]])
}
block_level <- function(s, col, log_scale) {
  med <- tapply(s[[col]], list(s$Image_Name, s$Invasion), median)
  med <- med[complete.cases(med), , drop = FALSE]
  dd  <- if (log_scale) log2(med[, "Invaded"]) - log2(med[, "Uninvaded"]) else med[, "Invaded"] - med[, "Uninvaded"]
  list(med = med, diff = dd,
       p = if (length(dd) >= 2 && sd(dd) > 0) t.test(dd)$p.value else NA)
}

stats <- NULL; per_block <- NULL; fits <- list()
for (m in markers$name) {
  for (measure in c("pct", "int")) {
    f <- if (measure == "pct") fit_pct(m) else fit_int(m)
    if (is.null(f)) { message(m, " (", measure, "): fewer than 2 blocks with data - skipped"); next }
    col <- if (measure == "pct") paste0("Pct_", m, "_pos") else paste0("Median_", m, "_intensity")
    bl  <- block_level(f$data, col, log_scale = measure == "int")
    fits[[paste(m, measure)]] <- c(f, list(col = col, marker = m, measure = measure, block = bl))
    per_block <- rbind(per_block, data.frame(
      Marker = m, Measure = if (measure == "pct") "% positive cells" else "median cell intensity",
      Image_Name = rownames(bl$med), Patient = block_labs[rownames(bl$med)],
      Median_uninvaded = signif(bl$med[, "Uninvaded"], 3), Median_invaded = signif(bl$med[, "Invaded"], 3),
      row.names = NULL))
    stats <- rbind(stats, data.frame(
      Marker = m,
      Measure = if (measure == "pct") "% positive cells" else "median cell intensity",
      Effect = if (measure == "pct") "odds ratio (invaded / uninvaded)" else "fold change (invaded / uninvaded)",
      Estimate = round(f$est, 2), CI95 = paste(round(f$lo, 2), "to", round(f$hi, 2)),
      p = signif(f$p, 2),
      Model_pct_uninvaded = if (measure == "pct") round(f$pct_un, 1) else NA,
      Model_pct_invaded   = if (measure == "pct") round(f$pct_inv, 1) else NA,
      Nerves_uninvaded = sum(f$data$Invasion == "Uninvaded"), Nerves_invaded = sum(f$data$Invasion == "Invaded"),
      Blocks = length(unique(f$data$Image_Name)), Patients = length(unique(f$data$Patient_ID)),
      Blocks_higher_in_invaded = paste0(sum(bl$diff > 0), "/", length(bl$diff)),
      Block_level_paired_t_p = signif(bl$p, 2),
      Random_slope_estimate = round(f$rs_est, 2),
      Random_slope_CI95 = paste(round(f$rs_lo, 2), "to", round(f$rs_hi, 2)),
      Random_slope_p = signif(f$rs_p, 2),
      Random_slope_note = if (nzchar(f$rs_note)) f$rs_note else "none",
      Random_effects = f$re, Singular_fit = f$singular,
      Estimate_other_optimizer = round(f$est_alt, 2),
      Fit_note = if (nzchar(f$note)) f$note else "none", row.names = NULL))
  }
}
stats$p_Holm <- signif(ave(stats$p, stats$Measure, FUN = function(p) p.adjust(p, "holm")), 2)
stats <- stats[, c(1:6, ncol(stats), 7:(ncol(stats) - 1))]
write.csv(stats, file.path(out_dir, "Statistics.csv"), row.names = FALSE)
if (any(stats$Fit_note != "none")) message("Fit notes for some models (see Fit_note in Statistics.csv; ",
                                           "compare Estimate with Estimate_other_optimizer - if they agree, the note does not affect the result)")
write.csv(per_block, file.path(out_dir, "Per_block_values.csv"), row.names = FALSE)
cat("\nInvaded vs uninvaded nerves (within blocks):\n")
print(stats[, c("Marker", "Measure", "Estimate", "CI95", "p", "p_Holm", "Model_pct_uninvaded",
                "Model_pct_invaded", "Blocks_higher_in_invaded", "Block_level_paired_t_p",
                "Random_slope_estimate", "Random_slope_p")], row.names = FALSE)

# ---- 4. Figures (SuperPlot, same format as the QuPath tenascin / p75NTR figures) ----
fmt_p <- function(p) if (p < 0.001) "p < 0.001" else paste0("p = ", trimws(formatC(p, format = "fg", digits = 2)))
plot_data <- function(keys) {
  df <- do.call(rbind, lapply(keys, function(k) {
    f <- fits[[k]]
    data.frame(panel = f$marker, group = as.character(f$data$Invasion),
               value = f$data[[f$col]], block = f$data$Image_Name)
  }))
  st <- do.call(rbind, lapply(keys, function(k) {
    f <- fits[[k]]
    msr <- if (f$measure == "pct") "% positive cells" else "median cell intensity"
    data.frame(panel = f$marker, p_label = fmt_p(stats$p_Holm[stats$Marker == f$marker & stats$Measure == msr]))
  }))
  list(df = df, stats = st)
}
make_figs <- function(keys, ylab, log_y, file) {
  if (!length(keys)) return(invisible())
  pdat <- plot_data(keys)
  for (k in keys) {                                  # one file per marker
    one <- pdat$df$panel == fits[[k]]$marker
    p <- superplot_2group(pdat$df[one, ], c("Uninvaded", "Invaded"), group_labels, block_cols, block_labs,
                          pdat$stats, ylab, log_y = log_y)
    save_superplot(p, file.path(out_dir, paste0(file, "_", fits[[k]]$marker)), 1)
  }
  p <- superplot_2group(pdat$df, c("Uninvaded", "Invaded"), group_labels, block_cols, block_labs,
                        pdat$stats, ylab, log_y = log_y)
  save_superplot(p, file.path(out_dir, file), length(keys))         # all markers side by side
}

legends <- character(0)
legend_common <- paste0(
  "Each small dot is one nerve (>= ", min_cells, " cells) in a tumour-bearing section, coloured by patient (block); ",
  "large dots: median per block, connected between nerve types; black bars: median of block medians. ",
  "Uninvaded (NI-TB) and invaded (InvN) nerves are compared within the same sections. ")

pct_keys <- intersect(paste(markers$name, "pct"), names(fits))
if (length(pct_keys)) {
  make_figs(pct_keys, "% positive cells", log_y = FALSE, "Fig_Nerve_NESTIN_ITGB4_pct_positive")
  legends["Fig_Nerve_NESTIN_ITGB4_pct_positive"] <- paste0(
    "Share of ", paste(sub(" pct", "", pct_keys), collapse = " / "), "-positive cells per nerve ",
    "(cells above a per-block intensity threshold). ", legend_common,
    "p values: binomial mixed model of positive / negative cells per nerve with patient and block as random ",
    "intercepts and a nerve-level random effect, Holm-adjusted over the markers.")
}
int_keys <- intersect(paste(markers$name, "int"), names(fits))
if (length(int_keys)) {
  make_figs(int_keys, "Median cell intensity (log scale)", log_y = TRUE, "Fig_Supp_Nerve_NESTIN_ITGB4_intensity")
  legends["Fig_Supp_Nerve_NESTIN_ITGB4_intensity"] <- paste0(
    "Median cell intensity per nerve (threshold-free measure). ", legend_common,
    "p values: linear mixed model of log2 median intensity with patient and block as random intercepts, ",
    "Holm-adjusted over the markers.")
}
writeLines(paste0(names(legends), ":\n", legends, "\n"), file.path(out_dir, "Figure_legends.txt"))
message("Done. Results in: ", out_dir)
