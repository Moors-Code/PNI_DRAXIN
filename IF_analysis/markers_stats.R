# =============================================================================
# Nerve markers – normalised intensity statistics (within tumour-bearing sections)
#
# Analyses (marker in region):
#   GAP43    in the nerve            (original_nerve)
#   Tenascin in the 0-50 um ring     (ring_immune_cells)
#   p75NTR   in the nerve            (original_nerve)
#   p75NTR   in the perineural ring  (ring_perineural)
#
# Normalisation (within each section / image, used for plotting across patients):
#   normalised = raw mean intensity of the region
#                / median raw intensity of that marker over the image's ring_normalization rings
#   1 = same as the section's surrounding-tissue background.
#
# Comparison: InvN vs NI-TB only. Both groups lie in the same tumour-bearing section and share
# staining and imaging, so the per-section normalisation factor cancels out of this comparison.
# NI-Adj (separate section) is kept in the per-nerve table for the descriptive supplementary plot.
#
# Model per analysis (tumour-bearing sections only), linear mixed model on log intensity (glmmTMB):
#   log(normalised) ~ group + (1 | patient)        (one tumour-bearing section per patient)
#   -> fold difference InvN / NI-TB (ratio of geometric means), Wald 95% CI, p value.
# Patient-level check: paired t-test on log patient medians (one value per patient).
#
# Run from the project folder:  Rscript scripts/markers_stats.R
#   in : input/plotting_data.csv
#   out: output/markers_per_nerve.csv            one row per nerve and analysis (all three groups)
#        output/markers_normalisation_factors.csv
#        output/markers_model_results.csv
#        output/markers_patient_ratios.csv
#        output/markers_stats_report.txt
# Packages: install.packages("glmmTMB")
# =============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(x, y) if (is.null(x)) y else x
library(glmmTMB)

if (basename(getwd()) == "scripts") setwd("..")
dir.create("output", showWarnings = FALSE)

ANALYSES <- data.frame(
  id     = c("gap43_nerve",    "tenascin_immune_ring",   "p75ntr_nerve",   "p75ntr_perineural"),
  marker = c("gap43",          "tenascin",               "p75ntr",         "p75ntr"),
  region = c("original_nerve", "ring_immune_cells",      "original_nerve", "ring_perineural"),
  label  = c("GAP43, nerve",   "Tenascin, 0-50 um ring", "p75NTR, nerve",  "p75NTR, perineural ring"))

GROUPS <- c("NI-Adj", "NI-TB", "InvN")
PATIENT_LABELS <- c("70865" = "Patient 3", "56945" = "Patient 4",
                    "80830" = "Patient 5", "23885" = "Patient 6", "27777" = "Patient 7")

# ---- 1. data & normalisation factors ----------------------------------------
raw <- read.csv("input/plotting_data.csv", check.names = FALSE)
names(raw) <- tolower(names(raw))
raw$patient  <- as.character(raw$patient)
raw$nerve_id <- sub("^(nerve_[0-9]+).*$", "\\1", raw$annotation_name)
raw$group <- with(raw, ifelse(nerve_classification == "Uninvaded_nerve" & block_type == "Adjacent",      "NI-Adj",
                       ifelse(nerve_classification == "Uninvaded_nerve" & block_type == "Tumor_bearing", "NI-TB",
                       ifelse(nerve_classification == "Invaded_nerve"   & block_type == "Tumor_bearing", "InvN", NA))))

ref_rings <- raw[raw$annotation_type == "ring_normalization", ]
markers <- unique(ANALYSES$marker)
factors <- aggregate(ref_rings[, paste0("intensity_raw_", markers)], by = list(image_name = ref_rings$image_name),
                     FUN = median, na.rm = TRUE)
names(factors) <- c("image_name", paste0("factor_", markers))
write.csv(factors, "output/markers_normalisation_factors.csv", row.names = FALSE)

nerve_area <- raw[raw$annotation_type == "original_nerve", c("image_name", "nerve_id", "area")]
names(nerve_area)[3] <- "nerve_area"

d <- do.call(rbind, lapply(seq_len(nrow(ANALYSES)), function(i) {
  a <- ANALYSES[i, ]
  x <- raw[raw$annotation_type == a$region & !is.na(raw$group), ]   # drops invaded nerves in adjacent blocks
  x <- merge(x, factors[, c("image_name", paste0("factor_", a$marker))], by = "image_name")
  x <- merge(x, nerve_area, by = c("image_name", "nerve_id"))
  data.frame(analysis = a$id, label = a$label, marker = a$marker, region = a$region,
             image_name = x$image_name, nerve_id = x$nerve_id, patient = x$patient,
             patient_label = unname(PATIENT_LABELS[x$patient]), block_type = x$block_type,
             group = x$group, nerve_area = x$nerve_area,
             raw = x[[paste0("intensity_raw_", a$marker)]],
             normalised = x[[paste0("intensity_raw_", a$marker)]] / x[[paste0("factor_", a$marker)]])
}))
d <- d[!is.na(d$normalised) & d$normalised > 0, ]
d$group <- factor(d$group, levels = GROUPS)
write.csv(d, "output/markers_per_nerve.csv", row.names = FALSE)

# ---- 2. mixed models: InvN vs NI-TB within tumour-bearing sections -----------
results <- NULL; summaries <- NULL; pat <- NULL
for (i in seq_len(nrow(ANALYSES))) {
  a  <- ANALYSES[i, ]
  da <- d[d$analysis == a$id & d$group %in% c("NI-TB", "InvN"), ]
  da$grp     <- factor(da$group, levels = c("NI-TB", "InvN"))
  da$patient <- factor(da$patient)
  da$y       <- log(da$normalised)
  m  <- glmmTMB(y ~ grp + (1 | patient), family = gaussian, data = da)
  co <- summary(m)$coefficients$cond["grpInvN", ]

  # patient level: log median per patient and group, paired across patients
  med  <- tapply(da$y, list(as.character(da$patient_label), da$grp), median)
  diff <- med[, "InvN"] - med[, "NI-TB"]
  tt   <- t.test(diff)
  pat  <- rbind(pat, data.frame(analysis = a$id, label = a$label, patient = rownames(med),
                                median_NI_TB = exp(med[, "NI-TB"]), median_InvN = exp(med[, "InvN"]),
                                ratio_InvN_vs_NI_TB = exp(diff), row.names = NULL))

  results <- rbind(results, data.frame(
    analysis = a$id, label = a$label, comparison = "InvN vs NI-TB",
    fold = exp(co["Estimate"]), ci_low = exp(co["Estimate"] - 1.96 * co["Std. Error"]),
    ci_high = exp(co["Estimate"] + 1.96 * co["Std. Error"]), p_value = co["Pr(>|z|)"],
    patients_higher = sum(diff > 0), patients_lower = sum(diff < 0), n_patients = length(diff),
    patient_level_fold = exp(mean(diff)), patient_level_p = tt$p.value,
    n_NI_TB = sum(da$grp == "NI-TB"), n_InvN = sum(da$grp == "InvN"),
    converged = isTRUE(m$sdr$pdHess), row.names = NULL))
  summaries <- c(summaries, paste("=====", a$label, "====="), capture.output(summary(m)), "")
}
write.csv(results, "output/markers_model_results.csv", row.names = FALSE)
write.csv(pat, "output/markers_patient_ratios.csv", row.names = FALSE)

# ---- 3. report ----------------------------------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.2g", p))
report <- c(
  sprintf("Nerve marker statistics (InvN vs NI-TB, tumour-bearing sections) – %s", Sys.Date()), "",
  "RESULTS (fold = InvN / NI-TB, ratio of geometric means)",
  sprintf("  %-26s mixed model: fold %.2f [%.2f–%.2f], %s | patients: %d higher, %d lower; patient-level paired t: fold %.2f, %s | nerves NI-TB %d, InvN %d",
          results$label, results$fold, results$ci_low, results$ci_high, fmt_p(results$p_value),
          results$patients_higher, results$patients_lower, results$patient_level_fold,
          fmt_p(results$patient_level_p), results$n_NI_TB, results$n_InvN), "",
  "PER-PATIENT RATIOS (InvN median / NI-TB median)",
  sprintf("  %-26s %s", pat$label, sprintf("%s: %.2f", pat$patient, pat$ratio_InvN_vs_NI_TB)), "",
  "METHODS",
  paste("Mean marker intensity was measured in the nerve, the perineural ring or the 0-50 um ring and normalised to the",
        "median intensity of the same marker in the normalisation rings of the same section. Because marker",
        "intensities varied between tissue sections, comparisons were restricted to nerves within the same",
        "tumour-bearing section (NI-TB vs InvN). Log-transformed normalised intensities were analysed with a",
        "linear mixed-effects model (glmmTMB, R) with nerve group as a fixed effect and patient as a random",
        "intercept; differences are reported as fold differences (ratio of geometric means) with Wald 95%",
        "confidence intervals. Consistency across patients was assessed with a paired t-test on log-transformed",
        "patient medians."), "",
  "MODEL SUMMARIES", summaries,
  "SESSION INFO", capture.output(sessionInfo()))
writeLines(report, "output/markers_stats_report.txt")
cat(report[1:(grep("^METHODS", report) - 1)], sep = "\n")
