# =============================================================================
# Figure 2E – statistics
#
# CD45+ density around nerves: InvN vs NI-TB and InvN vs NI-Adj, using
#   (a) the full perineural area        (ring_immune_cells)
#   (b) the non-tumour perineural area  (ring_immune_cells_no_tumor)
# The numerator is always the CD45+ count of the full ring_immune_cells.
#
# Model: negative-binomial mixed model (glmmTMB)
#   count ~ group + (1 | patient) + (1 | image) + offset(log(area))
#
# Run from the project folder:  Rscript scripts/fig2E_stats.R
#   in : input/plotting_data.csv
#   out: output/fig2E_per_nerve.csv        one row per nerve (used by fig2E_plot.R)
#        output/fig2E_model_results.csv    density ratios, 95% CI, p values (raw + Holm)
#        output/fig2E_patient_ratios.csv   InvN / reference ratio per patient
#        output/fig2E_stats_report.txt     model summaries + results text
# Packages: install.packages("glmmTMB")
# =============================================================================

# `%||%` is only in base R >= 4.4; some package versions use it, so define it for older R
if (!exists("%||%", mode = "function")) `%||%` <- function(x, y) if (is.null(x)) y else x

library(glmmTMB)

if (basename(getwd()) == "scripts") setwd("..")
dir.create("output", showWarnings = FALSE)

GROUPS <- c("NI-Adj", "NI-TB", "InvN")
EXCLUDE_PATIENTS <- "27777"                       # no usable CD45 counts
PATIENT_LABELS <- c("70865" = "Patient 3", "56945" = "Patient 4",
                    "80830" = "Patient 5", "23885" = "Patient 6", "27777" = "Patient 7")

# ---- 1. one row per nerve ---------------------------------------------------
raw <- read.csv("input/plotting_data.csv", check.names = FALSE)
names(raw) <- tolower(names(raw))
raw$nerve_id <- sub("^(nerve_[0-9]+).*$", "\\1", raw$annotation_name)

ring <- raw[raw$annotation_type == "ring_immune_cells",
            c("image_name", "nerve_id", "patient", "block_type", "nerve_classification",
              "area", "cell_count_cd45")]
no_tumour <- raw[raw$annotation_type == "ring_immune_cells_no_tumor",
                 c("image_name", "nerve_id", "area")]
names(no_tumour)[3] <- "area_tumorfree"
d <- merge(ring, no_tumour, by = c("image_name", "nerve_id"))

d$group <- with(d, ifelse(nerve_classification == "Uninvaded_nerve" & block_type == "Adjacent",      "NI-Adj",
                   ifelse(nerve_classification == "Uninvaded_nerve" & block_type == "Tumor_bearing", "NI-TB",
                   ifelse(nerve_classification == "Invaded_nerve"   & block_type == "Tumor_bearing", "InvN", NA))))
d <- subset(d, !is.na(group) &                          # drops invaded nerves in adjacent blocks
               !patient %in% EXCLUDE_PATIENTS &
               !is.na(cell_count_cd45) & cell_count_cd45 >= 0 &
               area_tumorfree > 0)

d$count             <- as.integer(d$cell_count_cd45)
d$patient_label     <- PATIENT_LABELS[as.character(d$patient)]
d$density_full      <- d$count / d$area           * 1e8   # cells / cm² (area in µm²)
d$density_tumorfree <- d$count / d$area_tumorfree * 1e8
d$group <- factor(d$group, levels = GROUPS)
write.csv(d, "output/fig2E_per_nerve.csv", row.names = FALSE)

# ---- 2. mixed models --------------------------------------------------------
d$grp     <- relevel(d$group, ref = "InvN")           # InvN = reference level
d$patient <- factor(d$patient)
d$image   <- factor(d$image_name)

versions <- c(full = "Full perineural area", tumorfree = "Non-tumour perineural area")
area_col <- c(full = "area", tumorfree = "area_tumorfree")

# All three pairwise comparisons, written as "first vs second" (ratio = first / second)
COMPARISONS <- list(c("InvN", "NI-TB"), c("InvN", "NI-Adj"), c("NI-TB", "NI-Adj"))

results <- NULL; summaries <- NULL
for (v in names(versions)) {
  d$log_area <- log(d[[area_col[v]]])
  m <- glmmTMB(count ~ grp + (1 | patient) + (1 | image) + offset(log_area),
               family = nbinom2, data = d)
  b <- fixef(m)$cond; V <- vcov(m)$cond
  for (cmp in COMPARISONS) {
    # log density of each group relative to InvN (the reference level, coefficient 0)
    w <- setNames(rep(0, length(b)), names(b))
    if (cmp[1] != "InvN") w[paste0("grp", cmp[1])] <- w[paste0("grp", cmp[1])] + 1
    if (cmp[2] != "InvN") w[paste0("grp", cmp[2])] <- w[paste0("grp", cmp[2])] - 1
    est <- sum(w * b); se <- sqrt(drop(t(w) %*% V %*% w))
    results <- rbind(results, data.frame(
      version = v, panel = versions[[v]], comparison = paste(cmp[1], "vs", cmp[2]),
      density_ratio = exp(est), ci_low = exp(est - 1.96 * se), ci_high = exp(est + 1.96 * se),
      p_value = 2 * pnorm(-abs(est / se)), converged = isTRUE(m$sdr$pdHess)))
  }
  summaries <- c(summaries, paste("=====", versions[[v]], "====="), capture.output(summary(m)), "")
}
# Holm correction for the three comparisons within each model
results$p_holm <- ave(results$p_value, results$version, FUN = function(p) p.adjust(p, method = "holm"))
write.csv(results, "output/fig2E_model_results.csv", row.names = FALSE)

# ---- 3. patient-level ratios (InvN median / reference median) --------------
pat <- NULL
for (v in names(versions)) {
  med <- tapply(d[[paste0("density_", v)]], list(d$patient_label, d$group), median)
  for (cmp in COMPARISONS)
    pat <- rbind(pat, data.frame(version = v, panel = versions[[v]], comparison = paste(cmp[1], "vs", cmp[2]),
                                 patient = rownames(med), ratio = med[, cmp[1]] / med[, cmp[2]],
                                 first_lower = med[, cmp[1]] < med[, cmp[2]], row.names = NULL))
}
write.csv(pat, "output/fig2E_patient_ratios.csv", row.names = FALSE)

# ---- 4. report --------------------------------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.2g", p))
fmt <- function(v, cmp) with(results[results$version == v & results$comparison == cmp, ],
  sprintf("ratio %.2f, 95%% CI %.2f–%.2f, Holm-adjusted %s", density_ratio, ci_low, ci_high, fmt_p(p_holm)))
n_lower <- function(v, cmp) with(pat[pat$version == v & pat$comparison == cmp, ],
  sprintf("%d/%d (%s)", sum(first_lower), length(first_lower), paste(patient[first_lower], collapse = ", ")))
cmp_names <- sapply(COMPARISONS, paste, collapse = " vs ")

report <- c(
  sprintf("Figure 2E statistics – %s – %d nerves, %d patients, %d images", Sys.Date(),
          nrow(d), nlevels(d$patient), nlevels(d$image)), "",
  "MODEL RESULTS (density ratio = first group / second group; p values Holm-adjusted within each model)",
  sprintf("  %-27s %-16s %s", results$panel, results$comparison,
          mapply(fmt, results$version, results$comparison)), "",
  "PATIENTS WHERE THE FIRST GROUP HAS THE LOWER MEDIAN DENSITY",
  unlist(lapply(names(versions), function(v)
    sprintf("  %-27s %-16s %s", versions[[v]], cmp_names, sapply(cmp_names, function(c) n_lower(v, c))))),
  "  Per-patient ratios: output/fig2E_patient_ratios.csv (a ratio of 0.98 counts as 'lower')", "",
  "METHODS",
  paste("CD45+ cell counts in the perineural area were modelled with a negative-binomial mixed-effects model",
        "(glmmTMB, R), with nerve group (NI-Adj, NI-TB, InvN) as a fixed effect, patient and image as random",
        "intercepts, and log perineural area (full, or non-tumour area only) as an offset. All pairwise group",
        "comparisons are reported as density ratios with Wald 95% confidence intervals; p values were adjusted",
        "for the three comparisons using the Holm method."), "",
  "RESULTS",
  paste0("CD45+ density in the perineural area did not differ between NI-TB and NI-Adj (",
         fmt("full", "NI-TB vs NI-Adj"), "), but was lower around InvN than around NI-TB (",
         fmt("full", "InvN vs NI-TB"), ") and NI-Adj (", fmt("full", "InvN vs NI-Adj"), ") (Figure 2E). ",
         "After excluding tumour-occupied area, density around InvN remained lower than around NI-TB (",
         fmt("tumorfree", "InvN vs NI-TB"), ") and NI-Adj (", fmt("tumorfree", "InvN vs NI-Adj"), ")."), "",
  "MODEL SUMMARIES", summaries,
  "SESSION INFO", capture.output(sessionInfo()))
writeLines(report, "output/fig2E_stats_report.txt")
cat(report[1:20], sep = "\n")

