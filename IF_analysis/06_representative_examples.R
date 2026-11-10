###############################################################################
# Suggest representative ROIs / nerves for example images, one set per finding
#
# Run after 03_cancer_analysis.R, 04_nerve_DRAXIN_analysis.R and 05_nerve_NESTIN_ITGB4_analysis.R.
# Inputs: output/Cancer_final/Per_ROI_values.csv
#         output/Nerve_final/Per_nerve_values.csv
#         output/Nerve_NESTIN_ITGB4/Per_nerve_values.csv
#
# Group findings (e.g. PFKFB3 lower in tumour cells next to nerves):
#   1. per block: effect = median log value in group A - median log value in group B
#      (blocks need >= min_per_group ROIs / nerves in each group)
#   2. typical effect = median of the block effects; the suggested block is the one in
#      the expected direction whose effect is closest to the typical effect
#   3. in that block, the ROIs / nerves closest to their group median are suggested
#      (pct = percentile within block and group; 50 = exactly typical)
# NESTIN / ITGB4 in nerves: median cell intensity per nerve (no threshold; same measure as
#   Fig_Supp_Nerve_NESTIN_ITGB4_intensity in 05)
# Correlation finding (DRAXIN with GAP43 and P75NGF):
#   1. per block: Spearman rho of DRAXIN with each marker (>= min_nerves_corr nerves)
#   2. suggested block: both rho in the expected direction, largest smaller rho
#   3. in that block, z-scores of log DRAXIN and of the markers (mean of z log GAP43 and
#      z % P75NGF+); a "high" and a "low" nerve that lie close to the trend
#      (|z DRAXIN - z markers| small) and are far apart along it
# Output (output/Representative_examples): Representative_examples.csv (all suggestions),
#   Block_ranking.csv (every block per finding), Representative_examples.txt (summary)
# Always check the suggested ROIs in the images: the script knows numbers, not image quality.
###############################################################################
source(file.path(if (basename(getwd()) == "scripts") "." else "scripts", "00_settings.R"))

out_dir        <- file.path("output", "Representative_examples")
min_cells      <- 10      # same as the analysis scripts
min_per_group  <- 3       # ROIs / nerves per group a block needs to be suggested
min_nerves_corr <- 6      # nerves a block needs for the correlation finding
n_suggest      <- 3       # candidates listed per group
exclude_blocks <- c("70865 bl15")   # e.g. the discovery patient (shown in its own figure)

cancer_file <- file.path("output", "Cancer_final", "Per_ROI_values.csv")
draxin_file <- file.path("output", "Nerve_final", "Per_nerve_values.csv")
nestin_file <- file.path("output", "Nerve_NESTIN_ITGB4", "Per_nerve_values.csv")

# Findings: data, value column, groups (A vs B), expected direction of A vs B
group_findings <- read.table(header = TRUE, sep = "|", strip.white = TRUE, text = "
finding                                   | file   | column         | group_col | A          | B         | direction
PFKFB3 lower in tumour cells next to nerves | cancer | PFKFB3         | Region    | Perineural | Distal    | lower
ITGB4 higher in tumour cells next to nerves | cancer | ITGB4          | Region    | Perineural | Distal    | higher
TOMM20 higher in tumour cells next to nerves| cancer | TOMM20         | Region    | Perineural | Distal    | higher
NESTIN higher in invaded nerves             | nestin | Median_NESTIN_intensity | Invasion  | Invaded    | Uninvaded | higher
ITGB4 lower in invaded nerves               | nestin | Median_ITGB4_intensity  | Invasion  | Invaded    | Uninvaded | lower")
corr_finding <- list(finding = "DRAXIN goes with GAP43 and P75NGF", y = "DRAXIN_per_100um2",
                     x_log = "GAP43_per_100um2", x_lin = "Pct_P75NGF_pos", direction = "positive")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read the analysis outputs ---------------------------------------------------
read_if <- function(f, cells_col) {
  if (!file.exists(f)) { message("Not found (run the analysis script first): ", f); return(NULL) }
  x <- read.csv(f, check.names = FALSE)
  x <- x[x[[cells_col]] >= min_cells & !x$Image_Name %in% exclude_blocks, ]
  x
}
data <- list(cancer = read_if(cancer_file, "N_PanCK_cells"),
             nestin = read_if(nestin_file, "N_cells"),
             draxin = read_if(draxin_file, "N_cells"))
lg <- function(v) { pc <- if (any(v == 0, na.rm = TRUE)) min(v[v > 0], na.rm = TRUE) / 2 else 0; log2(v + pc) }
lab_of <- function(x) {
  b <- unique(x[, c("Image_Name", "Patient_ID")])
  block_labels(b$Image_Name, b$Patient_ID)
}

suggestions <- NULL; ranking <- NULL; summary_txt <- character(0)

# ---- Group findings -------------------------------------------------------------------
for (k in seq_len(nrow(group_findings))) {
  gf <- group_findings[k, ]; x <- data[[gf$file]]
  if (is.null(x) || !gf$column %in% names(x)) { message("Skipped: ", gf$finding); next }
  x <- x[!is.na(x[[gf$column]]) & x[[gf$group_col]] %in% c(gf$A, gf$B), ]
  x$lv <- lg(x[[gf$column]]); labs <- lab_of(x)
  eff <- do.call(rbind, lapply(split(x, x$Image_Name), function(s) {
    a <- s$lv[s[[gf$group_col]] == gf$A]; b <- s$lv[s[[gf$group_col]] == gf$B]
    data.frame(Image_Name = s$Image_Name[1], n_A = length(a), n_B = length(b),
               effect_log2 = if (length(a) && length(b)) median(a) - median(b) else NA)
  }))
  eff <- eff[!is.na(eff$effect_log2), ]
  typical <- median(eff$effect_log2)
  sgn <- if (gf$direction == "higher") 1 else -1
  eff$expected_direction <- sign(eff$effect_log2) == sgn
  eff$enough <- eff$n_A >= min_per_group & eff$n_B >= min_per_group
  eff$distance_to_typical <- abs(eff$effect_log2 - typical)
  eff <- eff[order(!(eff$expected_direction & eff$enough), eff$distance_to_typical), ]
  ranking <- rbind(ranking, data.frame(Finding = gf$finding, Image_Name = eff$Image_Name,
                                       Patient = labs[eff$Image_Name], ROIs_A = eff$n_A, ROIs_B = eff$n_B,
                                       Ratio_A_vs_B = round(2^eff$effect_log2, 2), Typical_ratio = round(2^typical, 2),
                                       Expected_direction = eff$expected_direction, Enough_ROIs = eff$enough, row.names = NULL))
  best <- eff$Image_Name[1]
  if (!(eff$expected_direction[1] && eff$enough[1]))
    message(gf$finding, ": no block with the expected direction and >= ", min_per_group,
            " per group; showing the closest block")
  s <- x[x$Image_Name == best, ]
  for (g in c(gf$A, gf$B)) {
    sg <- s[s[[gf$group_col]] == g, ]
    sg$dev <- abs(sg$lv - median(sg$lv))
    sg$pct <- round(100 * (rank(sg$lv) - 0.5) / nrow(sg))
    sg <- head(sg[order(sg$dev), ], n_suggest)
    suggestions <- rbind(suggestions, data.frame(
      Finding = gf$finding, Role = paste0(g, " example"), Rank = seq_len(nrow(sg)),
      Image_Name = best, Patient = labs[best], ROI = sg$ROI,
      Value = signif(sg[[gf$column]], 3), Measure = gf$column,
      Percentile_in_group = sg$pct, Block_ratio_A_vs_B = round(2^eff$effect_log2[1], 2), row.names = NULL))
  }
  top <- suggestions[suggestions$Finding == gf$finding & suggestions$Rank == 1, ]
  summary_txt <- c(summary_txt, sprintf(
    "%s\n  %s: %s/%s ratio %.2f (typical across blocks %.2f)\n  examples: %s %s (%s = %s) and %s %s (%s = %s)\n",
    gf$finding, labs[best], gf$A, gf$B, 2^eff$effect_log2[1], 2^typical,
    gf$A, top$ROI[1], gf$column, top$Value[1], gf$B, top$ROI[2], gf$column, top$Value[2]))
}

# ---- Correlation finding (DRAXIN vs GAP43 and P75NGF) ----------------------------------
x <- data$draxin
cf <- corr_finding
if (!is.null(x) && all(c(cf$y, cf$x_log, cf$x_lin) %in% names(x))) {
  x <- x[complete.cases(x[, c(cf$y, cf$x_log, cf$x_lin)]), ]
  x$ly <- lg(x[[cf$y]]); x$lx1 <- lg(x[[cf$x_log]]); x$lx2 <- x[[cf$x_lin]]
  labs <- lab_of(x)
  rk <- do.call(rbind, lapply(split(x, x$Image_Name), function(s) {
    r1 <- if (nrow(s) >= 3) suppressWarnings(cor(s$ly, s$lx1, method = "spearman")) else NA
    r2 <- if (nrow(s) >= 3) suppressWarnings(cor(s$ly, s$lx2, method = "spearman")) else NA
    data.frame(Image_Name = s$Image_Name[1], Nerves = nrow(s), rho_GAP43 = round(r1, 2),
               rho_P75NGF = round(r2, 2))
  }))
  rk$min_rho <- pmin(rk$rho_GAP43, rk$rho_P75NGF)
  rk$enough  <- rk$Nerves >= min_nerves_corr
  rk <- rk[order(!(rk$enough & rk$min_rho > 0), -rk$min_rho), ]
  ranking <- rbind(ranking, data.frame(Finding = cf$finding, Image_Name = rk$Image_Name,
                                       Patient = labs[rk$Image_Name], ROIs_A = rk$Nerves, ROIs_B = NA,
                                       Ratio_A_vs_B = NA, Typical_ratio = NA,
                                       Expected_direction = rk$min_rho > 0, Enough_ROIs = rk$enough, row.names = NULL))
  best <- rk$Image_Name[1]
  s <- x[x$Image_Name == best, ]
  z <- function(v) if (sd(v) > 0) (v - mean(v)) / sd(v) else v * 0
  s$zy <- z(s$ly); s$zx <- (z(s$lx1) + z(s$lx2)) / 2
  s$off_trend <- abs(s$zy - s$zx); s$along <- (s$zy + s$zx) / 2
  on <- s[s$off_trend <= quantile(s$off_trend, 0.5), ]           # closest half to the trend
  pick <- function(df, role) data.frame(
    Finding = cf$finding, Role = role, Rank = seq_len(nrow(df)), Image_Name = best,
    Patient = labs[best], ROI = df$ROI,
    Value = paste0("DRAXIN ", signif(df[[cf$y]], 3), "; GAP43 ", signif(df[[cf$x_log]], 3),
                   "; % P75NGF+ ", round(df[[cf$x_lin]], 1)),
    Measure = "dots per 100 um2 / % positive", Percentile_in_group = round(100 * rank(s$along)[match(df$ROI, s$ROI)] / nrow(s)),
    Block_ratio_A_vs_B = NA, row.names = NULL)
  hi <- head(on[order(-on$along), ], n_suggest); lo <- head(on[order(on$along), ], n_suggest)
  suggestions <- rbind(suggestions, pick(hi, "High DRAXIN, GAP43 and P75NGF"),
                       pick(lo, "Low DRAXIN, GAP43 and P75NGF"))
  summary_txt <- c(summary_txt, sprintf(
    "%s\n  %s: Spearman rho GAP43 %.2f, P75NGF %.2f (%d nerves)\n  high example: %s (%s)\n  low example:  %s (%s)\n",
    cf$finding, labs[best], rk$rho_GAP43[1], rk$rho_P75NGF[1], rk$Nerves[1],
    hi$ROI[1], pick(hi[1, ], "")$Value, lo$ROI[1], pick(lo[1, ], "")$Value))
} else message("Skipped: ", cf$finding)

write.csv(suggestions, file.path(out_dir, "Representative_examples.csv"), row.names = FALSE)
write.csv(ranking, file.path(out_dir, "Block_ranking.csv"), row.names = FALSE)
writeLines(c("Suggested example ROIs / nerves (rank 1 per group; alternatives in Representative_examples.csv)",
             "Check them in the images before use.", "", summary_txt),
           file.path(out_dir, "Representative_examples.txt"))
cat(summary_txt, sep = "\n")
message("Done. Results in: ", out_dir)
