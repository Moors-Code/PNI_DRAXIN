# =============================================================================
# Suggest representative nerves for example images (NI-TB and InvN)
#
# For each tumour-bearing section (= patient), nerves are scored by how close they are to the
# patient's median for their group, in every analysis at once (GAP43 nerve, tenascin 0-50 um ring,
# p75NTR nerve, p75NTR perineural). A nerve that is "typical" for all markers can be shown in one
# multi-channel image panel. Very small nerves are skipped so the example is visible.
#
# Score = mean over analyses of |log(value) - log(patient group median)| / spread of that group
#         (spread = median absolute deviation of log values). Lower = more representative.
#
# Also ranks patients by how close their InvN/NI-TB ratios are to the overall effect (for the
# analyses with a significant difference), to suggest which section to take the examples from.
#
# Run from the project folder, after markers_stats.R:  Rscript scripts/markers_representative_nerves.R
#   in : output/markers_per_nerve.csv, output/markers_model_results.csv, output/markers_patient_ratios.csv
#   out: output/representative_nerves.csv     top candidates per patient and group
#        output/representative_patients.csv   how typical each patient's effect is
# =============================================================================

if (basename(getwd()) == "scripts") setwd("..")

N_PER_GROUP    <- 5      # candidates listed per patient and group
MIN_AREA_QUANT <- 0.25   # skip nerves smaller than this quantile of nerve area (per section)
MIN_NERVES     <- 5      # a patient needs at least this many NI-TB and InvN nerves to be suggested
# Analyses the example should illustrate; NULL = all with p < 0.05 in markers_model_results.csv
FOCUS_ANALYSES <- NULL

d   <- read.csv("output/markers_per_nerve.csv")
res <- read.csv("output/markers_model_results.csv")
pat <- read.csv("output/markers_patient_ratios.csv")
d   <- d[d$group %in% c("NI-TB", "InvN"), ]
d$logv <- log(d$normalised)

# deviation from the patient's group median, scaled by the group's spread
key <- paste(d$analysis, d$patient_label, d$group)
d$dev <- abs(d$logv - ave(d$logv, key, FUN = median)) / pmax(ave(d$logv, key, FUN = mad), 1e-6)
d$pct <- ave(d$logv, key, FUN = function(v) rank(v) / length(v))   # within-group percentile

# one row per nerve: mean deviation over analyses + the values for checking
score <- aggregate(dev ~ image_name + nerve_id + patient_label + group + nerve_area, data = d, FUN = mean)
names(score)[names(score) == "dev"] <- "score"
wide <- reshape(d[, c("image_name", "nerve_id", "analysis", "normalised", "pct")],
                idvar = c("image_name", "nerve_id"), timevar = "analysis", direction = "wide")
score <- merge(score, wide, by = c("image_name", "nerve_id"))

# skip small nerves
min_area <- ave(score$nerve_area, score$image_name, FUN = function(a) quantile(a, MIN_AREA_QUANT))
score <- score[score$nerve_area >= min_area, ]

score <- score[order(score$patient_label, score$group, score$score), ]
score$rank <- ave(score$score, score$patient_label, score$group, FUN = seq_along)
top <- score[score$rank <= N_PER_GROUP, ]
num <- sapply(top, is.numeric)
top[num] <- lapply(top[num], round, 3)
write.csv(top, "output/representative_nerves.csv", row.names = FALSE)

# which patient best reflects the overall effect (log distance to the model fold, averaged over
# the focus analyses, i.e. the effects the example images should show)
if (is.null(FOCUS_ANALYSES)) FOCUS_ANALYSES <- res$analysis[res$p_value < 0.05]
cat("Patient choice based on:", paste(FOCUS_ANALYSES, collapse = ", "), "\n\n")
pat <- merge(pat[pat$analysis %in% FOCUS_ANALYSES, ], res[, c("analysis", "fold")], by = "analysis")
pat$dist <- abs(log(pat$ratio_InvN_vs_NI_TB) - log(pat$fold))
rp <- aggregate(dist ~ patient, data = pat, FUN = mean)
rp <- rp[order(rp$dist), ]
names(rp)[2] <- "mean_log_distance_to_overall_effect"
n_nerves <- aggregate(nerve_id ~ patient_label + group, data = unique(d[, c("patient_label", "group", "nerve_id", "image_name")]), FUN = length)
rp$n_NI_TB <- n_nerves$nerve_id[match(paste(rp$patient, "NI-TB"), paste(n_nerves$patient_label, n_nerves$group))]
rp$n_InvN  <- n_nerves$nerve_id[match(paste(rp$patient, "InvN"),  paste(n_nerves$patient_label, n_nerves$group))]
rp$enough_nerves <- rp$n_NI_TB >= MIN_NERVES & rp$n_InvN >= MIN_NERVES
rp <- rp[order(!rp$enough_nerves, rp$mean_log_distance_to_overall_effect), ]
write.csv(rp, "output/representative_patients.csv", row.names = FALSE)

cat("Patients ranked by how typical their InvN/NI-TB effect is (lower = more typical):\n")
print(rp, row.names = FALSE)
best <- rp$patient[1]
cat(sprintf("\nTop candidate nerves for %s (lower score = more typical; pct = percentile within patient & group):\n", best))
show <- top[top$patient_label == best, c("group", "image_name", "nerve_id", "nerve_area", "score",
                                         grep("^pct\\.", names(top), value = TRUE))]
print(show, row.names = FALSE)
cat("\nFull list: output/representative_nerves.csv\n")
