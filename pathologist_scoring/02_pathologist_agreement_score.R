# ======================================================================
# Agreement between computational (NMF-based) and pathologist nerve /
# invasion annotations.
#
# Join key: NerveID, confirmed identical in both files via the
# Legend_Instructions sheet ("each row corresponds to an image file
# named exactly as the NerveID") -- both files independently have
# exactly 107 nerves.
#
# Reports THREE separate things, deliberately not collapsed into one
# number:
#   1) Pathologist-rejected computational nerve detections ("No_Nerve").
#      One-directional: only computational candidates were reviewed, so
#      missed nerves (detection false negatives) are not measured.
#   2) Invasion-status agreement (Invaded/NonInvaded), restricted to
#      nerves BOTH sides gave a valid Invaded/NonInvaded call for
#   3) Subtype-level detail (PNI/ENA/ENI), for context on where
#      disagreements concentrate
#
# LABEL SAFETY (revision fix)
# ---------------------------
# Earlier versions computed accuracy as sum(diag(table(comp, path))).
# That is only correct when both label vectors contain exactly the same
# categories in the same order. The export script (07, line 142) can
# emit "AMBIGUOUS:<values>" when a nerve's spots disagree, and the
# pathologist sheet may contain typos or blanks. Either one shifts the
# table so the diagonal pairs up the wrong categories; in the extreme
# case (no overlap at all) the table collapses and reports perfect
# agreement.
#
# This version:
#   - maps an explicit whitelist of spellings (case-insensitive, at most
#     one space/_/- between words) to Invaded, NonInvaded, No_Nerve;
#     anything else (e.g. "Invaded?", "In-vaded") stays unexpected and
#     is excluded
#   - rejects missing/blank/duplicate NerveIDs and warns if the join
#     does not contain exactly 107 matched nerves
#   - builds the confusion matrix from factors with FIXED levels
#     c("Invaded","NonInvaded"), so it is always 2x2 and aligned
#   - computes accuracy by direct element-wise comparison and
#     cross-checks it against the diagonal
#   - never silently drops rows: AMBIGUOUS / unexpected / missing labels
#     are counted, listed in excluded_from_agreement.csv, and a "strict"
#     accuracy (excluded nerves counted as disagreements) is reported
#     alongside the main figure
#
# SAFETY: read-only on both input files.
#
# WHAT THIS PRODUCES
# --------------------
# Tables (output/08_pathologist_agreement/tables/):
#   - nerve_annotation_agreement_merged.csv
#       Full joined table, one row per nerve both sides recorded, with
#       an agreement_class column.
#   - confusion_matrix.csv
#       Fixed 2x2 computational vs pathologist Invaded/NonInvaded.
#   - observed_label_crosstab.csv
#       Every RAW label combination actually observed (incl. AMBIGUOUS,
#       typos, NA), long format -- an audit, not used in any metric.
#   - nerve_accounting.csv
#       Every NerveID from either file counted exactly once: agree,
#       disagree, excluded (by side), No_Nerve, unmatched (by file).
#   - agreement_stats_summary.csv
#       THE MAIN OUTPUT TABLE. One row: all headline numbers.
#   - disagreements.csv
#       Evaluable nerves where the calls differ.
#   - excluded_from_agreement.csv
#       Nerves excluded from metric 2 because of an ambiguous,
#       unexpected or missing label (saved even if empty).
#   - subtype_breakdown.csv
#   - pathologist_internal_inconsistencies.csv (saved even if empty)
#   - agreement_report.txt
#       Full human-readable report.
# ======================================================================

library(readxl)
library(dplyr)
library(irr)      # for Cohen's kappa -- install.packages("irr") if missing

# ---- Namespace pinning ---------------------------------------------------
# Defensive, project-wide standard practice (see other scripts in this
# project): pins dplyr verbs explicitly in case another loaded package
# (e.g. plyr, or Bioconductor packages in a shared session) masks them.
rename    <- dplyr::rename
count     <- dplyr::count
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
first     <- dplyr::first
relocate  <- dplyr::relocate

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input"))
  stop("Run this script from the project folder (the one containing input/).")

PROJECT_DIR   <- "."
INPUT_DIR     <- file.path(PROJECT_DIR, "input")

# Where the computational nerve annotation export lives (produced by
# the 07 pathologist_nerve_annotation_export script). Read-only here.
PATHOLOGIST_EXPORT_TABLES_DIR <- file.path(PROJECT_DIR, "output", "01_pathologist_annotation_export", "tables")

ANALYSIS_NAME <- "02_pathologist_agreement"
TABLES_DIR    <- file.path(PROJECT_DIR, "output", ANALYSIS_NAME, "tables")
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Label vocabulary ----------------------------------------------------
INVASION_LEVELS <- c("Invaded", "NonInvaded")   # fixed order for every table
NO_NERVE_LABEL  <- "No_Nerve"
SUBTYPES        <- c("PNI", "ENA", "ENI")

EXPECTED_N_NERVES <- 107   # both files are expected to contain exactly these nerves

# Map an explicit whitelist of spelling variants onto the canonical
# vocabulary. Matching is case-insensitive, and at most ONE separator
# (space, underscore or hyphen) is allowed, and only between the words:
#   Invaded    <- "invaded"
#   NonInvaded <- "NonInvaded", "Non-Invaded", "non_invaded", "non invaded"
#   No_Nerve   <- "No_Nerve", "No Nerve", "no-nerve", "nonerve"
# Anything else ("Invaded?", "In-vaded", "Non--Invaded", "AMBIGUOUS:...")
# is kept verbatim so it is reported as unexpected, never silently re-coded.
LABEL_PATTERNS <- c(
  Invaded    = "^invaded$",
  NonInvaded = "^non[ _-]?invaded$",
  No_Nerve   = "^no[ _-]?nerve$"
)
canonical_label <- function(x) {
  x   <- trimws(as.character(x))
  out <- x
  for (lab in names(LABEL_PATTERNS)) {
    out[which(grepl(LABEL_PATTERNS[[lab]], x, ignore.case = TRUE))] <- lab
  }
  out[which(is.na(x) | x == "")] <- NA_character_
  out
}
stopifnot(identical(names(LABEL_PATTERNS), c(INVASION_LEVELS, NO_NERVE_LABEL)))

safe_div <- function(a, b) if (is.na(b) || b == 0) NA_real_ else a / b

# ---- 1. Load both files --------------------------------------------------
# NOTE: place nerve_annotation_for_pathologist_Daniela.xlsx (the
# pathologist's completed sheet) into input/ before
# running.
path_annot <- read_excel(file.path(INPUT_DIR, "nerve_annotation_for_pathologist_Daniela.xlsx"),
                         sheet = "Nerve_Annotations") %>%
  mutate(Invasion_Status_raw = Invasion_Status,
         Invasion_Status     = canonical_label(Invasion_Status),
         Subtype             = toupper(trimws(as.character(Subtype))))

comp_annot <- read.csv(file.path(PATHOLOGIST_EXPORT_TABLES_DIR, "nerve_annotations.csv"),
                       stringsAsFactors = FALSE) %>%
  mutate(Invaded_label_raw = Invaded_label,
         Invaded_label     = canonical_label(Invaded_label))

# ---- 1a. Reject missing or blank join keys --------------------------------
comp_annot$NerveID <- trimws(as.character(comp_annot$NerveID))
path_annot$NerveID <- trimws(as.character(path_annot$NerveID))
n_bad_comp_ids <- sum(is.na(comp_annot$NerveID) | comp_annot$NerveID == "")
n_bad_path_ids <- sum(is.na(path_annot$NerveID) | path_annot$NerveID == "")
if (n_bad_comp_ids > 0 || n_bad_path_ids > 0) {
  stop("Missing or blank NerveID values found (computational: ", n_bad_comp_ids,
       ", pathologist: ", n_bad_path_ids, "). Fix the join keys before computing agreement.")
}

# ---- 1b. Duplicate keys would silently inflate the join -----------------
dup_comp <- unique(comp_annot$NerveID[duplicated(comp_annot$NerveID)])
dup_path <- unique(path_annot$NerveID[duplicated(path_annot$NerveID)])
if (length(dup_comp) > 0 || length(dup_path) > 0) {
  stop("Duplicate NerveIDs found -- fix before computing agreement.\n",
       "  Computational: ", paste(dup_comp, collapse = ", "), "\n",
       "  Pathologist:   ", paste(dup_path, collapse = ", "))
}

# ---- 1c. Vocabulary check on each file ----------------------------------
unexpected_path <- path_annot %>%
  filter(!(Invasion_Status %in% c(INVASION_LEVELS, NO_NERVE_LABEL)))
unexpected_comp <- comp_annot %>%
  filter(!(Invaded_label %in% INVASION_LEVELS))   # includes AMBIGUOUS:... and NA

# ---- 2. Internal consistency check on the pathologist's own sheet ------
# (recomputes what the red-highlighting in the xlsx was flagging, since
# cell colour doesn't survive plain extraction). %in% is used so NA
# values are compared safely rather than dropped by filter().
inconsistent <- path_annot %>%
  filter(
    (Invasion_Status %in% "Invaded" & !(Subtype %in% SUBTYPES)) |
      (Invasion_Status %in% c("NonInvaded", NO_NERVE_LABEL) & !is.na(Subtype) & Subtype != "NA")
  )

# ---- 3. Join on NerveID, check for unmatched rows on either side -------
merged <- comp_annot %>% inner_join(path_annot, by = "NerveID")
unmatched_comp <- anti_join(comp_annot, path_annot, by = "NerveID")
unmatched_path <- anti_join(path_annot, comp_annot, by = "NerveID")

# ---- 3b. Expected-size check ----------------------------------------------
# Equal row counts do not guarantee equal ID sets, so the join itself is
# checked. Warnings (not stops) so the outputs can still be inspected.
size_warnings <- character(0)
if (nrow(comp_annot) != EXPECTED_N_NERVES || nrow(path_annot) != EXPECTED_N_NERVES) {
  size_warnings <- c(size_warnings, sprintf(
    "Input files do not both contain %d rows (computational: %d, pathologist: %d).",
    EXPECTED_N_NERVES, nrow(comp_annot), nrow(path_annot)))
}
if (nrow(merged) != EXPECTED_N_NERVES || nrow(unmatched_comp) > 0 || nrow(unmatched_path) > 0) {
  size_warnings <- c(size_warnings, sprintf(
    "Join contains %d matched nerves, not %d. Review unmatched NerveIDs before interpreting agreement.",
    nrow(merged), EXPECTED_N_NERVES))
}
for (w in size_warnings) warning(w, call. = FALSE)

# Every merged nerve gets exactly one class -- nothing is dropped silently
merged <- merged %>%
  mutate(agreement_class = case_when(
    Invasion_Status %in% NO_NERVE_LABEL     ~ "No_Nerve (pathologist)",
    !(Invasion_Status %in% INVASION_LEVELS) ~ "Excluded: unexpected/missing pathologist status",
    !(Invaded_label %in% INVASION_LEVELS)   ~ "Excluded: ambiguous/unexpected computational label",
    Invaded_label == Invasion_Status        ~ "Agree",
    TRUE                                    ~ "Disagree"
  ))

# ======================================================================
# METRIC 1: Pathologist-rejected computational nerve detections
# ======================================================================
n_no_nerve <- sum(merged$Invasion_Status %in% NO_NERVE_LABEL)
no_nerve_rows <- merged %>%
  filter(Invasion_Status %in% NO_NERVE_LABEL) %>%
  select(any_of(c("NerveID", "sample_id", "n_spots", "Invaded_label_raw", "Pathologist_Notes")))

# ======================================================================
# METRIC 2: Invasion-status agreement
# ======================================================================
# Confirmed = pathologist gave a valid Invaded/NonInvaded call
confirmed <- merged %>% filter(Invasion_Status %in% INVASION_LEVELS)
# Evaluable = confirmed AND computational call is also valid
evaluable <- confirmed %>% filter(Invaded_label %in% INVASION_LEVELS)

excluded <- merged %>%
  filter(startsWith(agreement_class, "Excluded")) %>%
  select(any_of(c("NerveID", "sample_id", "n_spots", "Invaded_label_raw",
                  "Invasion_Status_raw", "agreement_class", "Subtype", "Pathologist_Notes")))

n_confirmed          <- nrow(confirmed)
n_evaluable          <- nrow(evaluable)
n_excluded_path <- sum(merged$agreement_class == "Excluded: unexpected/missing pathologist status")

# Confirmed nerves whose computational label is not Invaded/NonInvaded,
# broken down by reason
excl_comp_labels        <- confirmed$Invaded_label[!(confirmed$Invaded_label %in% INVASION_LEVELS)]
n_excluded_comp         <- length(excl_comp_labels)   # all invalid computational labels
n_excluded_ambiguous    <- sum(!is.na(excl_comp_labels) & startsWith(excl_comp_labels, "AMBIGUOUS:"))
n_excluded_missing_comp <- sum(is.na(excl_comp_labels))
n_excluded_other_comp   <- n_excluded_comp - n_excluded_ambiguous - n_excluded_missing_comp
stopifnot(n_excluded_comp == n_confirmed - n_evaluable)

# Fixed-level factors => table is ALWAYS 2x2 with rows and columns aligned
comp_f <- factor(evaluable$Invaded_label,   levels = INVASION_LEVELS)
path_f <- factor(evaluable$Invasion_Status, levels = INVASION_LEVELS)
conf_matrix <- table(Computational = comp_f, Pathologist = path_f)
stopifnot(identical(rownames(conf_matrix), INVASION_LEVELS),
          identical(colnames(conf_matrix), INVASION_LEVELS))

# Accuracy by direct comparison, cross-checked against the (now valid) diagonal
n_agree <- sum(evaluable$Invaded_label == evaluable$Invasion_Status)
stopifnot(n_agree == sum(diag(conf_matrix)))
accuracy        <- safe_div(n_agree, n_evaluable)
# Strict: among pathologist-confirmed nerves, invalid computational calls
# count as disagreements. Nerves with an invalid/missing PATHOLOGIST status
# are not in either denominator (no reference call to compare against).
accuracy_strict <- safe_div(n_agree, n_confirmed)

# ---- Full accounting: every nerve (from either file) counted exactly once --
ACCOUNTING_ORDER <- c(
  "Agree", "Disagree",
  "Excluded: ambiguous/unexpected computational label",
  "Excluded: unexpected/missing pathologist status",
  "No_Nerve (pathologist)",
  "Unmatched: computational file only",
  "Unmatched: pathologist file only"
)
nerve_accounting <- data.frame(
  category = ACCOUNTING_ORDER,
  n = c(vapply(ACCOUNTING_ORDER[1:5], function(k) sum(merged$agreement_class == k), integer(1)),
        nrow(unmatched_comp), nrow(unmatched_path)),
  row.names = NULL
)
n_total_ids <- length(union(comp_annot$NerveID, path_annot$NerveID))
stopifnot(sum(nerve_accounting$n) == n_total_ids,
          all(merged$agreement_class %in% ACCOUNTING_ORDER))
nerve_accounting <- rbind(nerve_accounting,
                          data.frame(category = "TOTAL distinct NerveIDs", n = n_total_ids))

# Every RAW label combination actually observed (audit only; all
# agreement calculations use the canonical columns)
observed_crosstab <- as.data.frame(
  table(Computational = merged$Invaded_label_raw,
        Pathologist   = merged$Invasion_Status_raw, useNA = "ifany"),
  responseName = "n") %>%
  filter(n > 0)

# Cohen's kappa on the aligned factors; undefined when there is no
# variation (e.g. everything in one category) -> reported as NA
kappa_value <- NA_real_
kappa_p     <- NA_real_
if (n_evaluable >= 2) {
  kr <- tryCatch(kappa2(data.frame(comp_f, path_f)), error = function(e) NULL)
  if (!is.null(kr) && is.finite(kr$value)) {
    kappa_value <- kr$value
    kappa_p     <- kr$p.value
  }
}

# Sensitivity/specificity/PPV/NPV, pathologist as ground truth,
# read from named cells (not positions)
TP <- as.integer(conf_matrix["Invaded",    "Invaded"])
TN <- as.integer(conf_matrix["NonInvaded", "NonInvaded"])
FP <- as.integer(conf_matrix["Invaded",    "NonInvaded"])
FN <- as.integer(conf_matrix["NonInvaded", "Invaded"])

sensitivity <- safe_div(TP, TP + FN)
specificity <- safe_div(TN, TN + FP)
ppv         <- safe_div(TP, TP + FP)
npv         <- safe_div(TN, TN + FN)

disagreements <- evaluable %>%
  filter(Invaded_label != Invasion_Status) %>%
  select(any_of(c("NerveID", "sample_id", "n_spots", "Invaded_label",
                  "Invasion_Status", "Subtype", "Pathologist_Notes")))

# ======================================================================
# METRIC 3: Subtype breakdown among pathologist-Invaded nerves
# (includes ambiguous computational labels, as context)
# ======================================================================
subtype_breakdown <- confirmed %>%
  filter(Invasion_Status == "Invaded") %>%
  dplyr::count(Subtype, Invaded_label_raw)

# ======================================================================
# REPORT (printed to console and written to agreement_report.txt)
# ======================================================================
r3 <- function(x) round(x, 3)

print_report <- function() {
  cat("======================================================================\n")
  cat("Agreement between computational (NMF-based) and pathologist annotations\n")
  cat("======================================================================\n\n")
  
  cat("Files:\n")
  cat("  Pathologist file:  ", nrow(path_annot), "rows\n")
  cat("  Computational file:", nrow(comp_annot), "rows\n")
  cat("  Successfully joined:", nrow(merged), "of", nrow(comp_annot),
      "(expected", EXPECTED_N_NERVES, ")\n")
  for (w in size_warnings) cat("WARNING:", w, "\n")
  if (nrow(unmatched_comp) > 0) {
    cat("WARNING:", nrow(unmatched_comp), "computational nerves have no pathologist annotation:\n")
    print(unmatched_comp$NerveID)
  }
  if (nrow(unmatched_path) > 0) {
    cat("WARNING:", nrow(unmatched_path), "pathologist rows have no matching computational nerve:\n")
    print(unmatched_path$NerveID)
  }
  
  cat("\nLabel vocabulary check:\n")
  cat("  Computational labels outside {Invaded, NonInvaded}:", nrow(unexpected_comp), "\n")
  if (nrow(unexpected_comp) > 0) print(table(unexpected_comp$Invaded_label_raw, useNA = "ifany"))
  cat("  Pathologist statuses outside {Invaded, NonInvaded, No_Nerve}:", nrow(unexpected_path), "\n")
  if (nrow(unexpected_path) > 0) print(table(unexpected_path$Invasion_Status_raw, useNA = "ifany"))
  
  cat("\nInternally inconsistent pathologist rows (Invasion_Status vs Subtype):", nrow(inconsistent), "\n")
  if (nrow(inconsistent) > 0) {
    cat("Review these before trusting the agreement numbers below:\n")
    print(inconsistent %>% select(any_of(c("NerveID", "Invasion_Status", "Subtype", "Pathologist_Notes"))))
  }
  
  cat("\n--- Nerve accounting (each NerveID counted exactly once) ---\n")
  print(nerve_accounting, row.names = FALSE)
  
  cat("\n--- METRIC 1: Pathologist-rejected computational nerve detections ---\n")
  cat(n_no_nerve, "of", nrow(merged), "computational nerve candidates (",
      round(100 * safe_div(n_no_nerve, nrow(merged)), 1),
      "%) were judged by the pathologist NOT to be real nerve structures.\n")
  cat("Note: one-directional. Only computational candidates were reviewed, so\n")
  cat("nerves the pipeline missed (detection false negatives) are not measured.\n")
  if (n_no_nerve > 0) print(no_nerve_rows)
  
  cat("\n--- METRIC 2: Invasion-status agreement ---\n")
  cat("Pathologist-confirmed nerves with a valid status:", n_confirmed, "\n")
  cat("  of which excluded (invalid computational label):", n_excluded_comp, "\n")
  cat("    AMBIGUOUS:", n_excluded_ambiguous,
      "| missing:", n_excluded_missing_comp,
      "| other unexpected:", n_excluded_other_comp, "\n")
  cat("  evaluable (both sides Invaded/NonInvaded):", n_evaluable, "\n")
  cat("Excluded for unexpected/missing pathologist status:", n_excluded_path, "\n")
  if (nrow(excluded) > 0) {
    cat("\nExcluded nerves:\n")
    print(excluded)
  }
  
  cat("\nConfusion matrix (rows = computational, cols = pathologist; fixed 2x2):\n")
  print(conf_matrix)
  cat("\nRaw agreement (evaluable nerves):", r3(accuracy), "(", n_agree, "/", n_evaluable, ")\n")
  cat("Strict agreement (invalid computational calls on pathologist-confirmed nerves\n")
  cat("  counted as disagreement; invalid pathologist statuses are not in the denominator):",
      r3(accuracy_strict), "(", n_agree, "/", n_confirmed, ")\n")
  cat("Cohen's kappa:", r3(kappa_value), "(p =", format.pval(kappa_p, digits = 3), ")\n")
  cat("(Kappa and the ratios below are computed on the", n_evaluable, "evaluable nerves only;\n")
  cat(" report them together with the exclusion counts and strict agreement.)\n")
  cat("Treating pathologist annotation as ground truth:\n")
  cat("  Sensitivity:", r3(sensitivity), "\n")
  cat("  Specificity:", r3(specificity), "\n")
  cat("  PPV:        ", r3(ppv), "\n")
  cat("  NPV:        ", r3(npv), "\n")
  
  cat("\nDisagreements (", nrow(disagreements), "of", n_evaluable, "evaluable):\n")
  print(disagreements)
  
  cat("\n--- METRIC 3: Subtype breakdown among pathologist-Invaded nerves ---\n")
  print(subtype_breakdown)
  cat("\nNote: subtype categories (PNI/ENA/ENI) represent a biological gradient\n")
  cat("(epineural association -> perineural invasion -> endoneural invasion).\n")
  cat("If disagreements concentrate in ENA (the mildest/most ambiguous subtype),\n")
  cat("that's a different, more forgivable kind of miss than disagreement on\n")
  cat("clear-cut PNI/ENI cases -- worth distinguishing in how you report this.\n")
}

print_report()

# ======================================================================
# SAVE ALL RESULTS
# ======================================================================
write.csv(merged, file.path(TABLES_DIR, "nerve_annotation_agreement_merged.csv"), row.names = FALSE)

conf_matrix_df <- as.data.frame.matrix(conf_matrix)
conf_matrix_df$Computational <- rownames(conf_matrix_df)
conf_matrix_df <- conf_matrix_df %>% relocate(Computational)
write.csv(conf_matrix_df, file.path(TABLES_DIR, "confusion_matrix.csv"), row.names = FALSE)

write.csv(observed_crosstab, file.path(TABLES_DIR, "observed_label_crosstab.csv"), row.names = FALSE)
write.csv(nerve_accounting,  file.path(TABLES_DIR, "nerve_accounting.csv"),        row.names = FALSE)

stats_summary <- data.frame(
  n_computational_nerves             = nrow(comp_annot),
  n_pathologist_rows                 = nrow(path_annot),
  n_expected_nerves                  = EXPECTED_N_NERVES,
  n_joined                           = nrow(merged),
  n_unmatched_computational          = nrow(unmatched_comp),
  n_unmatched_pathologist            = nrow(unmatched_path),
  n_unexpected_computational_labels  = nrow(unexpected_comp),
  n_unexpected_pathologist_labels    = nrow(unexpected_path),
  n_pathologist_rejected_detections  = n_no_nerve,
  pct_pathologist_rejected_detections = round(100 * safe_div(n_no_nerve, nrow(merged)), 1),
  n_confirmed_nerves                 = n_confirmed,
  n_excluded_invalid_computational   = n_excluded_comp,
  n_excluded_ambiguous_computational = n_excluded_ambiguous,
  n_excluded_missing_computational   = n_excluded_missing_comp,
  n_excluded_other_computational     = n_excluded_other_comp,
  n_excluded_unexpected_pathologist  = n_excluded_path,
  n_evaluable                        = n_evaluable,
  n_agree                            = n_agree,
  n_disagreements                    = nrow(disagreements),
  raw_agreement_accuracy             = r3(accuracy),
  strict_agreement_accuracy          = r3(accuracy_strict),
  cohens_kappa                       = r3(kappa_value),
  kappa_p_value                      = kappa_p,
  sensitivity                        = r3(sensitivity),
  specificity                        = r3(specificity),
  ppv                                = r3(ppv),
  npv                                = r3(npv),
  TP = TP, TN = TN, FP = FP, FN = FN
)
write.csv(stats_summary, file.path(TABLES_DIR, "agreement_stats_summary.csv"), row.names = FALSE)

write.csv(disagreements,     file.path(TABLES_DIR, "disagreements.csv"),           row.names = FALSE)
write.csv(excluded,          file.path(TABLES_DIR, "excluded_from_agreement.csv"), row.names = FALSE)
write.csv(subtype_breakdown, file.path(TABLES_DIR, "subtype_breakdown.csv"),       row.names = FALSE)
write.csv(inconsistent,      file.path(TABLES_DIR, "pathologist_internal_inconsistencies.csv"), row.names = FALSE)

report_path <- file.path(TABLES_DIR, "agreement_report.txt")
sink(report_path)
tryCatch(print_report(), finally = sink())

cat("\n============================================================\n")
cat("All results saved to", TABLES_DIR, ":\n")
cat("  - nerve_annotation_agreement_merged.csv  (full joined table + agreement_class)\n")
cat("  - confusion_matrix.csv                   (fixed 2x2)\n")
cat("  - observed_label_crosstab.csv            (all observed raw label pairs)\n")
cat("  - nerve_accounting.csv                   (every NerveID counted once)\n")
cat("  - agreement_stats_summary.csv            (all headline numbers)\n")
cat("  - disagreements.csv\n")
cat("  - excluded_from_agreement.csv\n")
cat("  - subtype_breakdown.csv\n")
cat("  - pathologist_internal_inconsistencies.csv\n")
cat("  - agreement_report.txt                   (full human-readable report)\n")
cat("============================================================\n")
