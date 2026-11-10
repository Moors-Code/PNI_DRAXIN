library(tidyverse)

# =========================================================
# 01_build_Library_and_Screen.R
# Combines the 6 individual screen count files (S1-S6_counts.txt)
# with the library QC counts (merged_counts_final.csv) into a single
# Library_and_Screen.csv, matching the structure expected by 02_ScreenQC.Rmd.
#
# NOTE: this uses R's dplyr::full_join throughout, which matches rows by
# key value directly rather than assuming pre-sorted input - this avoids
# the sort/join locale bug found in the bash merge pipeline (see protocol
# notes on 05_bwa_alignment.sh for details).
# =========================================================

# ---------------------------------------------------------
# 1. Fill in before running: paths to your 6 screen count files
# ---------------------------------------------------------
screen_files <- c(
  S1 = "input_files/S1_counts.txt",
  S2 = "input_files/S2_counts.txt",
  S3 = "input_files/S3_counts.txt",
  S4 = "input_files/S4_counts.txt",
  S5 = "input_files/S5_counts.txt",
  S6 = "input_files/S6_counts.txt"
)

library_qc_file <- "input_files/merged_counts_final_libraryQC.csv"
output_file <- "input_files/Library_and_Screen.csv"

# Sample name mapping: S1-S6 -> the column names 02_ScreenQC.Rmd expects
# (confirmed against the original matrixScreen20250712.csv column order)
sample_rename <- c(
  S1 = "Ctr2",
  S2 = "Draxin2",
  S3 = "Before1",
  S4 = "Ctr1",
  S5 = "Draxin1",
  S6 = "Before2"
)

# ---------------------------------------------------------
# 2. Read each S#_counts.txt file
# Format: "  <count> <guide_id>" (whitespace-separated, no header)
# ---------------------------------------------------------
read_counts_txt <- function(path, sample_name) {
  read.table(path, header = FALSE, col.names = c("Count", "sgRNA"),
             stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    dplyr::select(sgRNA, !!sample_name := Count)
}

screen_list <- Map(read_counts_txt, screen_files, names(screen_files))

# ---------------------------------------------------------
# 3. Merge all 6 screen samples together (full_join keeps every guide seen
# in ANY sample; missing values filled with 0 afterward)
# ---------------------------------------------------------
screen_merged <- Reduce(function(x, y) full_join(x, y, by = "sgRNA"), screen_list)

# Rename S1-S6 columns to Ctr2/Draxin2/Before1/Ctr1/Draxin1/Before2
colnames(screen_merged)[match(names(sample_rename), colnames(screen_merged))] <- unname(sample_rename)

# ---------------------------------------------------------
# 4. Read the library QC merged counts
# Expected columns: "" (guide ID), rep1_counts, rep2_counts, rep3_counts, rep4_counts
# ---------------------------------------------------------
library_qc <- read.csv(library_qc_file, check.names = FALSE)
colnames(library_qc)[1] <- "sgRNA"

# ---------------------------------------------------------
# 5. Full join screen + library on sgRNA, fill missing with 0
# ---------------------------------------------------------
combined <- full_join(library_qc, screen_merged, by = "sgRNA")

count_cols <- setdiff(colnames(combined), "sgRNA")
combined[count_cols] <- lapply(combined[count_cols], function(x) {
  x[is.na(x)] <- 0
  x
})

# ---------------------------------------------------------
# 6. Add Gene column (strip trailing _<number> guide-index suffix)
# ---------------------------------------------------------
combined <- combined %>%
  mutate(Gene = sub("_[0-9]+$", "", sgRNA)) %>%
  dplyr::select(sgRNA, Gene, everything())

# ---------------------------------------------------------
# 7. Sanity checks before writing
# ---------------------------------------------------------
cat("Total guides in combined matrix:", nrow(combined), "\n")
cat("Expected (deduplicated reference size): 49351\n")
cat("Columns:", paste(colnames(combined), collapse = ", "), "\n")
cat("Any NA remaining?", any(is.na(combined)), "\n")

# ---------------------------------------------------------
# 8. Write output
# ---------------------------------------------------------
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write.csv(combined, output_file, row.names = FALSE, quote = FALSE)

cat("Written to:", output_file, "\n")
