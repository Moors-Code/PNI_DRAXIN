#!/bin/bash
#SBATCH --job-name=mageck_all_variants
#SBATCH --output=mageck_all_variants_%j.log
#SBATCH --error=mageck_all_variants_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --partition=normal

# Runs 7 MAGECK test variants: 3 input matrices (raw / T0_only-filtered /
# T0_plus_control-filtered) x normalization method (control / median / total),
# all comparing Draxin1,Draxin2 vs Ctr1,Ctr2, always --paired.
#
# Fill in before running:
#   <RAW_CSV>             -> path to Library_and_Screen.csv (unfiltered)
#   <T0_ONLY_CSV>         -> path to MAGECK_input_voom_T0_only_RAW_COUNTS.csv
#   <T0_PLUS_CONTROL_CSV> -> path to MAGECK_input_voom_T0_plus_control_RAW_COUNTS.csv
#   <NEG_CONTROLS_TXT>    -> path to neg_controls_clean.txt
#   <OUTPUT_DIR>          -> where all mageck output files should be written
# (conda env is mageck_env - already set below)

source ~/.bashrc
conda activate mageck_env

RAW_CSV="<RAW_CSV>"
T0_ONLY_CSV="<T0_ONLY_CSV>"
T0_PLUS_CONTROL_CSV="<T0_PLUS_CONTROL_CSV>"
NEG_CONTROLS="<NEG_CONTROLS_TXT>"
OUTPUT_DIR="<OUTPUT_DIR>"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# -----------------------------------------------------------
# 1) T0_ONLY_PAIRED_CONTROL
# -----------------------------------------------------------
mageck test \
  -k "$T0_ONLY_CSV" \
  -t Draxin1,Draxin2 \
  -c Ctr1,Ctr2 \
  --paired \
  --norm-method control \
  --control-sgrna "$NEG_CONTROLS" \
  --output-prefix T0_ONLY_PAIRED_CONTROL

# -----------------------------------------------------------
# 2) T0_ONLY_PAIRED_MEDIAN
# -----------------------------------------------------------
mageck test \
  -k "$T0_ONLY_CSV" \
  -t Draxin1,Draxin2 \
  -c Ctr1,Ctr2 \
  --paired \
  --norm-method median \
  --output-prefix T0_ONLY_PAIRED_MEDIAN

# -----------------------------------------------------------
# 3) T0_PLUS_CONTROL_PAIRED_CONTROL
# -----------------------------------------------------------
mageck test \
  -k "$T0_PLUS_CONTROL_CSV" \
  -t Draxin1,Draxin2 \
  -c Ctr1,Ctr2 \
  --paired \
  --norm-method control \
  --control-sgrna "$NEG_CONTROLS" \
  --output-prefix T0_PLUS_CONTROL_PAIRED_CONTROL

# -----------------------------------------------------------
# 4) T0_PLUS_CONTROL_PAIRED_MEDIAN
# -----------------------------------------------------------
mageck test \
  -k "$T0_PLUS_CONTROL_CSV" \
  -t Draxin1,Draxin2 \
  -c Ctr1,Ctr2 \
  --paired \
  --norm-method median \
  --output-prefix T0_PLUS_CONTROL_PAIRED_MEDIAN

# -----------------------------------------------------------
# 5) T0_PLUS_CONTROL_PAIRED_TOTAL
# -----------------------------------------------------------
mageck test \
  -k "$T0_PLUS_CONTROL_CSV" \
  -t Draxin1,Draxin2 \
  -c Ctr1,Ctr2 \
  --paired \
  --norm-method total \
  --output-prefix T0_PLUS_CONTROL_PAIRED_TOTAL

# -----------------------------------------------------------
# 6) NO_PRE_FILTERING_PAIRED_CONTROL
# -----------------------------------------------------------
mageck test \
  -k "$RAW_CSV" \
  -t Draxin1,Draxin2 \
  -c Ctr1,Ctr2 \
  --paired \
  --norm-method control \
  --control-sgrna "$NEG_CONTROLS" \
  --output-prefix NO_PRE_FILTERING_PAIRED_CONTROL

# -----------------------------------------------------------
# 7) NO_PRE_FILTERING_PAIRED_MEDIAN
# -----------------------------------------------------------
mageck test \
  -k "$RAW_CSV" \
  -t Draxin1,Draxin2 \
  -c Ctr1,Ctr2 \
  --paired \
  --norm-method median \
  --output-prefix NO_PRE_FILTERING_PAIRED_MEDIAN

echo "All 7 MAGECK runs complete. Output written to: $OUTPUT_DIR"
