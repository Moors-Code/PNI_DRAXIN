# ============================================================
# NMF model fitting (run on Euler HPC cluster)
# ============================================================
#
# Selects highly variable genes across samples (combining per-sample
# variance models), then fits sparse NMF (K=30, L1=0.5) on the resulting
# expression matrix. This is the actual model-fitting step referenced
# by the local visualization/interpretation script
# (nmf_visualization_and_export_cleaned.Rmd) -- run this first, on
# Euler, since it's computationally heavy; the visualization script
# only ever reads its output, never re-fits the model.
#
# WHAT THIS PRODUCES
# --------------------
# Processed (processed/):
#   - PNI_NMF_L105-k30.rds
#       THE FITTED NMF MODEL. An RcppML nmf S4 object (@w, @h, @d).
#       Read by nmf_visualization_and_export_cleaned.Rmd downstream.
# ============================================================

library(scater)
library(scran)
library(RcppML)
library(batchelor)

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("processed"))
  stop("Run this script from the project folder (the one containing processed/).")

PROJECT_DIR   <- "."
INPUT_DIR     <- file.path(PROJECT_DIR, "input")
PROCESSED_DIR <- file.path(PROJECT_DIR, "processed")
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Select highly variable genes and fit the sparse NMF model ------
# NOTE: reads PNI_SPE_QC.rds from processed/ -- this is the output of
# the QC/filtering script (01_QC_filtering), not a raw input file.
spe <- readRDS(file.path(PROCESSED_DIR, "PNI_SPE_QC.rds"))
colnames(spe) <- spe$barcode

# ---- Across Sample HVGs ----
set.seed(42) # RND seed
# Per Sample Var Dec
smps <- unique(spe$sample_id)
vars <- lapply(smps, function(SMP) {
                 out <- modelGeneVar(spe[,spe$sample_id==SMP], lowess=FALSE)
                 return(out)
})

combVar <- combineVar(vars)
hvgs <- getTopHVGs(combVar, n=5000)
m <- logcounts(spe)[hvgs,]

#--- Run NMF ---
nmfout <- nmf(m,k=30,tol=1e-6,seed=1:3,L1=0.5)
saveRDS(nmfout, file.path(PROCESSED_DIR, "PNI_NMF_L105-k30.rds"))
cat("Saved fitted NMF model to:", file.path(PROCESSED_DIR, "PNI_NMF_L105-k30.rds"), "\n")
