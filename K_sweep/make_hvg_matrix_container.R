# =============================================================================
# Highly variable gene matrix for the K-sensitivity NMF sweep
#
# Submitted by run_hvg_matrix_container.sh, inside the Singularity container
# that holds the package environment of the original NMF analysis.
#
# Per-sample gene variance (modelGeneVar, lowess disabled) is combined across
# samples (combineVar) and the top 5000 highly variable genes are retained;
# their log-count matrix is saved as the input to the sweep.
#
#   in : data/SPE/PNI_SPE_QC_kcsweep.rds   quality-controlled SpatialExperiment
#   out: data/coldata/PNI_HVG_matrix_container.rds
# =============================================================================

library(SpatialExperiment)
library(scater)
library(scran)
library(batchelor)

SPE_FILE <- file.path("data", "SPE", "PNI_SPE_QC_kcsweep.rds")
OUT_FILE <- file.path("data", "coldata", "PNI_HVG_matrix_container.rds")
N_HVG    <- 5000

if (!file.exists(SPE_FILE)) stop("SpatialExperiment object not found: ", SPE_FILE)
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

spe <- readRDS(SPE_FILE)
colnames(spe) <- spe$barcode

set.seed(42)
smps <- unique(spe$sample_id)
vars <- lapply(smps, function(SMP) {
  out <- modelGeneVar(spe[, spe$sample_id == SMP], lowess = FALSE)
  return(out)
})

combVar <- combineVar(vars)
hvgs <- getTopHVGs(combVar, n = N_HVG)
m <- logcounts(spe)[hvgs, ]

saveRDS(m, OUT_FILE)
cat("HVG matrix saved. Dimensions:", dim(m), "\n")
