#!/bin/bash
#SBATCH --job-name=nmf_ksweep
#SBATCH --array=0-5
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=ksweep_%A_%a.log

# =============================================================================
# K-sensitivity NMF sweep, step 2 of 3 (run on an HPC cluster with SLURM)
#
# SLURM array job, one task per K. Each task fits a sparse NMF (RcppML) to the
# HVG matrix from step 1 and saves the fit together with the relative Frobenius
# reconstruction error of that fit.
#
# The NMF settings reported in the Methods are set here: K as given below,
# tol = 1e-6, seed = 1:3, L1 = 0.5.
#
#   in : data/coldata/PNI_HVG_matrix_container.rds   (run_hvg_matrix_container.sh)
#   out: data/coldata/ksweep/PNI_NMF_L105-k<KK>.rds  one NMF fit per K
#        data/coldata/ksweep/recon_error_k<KK>.csv   one-row error file per K
#
# Runs inside the same container as step 1 (built from the Docker image
# bioconductor/bioconductor_docker:RELEASE_3_17, Bioconductor 3.17, R 4.3, with
# RcppML 0.5.5 installed), with R_LIBS_USER unset. Set CONTAINER to the image
# before submitting.
#
# Run order: run_hvg_matrix_container.sh, this script, then 04_Ksweep_analysis.R
# locally, which combines the per-K error files and compares the fits.
# =============================================================================

CONTAINER=<path to PNI_Image.sif>

if [ ! -f "$CONTAINER" ]; then
  echo "Set CONTAINER to the Singularity image (currently: $CONTAINER)" >&2
  exit 1
fi

unset R_LIBS_USER

KS=(15 20 25 30 35 40)
K=${KS[$SLURM_ARRAY_TASK_ID]}

mkdir -p data/coldata/ksweep

singularity exec "$CONTAINER" Rscript -e "
library(RcppML)
m <- readRDS('data/coldata/PNI_HVG_matrix_container.rds')
fit <- nmf(m, k=$K, tol=1e-6, seed=1:3, L1=0.5)
saveRDS(fit, sprintf('data/coldata/ksweep/PNI_NMF_L105-k%02d.rds', $K))
m_hat <- fit@w %*% diag(fit@d) %*% fit@h
rel_err <- sqrt(sum((m - m_hat)^2)) / sqrt(sum(m^2))
write.csv(data.frame(K=$K, rel_frob_error=rel_err),
          sprintf('data/coldata/ksweep/recon_error_k%02d.csv', $K), row.names=FALSE)
"
