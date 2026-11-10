#!/bin/bash
#SBATCH --job-name=hvg_matrix_container
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=2
#SBATCH --output=hvg_matrix_container_%j.log

# =============================================================================
# K-sensitivity NMF sweep, step 1 of 3 (run on an HPC cluster with SLURM)
#
# Computes the highly variable gene matrix the sweep is run on: per-sample
# HVGs from the Visium SpatialExperiment object (modelGeneVar), combined
# across samples (combineVar), top 5000 retained.
#
#   in : data/SPE/PNI_SPE_QC_kcsweep.rds
#   out: data/coldata/PNI_HVG_matrix_container.rds
#
# Runs inside a Singularity/Apptainer container built from the Docker image
# bioconductor/bioconductor_docker:RELEASE_3_17 (Bioconductor 3.17, R 4.3).
# R_LIBS_USER is unset so that only the container's own packages are used and a
# personal library cannot interfere. Set CONTAINER to the image before
# submitting.
#
# Run order: this script, then run_ksweep_container.sh, then
# 04_Ksweep_analysis.R locally.
# =============================================================================

CONTAINER=<path to PNI_Image.sif>

if [ ! -f "$CONTAINER" ]; then
  echo "Set CONTAINER to the Singularity image (currently: $CONTAINER)" >&2
  exit 1
fi

unset R_LIBS_USER
mkdir -p data/coldata

singularity exec "$CONTAINER" Rscript make_hvg_matrix_container.R
