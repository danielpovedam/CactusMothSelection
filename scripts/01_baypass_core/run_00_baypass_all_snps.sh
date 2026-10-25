#!/bin/bash
#SBATCH --job-name=baypass_all_snps
#SBATCH --output=logs/baypass_all_snps_%A_%a.out
#SBATCH --error=logs/baypass_all_snps_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-2

set -euo pipefail

# ============================================================
# BayPass v3.1
# Final all-SNP analysis for Cactoblastis cactorum
#
# Analyses:
#   - core model / population covariance
#   - XtX differentiation
#   - genotype-environment associations
#   - C2 host-use contrast
#
# Three independent MCMC chains.
# ============================================================

module purge
module load bioinfo/BayPass/3.1


# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

PROJECT="/home/dpoveda/micetral/Population_snps"

BAYPASS_DIR="${PROJECT}/baypass"
PREP_DIR="${PROJECT}/baypass_preparation"

GFILE="${PREP_DIR}/baypass/all_snps.baypass.geno"
COVFILE="${BAYPASS_DIR}/covariates/all.variables.cov"
C2FILE="${BAYPASS_DIR}/covariates/Host_C2.cov"

OUTBASE="${BAYPASS_DIR}/runs_10x/all_snps"

mkdir -p "${OUTBASE}" logs


# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

NPOP=28

SEEDS=(
    12345
    67890
    13579
)

NPILOT=20
PILOTLENGTH=500
BURNIN=10000
NVAL=2500
THIN=20


# ------------------------------------------------------------
# SELECT MCMC CHAIN
# ------------------------------------------------------------

SEED="${SEEDS[$SLURM_ARRAY_TASK_ID]}"

OUTDIR="${OUTBASE}/seed_${SEED}"
mkdir -p "${OUTDIR}"

PREFIX="${OUTDIR}/BP_all_snps_seed${SEED}"


# ------------------------------------------------------------
# INPUT CHECKS
# ------------------------------------------------------------

[[ -s "${GFILE}" ]] || {
    echo "ERROR: genotype file not found: ${GFILE}"
    exit 1
}

[[ -s "${COVFILE}" ]] || {
    echo "ERROR: covariate file not found: ${COVFILE}"
    exit 1
}

[[ -s "${C2FILE}" ]] || {
    echo "ERROR: C2 contrast file not found: ${C2FILE}"
    exit 1
}


# ------------------------------------------------------------
# RUN BAYPASS
# ------------------------------------------------------------

baypass \
    -countdatafile "${GFILE}" \
    -npop "${NPOP}" \
    -efile "${COVFILE}" \
    -contrastfile "${C2FILE}" \
    -npilot "${NPILOT}" \
    -pilotlength "${PILOTLENGTH}" \
    -burnin "${BURNIN}" \
    -nval "${NVAL}" \
    -thin "${THIN}" \
    -print_omega_samples \
    -seed "${SEED}" \
    -nthreads "${SLURM_CPUS_PER_TASK}" \
    -outprefix "${PREFIX}"

echo "BayPass completed successfully:"
echo "${OUTDIR}"
