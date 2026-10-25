#!/bin/bash
#SBATCH --job-name=baypass_cacto
#SBATCH --output=logs/baypass_%A_%a.out
#SBATCH --error=logs/baypass_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-5

set -euo pipefail


# ============================================================
# BayPass v3.1
# Cactoblastis cactorum final analysis
#
# Analyses:
#   - core model / population covariance
#   - XtX differentiation
#   - genotype-environment associations
#   - C2 host-use contrast
#
# Datasets:
#   - LD-pruned SNP dataset ("random_snp")
#   - all-SNP dataset ("all_snps")
#
# Three independent MCMC chains per dataset.
# ============================================================


# ------------------------------------------------------------
# SOFTWARE
# ------------------------------------------------------------

module purge
module load bioinfo/BayPass/3.1


# ------------------------------------------------------------
# PROJECT PATHS
#
# Change PROJECT if the repository/data are stored elsewhere.
# ------------------------------------------------------------

PROJECT="/home/dpoveda/micetral/Population_snps"

BAYPASS_DIR="${PROJECT}/baypass"
PREP_DIR="${PROJECT}/baypass_preparation"

OUTBASE="${BAYPASS_DIR}/runs_10x"
LOGDIR="${BAYPASS_DIR}/logs"

mkdir -p "${OUTBASE}" "${LOGDIR}"


# ------------------------------------------------------------
# POPULATIONS
# ------------------------------------------------------------

NPOP=28


# ------------------------------------------------------------
# DATASETS
# ------------------------------------------------------------

DATASETS=(
    "random_snp"
    "all_snps"
)

GFILES=(
    "${PREP_DIR}/baypass/random_snp.baypass.geno"
    "${PREP_DIR}/baypass/all_snps.baypass.geno"
)


# ------------------------------------------------------------
# ENVIRONMENTAL COVARIATES
#
# Rows:
#   1 TempPC1
#   2 PrecipPC1
#   3 AridityIndex
#   4 Elevation
#   5 Latitude
#   6 Longitude
#   7 BIO1
#   8 BIO12
#
# Columns correspond to the 28 populations in BayPass order.
# ------------------------------------------------------------

COVFILE="${BAYPASS_DIR}/covariates/all.variables.cov"


# ------------------------------------------------------------
# HOST-USE CONTRAST
#
# O. ficus-indica = +1
# Native Opuntia  = -1
# ------------------------------------------------------------

C2FILE="${BAYPASS_DIR}/covariates/Host_C2.cov"


# ------------------------------------------------------------
# INDEPENDENT MCMC CHAINS
# ------------------------------------------------------------

SEEDS=(
    12345
    67890
    13579
)


# ------------------------------------------------------------
# FINAL MCMC SETTINGS
#
# Pilot phase:    20 x 500 = 10,000 iterations
# Burn-in:        10,000 iterations
# Posterior:      2,500 retained values
# Thinning:       20
# ------------------------------------------------------------

NPILOT=20
PILOTLENGTH=500
BURNIN=10000
NVAL=2500
THIN=20


# ------------------------------------------------------------
# SLURM ARRAY MAPPING
#
# 0 = random_snp / 12345
# 1 = random_snp / 67890
# 2 = random_snp / 13579
# 3 = all_snps   / 12345
# 4 = all_snps   / 67890
# 5 = all_snps   / 13579
# ------------------------------------------------------------

NSEED=${#SEEDS[@]}

DATASET_IDX=$(( SLURM_ARRAY_TASK_ID / NSEED ))
SEED_IDX=$(( SLURM_ARRAY_TASK_ID % NSEED ))

DATASET="${DATASETS[$DATASET_IDX]}"
GFILE="${GFILES[$DATASET_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"


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
# OUTPUT
# ------------------------------------------------------------

OUTDIR="${OUTBASE}/${DATASET}/seed_${SEED}"

mkdir -p "${OUTDIR}"

PREFIX="${OUTDIR}/BP_${DATASET}_seed${SEED}"


# ------------------------------------------------------------
# REPORT
# ------------------------------------------------------------

NSNP=$(wc -l < "${GFILE}")

echo "======================================================"
echo "BayPass v3.1 - Cactoblastis cactorum"
echo "======================================================"
echo "Dataset:      ${DATASET}"
echo "Seed:         ${SEED}"
echo "Populations:  ${NPOP}"
echo "SNPs:         ${NSNP}"
echo "Threads:      ${SLURM_CPUS_PER_TASK}"
echo
echo "MCMC:"
echo "  npilot       = ${NPILOT}"
echo "  pilotlength  = ${PILOTLENGTH}"
echo "  burnin       = ${BURNIN}"
echo "  nval         = ${NVAL}"
echo "  thin         = ${THIN}"
echo
echo "Output:"
echo "  ${PREFIX}"
echo "======================================================"


# ------------------------------------------------------------
# RUN BAYPASS
#
# Covariate scaling is performed by default in BayPass v3.1.
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


echo
echo "BayPass completed successfully:"
echo "${OUTDIR}"
