#!/bin/bash
#SBATCH --job-name=bp_C2geo16
#SBATCH --output=sensitivity_C2_geographic/logs/bp_C2geo16_%A_%a.out
#SBATCH --error=sensitivity_C2_geographic/logs/bp_C2geo16_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-2

set -euo pipefail

############################################
# MODULE
############################################

module purge
module load bioinfo/BayPass/3.1

############################################
# CONFIG
############################################

BASE="/home/dpoveda/micetral/Population_snps/baypass"
ROOT="${BASE}/sensitivity_C2_geographic"

GFILE="${ROOT}/input/all_snps.C2_geographic_16pop.polymorphic_only.geno"

# 41A created this as 16 rows x 1 column.
# BayPass v3.1 expects the C2 contrast as 1 row x NPOP columns.
C2_VERTICAL="${ROOT}/input/Host_C2_geographic_16pop.txt"
C2FILE="${ROOT}/input/Host_C2_geographic_16pop.cov"

OUTBASE="${ROOT}/runs"
LOGDIR="${ROOT}/logs"

NPOP=16
EXPECTED_SNPS=130732

############################################
# THREE INDEPENDENT MCMC CHAINS
############################################

SEEDS=(
    12345
    67890
    13579
)

############################################
# SAME MCMC SETTINGS AS PRIMARY ANALYSIS
############################################

NPILOT=20
PILOTLENGTH=500
BURNIN=10000
NVAL=2500
THIN=20

############################################
# CREATE DIRECTORIES
############################################

mkdir -p "${OUTBASE}" "${LOGDIR}"

############################################
# BUILD BAYPASS C2 FILE
#
# Convert:
#   16 rows x 1 column
# to:
#   1 row x 16 columns
############################################

[[ -s "${C2_VERTICAL}" ]] || {
    echo "[ERROR] Missing host contrast source:"
    echo "${C2_VERTICAL}"
    exit 1
}

paste -sd' ' "${C2_VERTICAL}" > "${C2FILE}"

############################################
# TASK MAPPING
############################################

SEED="${SEEDS[$SLURM_ARRAY_TASK_ID]}"

OUTDIR="${OUTBASE}/seed_${SEED}"
mkdir -p "${OUTDIR}"

PREFIX="${OUTDIR}/BP_C2geo16_seed${SEED}"

############################################
# INPUT CHECKS
############################################

[[ -s "${GFILE}" ]] || {
    echo "[ERROR] Missing genotype file:"
    echo "${GFILE}"
    exit 1
}

[[ -s "${C2FILE}" ]] || {
    echo "[ERROR] Missing C2 file:"
    echo "${C2FILE}"
    exit 1
}

############################################
# VALIDATE GENOTYPE FILE
############################################

NSNP=$(wc -l < "${GFILE}")
GENO_NCOL=$(awk 'NR==1 {print NF}' "${GFILE}")
EXPECTED_GENO_COLS=$((2 * NPOP))

if [[ "${NSNP}" -ne "${EXPECTED_SNPS}" ]]; then
    echo "[ERROR] Genotype file has ${NSNP} SNP rows."
    echo "Expected ${EXPECTED_SNPS}."
    exit 1
fi

if [[ "${GENO_NCOL}" -ne "${EXPECTED_GENO_COLS}" ]]; then
    echo "[ERROR] Genotype file has ${GENO_NCOL} columns."
    echo "Expected ${EXPECTED_GENO_COLS} columns for ${NPOP} populations."
    exit 1
fi

############################################
# VALIDATE GENOTYPE VALUES
############################################

awk '
{
    for(i=1;i<=NF;i++) {
        if($i !~ /^[0-9]+$/) {
            print "[ERROR] Non-integer genotype count:",
                  "row=" NR,
                  "column=" i,
                  "value=" $i > "/dev/stderr"
            bad=1
            exit
        }
    }
}
END {
    if(bad)
        exit 1
}
' "${GFILE}"

############################################
# VALIDATE C2 FILE
############################################

C2_ROWS=$(wc -l < "${C2FILE}")
C2_COLS=$(awk 'NR==1 {print NF}' "${C2FILE}")

if [[ "${C2_ROWS}" -ne 1 ]]; then
    echo "[ERROR] Host C2 file must contain exactly one row."
    echo "Found ${C2_ROWS} rows."
    exit 1
fi

if [[ "${C2_COLS}" -ne "${NPOP}" ]]; then
    echo "[ERROR] Host C2 file contains ${C2_COLS} values."
    echo "Expected ${NPOP}."
    exit 1
fi

awk '
{
    nneg=0
    npos=0

    for(i=1;i<=NF;i++) {
        if($i != -1 && $i != 1) {
            print "[ERROR] Host C2 value must be -1 or 1:",
                  "column=" i,
                  "value=" $i > "/dev/stderr"
            bad=1
        }

        if($i == -1) nneg++
        if($i == 1)  npos++
    }

    if(nneg != 8 || npos != 8) {
        print "[ERROR] Expected 8 native (-1) and 8 exotic (+1).",
              "Observed:", nneg, "native and", npos, "exotic." > "/dev/stderr"
        bad=1
    }
}
END {
    if(bad)
        exit 1
}
' "${C2FILE}"

############################################
# REPORT SETTINGS
############################################

POSTERIOR_ITER=$((NVAL * THIN))
PILOT_ITER=$((NPILOT * PILOTLENGTH))

echo "======================================================"
echo "BayPass 3.1 — Cactoblastis cactorum"
echo "GEOGRAPHICALLY RESTRICTED HOST C2 SENSITIVITY"
echo "CORE + XtX + C2; NO GEA"
echo "======================================================"
echo
echo "SLURM task:       ${SLURM_ARRAY_TASK_ID}"
echo "Seed:             ${SEED}"
echo
echo "N populations:    ${NPOP}"
echo "N SNPs:           ${NSNP}"
echo "Geno columns:     ${GENO_NCOL}"
echo "Host groups:      8 O. ficus-indica / 8 native Opuntia"
echo
echo "MCMC SETTINGS"
echo "--------------------------------------------"
echo "npilot:           ${NPILOT}"
echo "pilotlength:      ${PILOTLENGTH}"
echo "pilot iterations: ${PILOT_ITER}"
echo "burnin:           ${BURNIN}"
echo "nval:             ${NVAL}"
echo "thin:             ${THIN}"
echo "sampling iter.:   ${POSTERIOR_ITER}"
echo
echo "Threads:          ${SLURM_CPUS_PER_TASK}"
echo
echo "GFILE:"
echo "${GFILE}"
echo
echo "C2FILE:"
echo "${C2FILE}"
echo
echo "C2 contrast:"
cat "${C2FILE}"
echo
echo "OUTPUT PREFIX:"
echo "${PREFIX}"
echo
echo "======================================================"

############################################
# RUN BAYPASS
#
# Same BayPass v3.1 command style as the
# primary full-range analysis.
#
# No -efile: this sensitivity analysis runs
# only the core model + XtX + C2 contrast.
############################################

baypass \
    -countdatafile "${GFILE}" \
    -npop "${NPOP}" \
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

############################################
# VERIFY CORE OUTPUTS
############################################

OMEGA="${PREFIX}_mat_omega.out"
XTX="${PREFIX}_summary_pi_xtx.out"
BETA="${PREFIX}_summary_beta_params.out"
BPLOG="${PREFIX}_baypass.log"

echo
echo "======================================================"
echo "OUTPUT VALIDATION"
echo "======================================================"

for f in \
    "${OMEGA}" \
    "${XTX}" \
    "${BETA}" \
    "${BPLOG}"
do
    if [[ ! -s "${f}" ]]; then
        echo "[ERROR] Expected output missing:"
        echo "${f}"
        exit 1
    fi

    echo "[OK] $(basename "${f}")"
done

############################################
# LOCATE C2 OUTPUT
#
# Keep this flexible because BayPass output
# suffixes can differ slightly by build.
############################################

C2OUT=$(find "${OUTDIR}" -maxdepth 1 -type f \
    \( -name "${PREFIX##*/}_summary_contrast.out" \
       -o -name "*summary_contrast*.out" \
       -o -name "*contrast*.out" \) \
    | head -n 1 || true)

if [[ -z "${C2OUT}" || ! -s "${C2OUT}" ]]; then
    echo "[ERROR] Could not locate the C2 contrast output."
    echo "Files produced in ${OUTDIR}:"
    ls -lh "${OUTDIR}"
    exit 1
fi

echo "[OK] C2 output: $(basename "${C2OUT}")"

############################################
# OMEGA DIMENSIONS
############################################

OMEGA_ROWS=$(wc -l < "${OMEGA}")
OMEGA_COLS=$(awk 'NR==1 {print NF}' "${OMEGA}")

if [[ "${OMEGA_ROWS}" -ne "${NPOP}" ||
      "${OMEGA_COLS}" -ne "${NPOP}" ]]
then
    echo "[ERROR] Omega matrix dimensions are incorrect:"
    echo "${OMEGA_ROWS} x ${OMEGA_COLS}"
    exit 1
fi

############################################
# XtX SNP COUNT
############################################

XTX_ROWS=$(awk 'END {print NR-1}' "${XTX}")

if [[ "${XTX_ROWS}" -ne "${NSNP}" ]]; then
    echo "[ERROR] XtX output contains ${XTX_ROWS} SNPs."
    echo "Expected ${NSNP}."
    exit 1
fi

############################################
# FINAL MESSAGE
############################################

echo
echo "======================================================"
echo "SUCCESS"
echo "======================================================"
echo
echo "Seed:           ${SEED}"
echo "SNPs:           ${NSNP}"
echo "Omega:          ${OMEGA_ROWS} x ${OMEGA_COLS}"
echo "XtX SNP rows:   ${XTX_ROWS}"
echo "C2 output:      ${C2OUT}"
echo
echo "Output:"
echo "${OUTDIR}"
echo
echo "======================================================"
