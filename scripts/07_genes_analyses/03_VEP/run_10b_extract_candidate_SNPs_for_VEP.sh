#!/bin/bash
#SBATCH --job-name=extract_vep_ccac
#SBATCH --output=logs/extract_vep_ccac_%j.out
#SBATCH --error=logs/extract_vep_ccac_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail

mkdir -p logs

module purge

############################################################
## EDIT ONLY THIS PATH IF NEEDED
############################################################

VCF_ALL="/home/dpoveda/micetral/Population_snps/final_vcf_filtering/all_snps/all_snps.final.filtered.vcf.gz"

############################################################

OUTDIR="/home/dpoveda/micetral/Population_snps/baypass/selection_scan/VEP_ALL_selection_candidates"

REGIONS="${OUTDIR}/ALL_selection_candidate_SNP_positions.tsv"

VCF_OUT="${OUTDIR}/ALL_selection_candidate_SNPs_for_VEP.vcf.gz"

echo "============================================================"
echo "Extract candidate SNPs for VEP"
echo "Date: $(date)"
echo "Input VCF: ${VCF_ALL}"
echo "Regions: ${REGIONS}"
echo "Output: ${VCF_OUT}"
echo "============================================================"

[[ -s "${VCF_ALL}" ]] || {
    echo "ERROR: missing VCF ${VCF_ALL}"
    exit 1
}

[[ -s "${REGIONS}" ]] || {
    echo "ERROR: missing regions ${REGIONS}"
    exit 1
}


# ============================================================
# CHECK / CREATE VCF INDEX
# ============================================================

if [[ ! -s "${VCF_ALL}.tbi" && ! -s "${VCF_ALL}.csi" ]]; then

    echo "Indexing source VCF..."

    bcftools index \
        --threads "${SLURM_CPUS_PER_TASK}" \
        "${VCF_ALL}"
fi


# ============================================================
# EXTRACT EXACT CANDIDATE POSITIONS
# ============================================================

bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -R "${REGIONS}" \
    -m2 \
    -M2 \
    -v snps \
    -Oz \
    -o "${VCF_OUT}" \
    "${VCF_ALL}"


bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    "${VCF_OUT}"


# ============================================================
# QC
# ============================================================

echo
echo "Candidate genomic positions requested:"
wc -l "${REGIONS}"

echo
echo "VCF SNP records extracted:"
bcftools view -H "${VCF_OUT}" | wc -l

echo
echo "VCF chromosomes:"
bcftools query -f '%CHROM\n' "${VCF_OUT}" \
    | sort -u \
    | wc -l

echo
echo "DONE"
