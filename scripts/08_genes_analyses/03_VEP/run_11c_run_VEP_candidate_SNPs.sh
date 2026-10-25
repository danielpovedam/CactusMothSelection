#!/bin/bash
#SBATCH --job-name=vep_ccac_final
#SBATCH --output=logs/vep_ccac_final_%j.out
#SBATCH --error=logs/vep_ccac_final_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

set -euo pipefail

############################################################
## FINAL VEP ANALYSIS
## Cactoblastis cactorum
##
## Candidate SNP set:
##   XtX + host-use C2 + GEA
##
## Workflow:
##   1. Load software
##   2. Identify final all-SNP VCF
##   3. Extract candidate SNP positions
##   4. Verify complete position recovery
##   5. Filter corrected GFF to FASTA-compatible contigs
##   6. Compress/index GFF
##   7. Sort/index candidate VCF
##   8. Run one-SNP VEP smoke test
##   9. Run full VEP
##  10. Generate consequence and QC summaries
##
## IMPORTANT:
##   - No Ensembl cache is used.
##   - No --offline is used.
##   - Annotation source = custom GFF + custom FASTA.
############################################################


# ============================================================
# SOFTWARE
# ============================================================

mkdir -p logs

module purge

module load devel/Miniconda/Miniconda3
module load bioinfo/Ensembl-VEP/111.0
module load bioinfo/samtools/1.23
module load bioinfo/Bcftools/1.21

source activate bio-env

export LC_ALL=C
export LANG=C


# ============================================================
# SOFTWARE QC
# ============================================================

for prog in vep samtools bcftools bgzip tabix python
do
    if ! command -v "${prog}" >/dev/null 2>&1
    then
        echo "ERROR: required executable not found: ${prog}"
        exit 1
    fi
done


echo "============================================================"
echo "SOFTWARE QC"
echo "============================================================"

echo "vep:      $(command -v vep)"
echo "samtools: $(command -v samtools)"
echo "bcftools: $(command -v bcftools)"
echo "bgzip:    $(command -v bgzip)"
echo "tabix:    $(command -v tabix)"
echo "python:   $(command -v python)"

echo
samtools --version | head -1
bcftools --version | head -1

echo "============================================================"


# ============================================================
# PATHS
# ============================================================

ROOT="/home/dpoveda/micetral/Population_snps"

BASE="${ROOT}/baypass"

OUTDIR="${BASE}/selection_scan/VEP_ALL_selection_candidates"

mkdir -p "${OUTDIR}"


# ============================================================
# REFERENCE
# ============================================================

FASTA="${ROOT}/Ccactorum.fcs_gx_2026_FINAL.fa"
FAI="${FASTA}.fai"


# ============================================================
# CORRECTED GFF
# ============================================================

GFF_ORIGINAL="${ROOT}/Ccac_functional-annotation/02_agat/Ccac/out_abril4/Ccac.annotations.corrected.gff3"


# ============================================================
# FINAL FILTERING SUMMARY
# ============================================================

FILTER_SUMMARY="${ROOT}/final_vcf_filtering/Table_S3_vcf_filtering_summary.tsv"


# ============================================================
# CANDIDATE INPUTS
# ============================================================

MASTER="${OUTDIR}/ALL_selection_candidate_SNPs_master.tsv"

REGIONS="${OUTDIR}/ALL_selection_candidate_SNP_positions.tsv"


# ============================================================
# OUTPUT VCFs
# ============================================================

VCF_EXTRACTED="${OUTDIR}/ALL_selection_candidate_SNPs_for_VEP.vcf.gz"

VCF_SORTED="${OUTDIR}/ALL_selection_candidate_SNPs_for_VEP.sorted.vcf.gz"


# ============================================================
# VEP-SPECIFIC GFF
# ============================================================

GFF_VEP="${OUTDIR}/Ccac.annotations.corrected.FASTAcompatible.sorted.gff3"

GFF_VEP_GZ="${GFF_VEP}.gz"


# ============================================================
# VEP OUTPUT
# ============================================================

VEP_OUT="${OUTDIR}/ALL_selection_candidate_SNPs.vep.tsv"


# ============================================================
# QC OUTPUTS
# ============================================================

TARGETS="${OUTDIR}/candidate_targets.tsv"

EXTRACTED_POS="${OUTDIR}/extracted_unique_positions.tsv"

MISSING_POS="${OUTDIR}/missing_candidate_positions.tsv"

GFF_QC="${OUTDIR}/GFF_FASTA_compatibility_QC.txt"

RAW_COUNTS="${OUTDIR}/VEP_raw_consequence_counts.txt"

QC_SUMMARY="${OUTDIR}/VEP_run_QC_summary.txt"

VCF_INFO="${OUTDIR}/selected_final_allSNP_VCF.txt"


# ============================================================
# SMOKE TEST OUTPUTS
# ============================================================

TEST_VCF="${OUTDIR}/VEP_test_one_variant.vcf.gz"

VEP_TEST_OUT="${OUTDIR}/VEP_test_one_variant.tsv"

VEP_TEST_ERR="${OUTDIR}/VEP_test.stderr.txt"


# ============================================================
# HEADER
# ============================================================

echo
echo "============================================================"
echo "FINAL VEP ANALYSIS - Cactoblastis cactorum"
echo "============================================================"

echo "Date: $(date)"
echo "Node: $(hostname)"

echo
echo "Reference FASTA:"
echo "${FASTA}"

echo
echo "Corrected GFF:"
echo "${GFF_ORIGINAL}"

echo
echo "Candidate positions:"
echo "${REGIONS}"

echo
echo "Filtering summary:"
echo "${FILTER_SUMMARY}"

echo
echo "Output directory:"
echo "${OUTDIR}"

echo "============================================================"


# ============================================================
# INPUT CHECKS
# ============================================================

for file in \
    "${FASTA}" \
    "${GFF_ORIGINAL}" \
    "${FILTER_SUMMARY}" \
    "${REGIONS}"
do

    if [[ ! -s "${file}" ]]
    then
        echo
        echo "ERROR: missing or empty input file:"
        echo "${file}"
        exit 1
    fi

done


if [[ -s "${MASTER}" ]]
then
    echo
    echo "Candidate SNP master table found:"
    echo "${MASTER}"
else
    echo
    echo "WARNING: candidate SNP master table not found:"
    echo "${MASTER}"
fi


# ============================================================
# FASTA INDEX
# ============================================================

echo
echo "============================================================"
echo "REFERENCE INDEX"
echo "============================================================"


if [[ ! -s "${FAI}" ]]
then

    echo "FASTA index absent. Creating..."
    samtools faidx "${FASTA}"

else

    echo "FASTA index already exists."

fi


if [[ ! -s "${FAI}" ]]
then
    echo "ERROR: FASTA index was not created."
    exit 1
fi


N_FASTA_CONTIGS=$(wc -l < "${FAI}")

echo "FASTA contigs: ${N_FASTA_CONTIGS}"


# ============================================================
# IDENTIFY FINAL ALL-SNP VCF
# ============================================================

echo
echo "============================================================"
echo "IDENTIFYING FINAL ALL-SNP VCF"
echo "============================================================"


python - "${FILTER_SUMMARY}" "${VCF_INFO}" <<'PY'
import csv
import sys

table = sys.argv[1]
outfile = sys.argv[2]

hits = []

with open(table, newline="") as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    required = {
        "Dataset",
        "Step",
        "VCF_file"
    }

    missing = required.difference(
        reader.fieldnames or []
    )

    if missing:
        raise RuntimeError(
            "Missing required columns: "
            + ", ".join(sorted(missing))
        )

    for row in reader:

        if row["Dataset"].strip() != "all_snps":
            continue

        try:
            step = int(float(row["Step"]))
        except Exception:
            continue

        vcf = row["VCF_file"].strip()

        if not vcf:
            continue

        hits.append(
            (
                step,
                vcf,
                row.get("Filtering_step", ""),
                row.get("N_SNPs_retained", "")
            )
        )


if not hits:
    raise RuntimeError(
        "No valid Dataset=all_snps VCF found"
    )


hits.sort(
    key=lambda x: x[0],
    reverse=True
)


step, vcf, filtering, nsnps = hits[0]


with open(outfile, "w") as out:
    out.write(vcf + "\n")


print(
    f"Dataset = all_snps | "
    f"Step = {step} | "
    f"Filtering_step = {filtering} | "
    f"N_SNPs_retained = {nsnps}"
)

print(
    f"VCF = {vcf}"
)
PY


VCF_ALL=$(head -1 "${VCF_INFO}")


if [[ -z "${VCF_ALL}" ]]
then
    echo "ERROR: could not identify final all-SNP VCF."
    exit 1
fi


if [[ ! -s "${VCF_ALL}" ]]
then
    echo "ERROR: final all-SNP VCF does not exist:"
    echo "${VCF_ALL}"
    exit 1
fi


echo
echo "Final all-SNP VCF:"
echo "${VCF_ALL}"


# ============================================================
# SOURCE VCF INDEX
# ============================================================

echo
echo "============================================================"
echo "SOURCE VCF INDEX"
echo "============================================================"


if [[ ! -s "${VCF_ALL}.tbi" && ! -s "${VCF_ALL}.csi" ]]
then

    echo "Creating source VCF index..."

    bcftools index \
        --threads "${SLURM_CPUS_PER_TASK}" \
        "${VCF_ALL}"

else

    echo "Source VCF index already exists."

fi


# ============================================================
# PREPARE CANDIDATE TARGETS
# ============================================================

echo
echo "============================================================"
echo "CANDIDATE POSITION QC"
echo "============================================================"


python - "${REGIONS}" "${TARGETS}" <<'PY'
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

positions = set()

with open(input_file) as handle:

    for line_number, line in enumerate(handle, start=1):

        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 2:
            raise RuntimeError(
                f"Invalid candidate line {line_number}"
            )

        chrom = fields[0].strip()

        try:
            pos = int(fields[1])
        except ValueError:
            raise RuntimeError(
                f"Invalid position at line {line_number}"
            )

        positions.add(
            (
                chrom,
                pos
            )
        )


if not positions:
    raise RuntimeError(
        "No candidate positions found"
    )


with open(output_file, "w") as out:

    for chrom, pos in sorted(
        positions,
        key=lambda x: (
            x[0],
            x[1]
        )
    ):

        out.write(
            f"{chrom}\t{pos}\t{pos}\n"
        )


print(
    "Unique requested positions:",
    len(positions)
)
PY


N_REQUESTED=$(wc -l < "${TARGETS}")

echo "Unique candidate positions: ${N_REQUESTED}"


# ============================================================
# EXTRACT CANDIDATE SNPs
# ============================================================

echo
echo "============================================================"
echo "EXTRACTING CANDIDATE SNPs"
echo "============================================================"


rm -f \
    "${VCF_EXTRACTED}" \
    "${VCF_EXTRACTED}.tbi" \
    "${VCF_EXTRACTED}.csi"


bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -T "${TARGETS}" \
    -m2 \
    -M2 \
    -v snps \
    -Oz \
    -o "${VCF_EXTRACTED}" \
    "${VCF_ALL}"


if [[ ! -s "${VCF_EXTRACTED}" ]]
then
    echo "ERROR: candidate VCF was not generated."
    exit 1
fi


bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    "${VCF_EXTRACTED}"


N_RECORDS=$(bcftools view -H "${VCF_EXTRACTED}" | wc -l)


if [[ "${N_RECORDS}" -eq 0 ]]
then
    echo "ERROR: zero candidate SNP records recovered."
    exit 1
fi


echo "Extracted candidate VCF records: ${N_RECORDS}"


# ============================================================
# POSITION RECOVERY QC
# ============================================================

echo
echo "============================================================"
echo "POSITION RECOVERY QC"
echo "============================================================"


bcftools query \
    -f '%CHROM\t%POS\n' \
    "${VCF_EXTRACTED}" \
    | sort -u \
    > "${EXTRACTED_POS}"


N_EXTRACTED_POS=$(wc -l < "${EXTRACTED_POS}")


python - \
    "${TARGETS}" \
    "${EXTRACTED_POS}" \
    "${MISSING_POS}" <<'PY'
import sys


def read_positions(path):

    result = set()

    with open(path) as handle:

        for line in handle:

            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            result.add(
                (
                    fields[0],
                    int(fields[1])
                )
            )

    return result


requested = read_positions(
    sys.argv[1]
)

observed = read_positions(
    sys.argv[2]
)

missing = requested - observed


with open(sys.argv[3], "w") as out:

    for chrom, pos in sorted(
        missing,
        key=lambda x: (
            x[0],
            x[1]
        )
    ):

        out.write(
            f"{chrom}\t{pos}\n"
        )


print(
    "Requested:",
    len(requested)
)

print(
    "Recovered:",
    len(observed)
)

print(
    "Missing:",
    len(missing)
)


if missing:
    raise RuntimeError(
        "Not all candidate positions were recovered"
    )
PY


N_MISSING=$(wc -l < "${MISSING_POS}")


echo
echo "Unique requested positions: ${N_REQUESTED}"
echo "Unique recovered positions: ${N_EXTRACTED_POS}"
echo "Missing positions:          ${N_MISSING}"


# ============================================================
# SNP-ONLY QC
# ============================================================

N_NON_SNP=$(
    bcftools view \
        -H \
        -e 'TYPE="snp"' \
        "${VCF_EXTRACTED}" \
    | wc -l
)


echo "Non-SNP records: ${N_NON_SNP}"


if [[ "${N_NON_SNP}" -ne 0 ]]
then
    echo "ERROR: non-SNP records detected."
    exit 1
fi


# ============================================================
# PREPARE FASTA-COMPATIBLE GFF
# ============================================================

echo
echo "============================================================"
echo "PREPARING FASTA-COMPATIBLE GFF"
echo "============================================================"


rm -f \
    "${GFF_VEP}" \
    "${GFF_VEP_GZ}" \
    "${GFF_VEP_GZ}.tbi"


python - \
    "${GFF_ORIGINAL}" \
    "${FAI}" \
    "${GFF_VEP}" \
    "${GFF_QC}" <<'PY'

import sys
from collections import Counter


gff_in = sys.argv[1]
fai = sys.argv[2]
gff_out = sys.argv[3]
qc_out = sys.argv[4]


# ------------------------------------------------------------
# FASTA contigs and order
# ------------------------------------------------------------

order = {}

with open(fai) as handle:

    for i, line in enumerate(handle):

        chrom = line.rstrip("\n").split("\t")[0]

        order[chrom] = i


# ------------------------------------------------------------
# Read and filter GFF
# ------------------------------------------------------------

features = []

removed_contigs = Counter()

retained_types = Counter()

removed_features = 0

retained_features = 0

original_index = 0

inside_fasta = False


with open(gff_in) as handle:

    for raw in handle:


        if raw.startswith("##FASTA"):

            inside_fasta = True
            continue


        if inside_fasta:

            continue


        if raw.startswith("#"):

            continue


        line = raw.rstrip("\n")


        if not line:

            continue


        fields = line.split("\t")


        if len(fields) < 9:

            continue


        chrom = fields[0]

        feature_type = fields[2]


        if chrom not in order:

            removed_contigs[chrom] += 1

            removed_features += 1

            continue


        start = int(fields[3])

        end = int(fields[4])


        features.append(
            (
                order[chrom],
                start,
                end,
                original_index,
                line
            )
        )


        retained_types[
            feature_type
        ] += 1


        retained_features += 1

        original_index += 1


# ------------------------------------------------------------
# Sort
# ------------------------------------------------------------

features.sort(
    key=lambda x: (
        x[0],
        x[1],
        x[2],
        x[3]
    )
)


# ------------------------------------------------------------
# Write clean GFF
# ------------------------------------------------------------

with open(gff_out, "w") as out:

    for _, _, _, _, line in features:

        out.write(
            line + "\n"
        )


# ------------------------------------------------------------
# QC report
# ------------------------------------------------------------

with open(qc_out, "w") as out:

    out.write(
        f"FASTA_contigs\t{len(order)}\n"
    )

    out.write(
        f"GFF_features_retained\t{retained_features}\n"
    )

    out.write(
        f"GFF_features_removed\t{removed_features}\n"
    )

    out.write(
        f"GFF_contigs_removed\t{len(removed_contigs)}\n"
    )


    out.write(
        "\nRemoved_contigs\n"
    )


    for chrom, n in sorted(
        removed_contigs.items()
    ):

        out.write(
            f"{chrom}\t{n}\n"
        )


    out.write(
        "\nRetained_feature_types\n"
    )


    for feature, n in retained_types.most_common():

        out.write(
            f"{feature}\t{n}\n"
        )


print(
    "FASTA contigs:",
    len(order)
)

print(
    "GFF retained features:",
    retained_features
)

print(
    "GFF removed features:",
    removed_features
)

print(
    "GFF removed contigs:",
    len(removed_contigs)
)


important = [
    "gene",
    "transcript",
    "mRNA",
    "exon",
    "CDS",
    "intron",
    "five_prime_UTR",
    "three_prime_UTR"
]


print()

print(
    "Retained important feature types:"
)


for feature in important:

    print(
        feature,
        retained_types.get(
            feature,
            0
        )
    )


n_transcripts = (
    retained_types.get(
        "transcript",
        0
    )
    +
    retained_types.get(
        "mRNA",
        0
    )
)


if n_transcripts == 0:

    raise RuntimeError(
        "No transcript/mRNA features retained"
    )


if retained_types.get(
    "exon",
    0
) == 0:

    raise RuntimeError(
        "No exon features retained"
    )


if retained_types.get(
    "CDS",
    0
) == 0:

    raise RuntimeError(
        "No CDS features retained"
    )
PY


if [[ ! -s "${GFF_VEP}" ]]
then
    echo "ERROR: FASTA-compatible GFF was not created."
    exit 1
fi


# ============================================================
# BGZIP + TABIX GFF
# ============================================================

echo
echo "Compressing GFF..."


bgzip \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -c "${GFF_VEP}" \
    > "${GFF_VEP_GZ}"


echo "Indexing GFF..."


tabix \
    -f \
    -p gff \
    "${GFF_VEP_GZ}"


if [[ ! -s "${GFF_VEP_GZ}" ]]
then
    echo "ERROR: compressed GFF missing."
    exit 1
fi


if [[ ! -s "${GFF_VEP_GZ}.tbi" ]]
then
    echo "ERROR: GFF tabix index missing."
    exit 1
fi


# ============================================================
# VCF / FASTA COMPATIBILITY
# ============================================================

echo
echo "============================================================"
echo "VCF / FASTA COMPATIBILITY"
echo "============================================================"


python - "${VCF_EXTRACTED}" "${FAI}" <<'PY'

import gzip
import sys


vcf = sys.argv[1]
fai = sys.argv[2]


fasta_contigs = set()


with open(fai) as handle:

    for line in handle:

        fasta_contigs.add(
            line.split("\t")[0]
        )


vcf_contigs = set()


with gzip.open(
    vcf,
    "rt"
) as handle:

    for line in handle:

        if line.startswith("#"):
            continue

        vcf_contigs.add(
            line.split("\t")[0]
        )


bad = sorted(
    vcf_contigs -
    fasta_contigs
)


print(
    "Candidate VCF contigs:",
    len(vcf_contigs)
)

print(
    "Candidate VCF contigs absent from FASTA:",
    len(bad)
)


if bad:

    print(
        "\n".join(bad)
    )

    raise RuntimeError(
        "Candidate VCF contains contigs absent from FASTA"
    )
PY


# ============================================================
# SORT CANDIDATE VCF
#
# IMPORTANT:
# bcftools sort -Oz creates BGZF-compressed output.
# ============================================================

echo
echo "============================================================"
echo "SORTING CANDIDATE VCF"
echo "============================================================"


rm -f \
    "${VCF_SORTED}" \
    "${VCF_SORTED}.tbi" \
    "${VCF_SORTED}.csi"


bcftools sort \
    "${VCF_EXTRACTED}" \
    -Oz \
    -o "${VCF_SORTED}"


if [[ ! -s "${VCF_SORTED}" ]]
then
    echo "ERROR: sorted VCF was not produced."
    exit 1
fi


N_SORTED=$(
    bcftools view \
        -H \
        "${VCF_SORTED}" \
    | wc -l
)


echo "Original records: ${N_RECORDS}"
echo "Sorted records:   ${N_SORTED}"


if [[ "${N_RECORDS}" -ne "${N_SORTED}" ]]
then
    echo "ERROR: VCF record count changed during sorting."
    exit 1
fi


# ============================================================
# INDEX SORTED VCF
# ============================================================

echo "Indexing sorted VCF..."


bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    "${VCF_SORTED}"


if [[ ! -s "${VCF_SORTED}.csi" && ! -s "${VCF_SORTED}.tbi" ]]
then
    echo "ERROR: sorted VCF index not created."
    exit 1
fi


N_INDEXED=$(
    bcftools index \
        -n \
        "${VCF_SORTED}"
)


echo "Indexed records:  ${N_INDEXED}"


if [[ "${N_INDEXED}" -ne "${N_SORTED}" ]]
then
    echo "ERROR: indexed record count differs from sorted VCF."
    exit 1
fi


# ============================================================
# ONE-VARIANT VEP SMOKE TEST
#
# IMPORTANT:
# Do NOT use:
#   bcftools view ... | head -1
#
# under set -o pipefail, because head can terminate bcftools
# with SIGPIPE.
# ============================================================

echo
echo "============================================================"
echo "VEP SMOKE TEST"
echo "============================================================"


ALL_VARIANTS_TMP="${OUTDIR}/all_variants.tmp"

FIRST_VARIANT_TMP="${OUTDIR}/first_variant.tmp"

FIRST_HEADER_TMP="${OUTDIR}/first_variant_header.tmp"


rm -f \
    "${TEST_VCF}" \
    "${TEST_VCF}.csi" \
    "${TEST_VCF}.tbi" \
    "${VEP_TEST_OUT}" \
    "${VEP_TEST_ERR}" \
    "${ALL_VARIANTS_TMP}" \
    "${FIRST_VARIANT_TMP}" \
    "${FIRST_HEADER_TMP}"


# ------------------------------------------------------------
# Write all records first.
# This avoids SIGPIPE caused by piping bcftools into head.
# ------------------------------------------------------------

bcftools view \
    -H \
    "${VCF_SORTED}" \
    > "${ALL_VARIANTS_TMP}"


if [[ ! -s "${ALL_VARIANTS_TMP}" ]]
then
    echo "ERROR: could not extract candidate records for smoke test."
    exit 1
fi


head -1 \
    "${ALL_VARIANTS_TMP}" \
    > "${FIRST_VARIANT_TMP}"


bcftools view \
    -h \
    "${VCF_SORTED}" \
    > "${FIRST_HEADER_TMP}"


if [[ ! -s "${FIRST_VARIANT_TMP}" ]]
then
    echo "ERROR: first test variant was not generated."
    exit 1
fi


if [[ ! -s "${FIRST_HEADER_TMP}" ]]
then
    echo "ERROR: smoke-test VCF header was not generated."
    exit 1
fi


cat \
    "${FIRST_HEADER_TMP}" \
    "${FIRST_VARIANT_TMP}" \
    | bgzip -c \
    > "${TEST_VCF}"


if [[ ! -s "${TEST_VCF}" ]]
then
    echo "ERROR: one-variant smoke-test VCF not produced."
    exit 1
fi


bcftools index \
    -f \
    "${TEST_VCF}"


N_TEST_VCF=$(
    bcftools view \
        -H \
        "${TEST_VCF}" \
    | wc -l
)


echo "Smoke-test VCF records: ${N_TEST_VCF}"


if [[ "${N_TEST_VCF}" -ne 1 ]]
then
    echo "ERROR: smoke-test VCF must contain exactly one SNP."
    exit 1
fi


rm -f \
    "${ALL_VARIANTS_TMP}" \
    "${FIRST_VARIANT_TMP}" \
    "${FIRST_HEADER_TMP}"


# ============================================================
# RUN VEP SMOKE TEST
#
# No --offline
# No --cache
#
# Annotation source:
#   custom GFF + custom FASTA
# ============================================================

set +e


vep \
    --input_file "${TEST_VCF}" \
    --output_file "${VEP_TEST_OUT}" \
    --gff "${GFF_VEP_GZ}" \
    --fasta "${FASTA}" \
    --tab \
    --symbol \
    --biotype \
    --numbers \
    --protein \
    --check_ref \
    --force_overwrite \
    --fork 1 \
    2> "${VEP_TEST_ERR}"


TEST_STATUS=$?


set -e


if [[ "${TEST_STATUS}" -ne 0 ]]
then

    echo
    echo "ERROR: VEP smoke test failed."
    echo
    echo "VEP smoke-test stderr:"
    echo "------------------------------------------------------------"

    cat "${VEP_TEST_ERR}"

    echo "------------------------------------------------------------"

    exit 1

fi


if [[ ! -s "${VEP_TEST_OUT}" ]]
then
    echo "ERROR: VEP smoke-test output missing or empty."
    exit 1
fi


N_TEST_ROWS=$(
    grep -vc '^#' \
        "${VEP_TEST_OUT}" \
    || true
)


echo "Smoke-test consequence rows: ${N_TEST_ROWS}"


if [[ "${N_TEST_ROWS}" -eq 0 ]]
then

    echo
    echo "ERROR: VEP smoke test produced zero annotation rows."

    echo
    echo "Smoke-test stderr:"
    cat "${VEP_TEST_ERR}"

    exit 1

fi


echo
echo "VEP SMOKE TEST PASSED"


# ============================================================
# RUN FULL VEP
#
# No --offline
# No --cache
# ============================================================

echo
echo "============================================================"
echo "RUNNING FULL ENSEMBL VEP v111"
echo "============================================================"


rm -f "${VEP_OUT}"


vep \
    --input_file "${VCF_SORTED}" \
    --output_file "${VEP_OUT}" \
    --gff "${GFF_VEP_GZ}" \
    --fasta "${FASTA}" \
    --tab \
    --symbol \
    --biotype \
    --numbers \
    --protein \
    --check_ref \
    --force_overwrite \
    --fork "${SLURM_CPUS_PER_TASK}"


# ============================================================
# VEP OUTPUT QC
# ============================================================

echo
echo "============================================================"
echo "VEP OUTPUT QC"
echo "============================================================"


if [[ ! -s "${VEP_OUT}" ]]
then
    echo "ERROR: full VEP output is missing or empty."
    exit 1
fi


HEADER_COUNT=$(
    grep -c '^#Uploaded_variation' \
        "${VEP_OUT}" \
    || true
)


echo "VEP header count: ${HEADER_COUNT}"


if [[ "${HEADER_COUNT}" -ne 1 ]]
then
    echo "ERROR: expected exactly one VEP header row."
    exit 1
fi


N_VEP_ROWS=$(
    grep -vc '^#' \
        "${VEP_OUT}" \
    || true
)


echo "Input VCF records:    ${N_SORTED}"
echo "VEP consequence rows: ${N_VEP_ROWS}"


if [[ "${N_VEP_ROWS}" -eq 0 ]]
then
    echo "ERROR: VEP produced zero consequence rows."
    exit 1
fi


# ============================================================
# UNIQUE VEP LOCATIONS
# ============================================================

N_VEP_LOCATIONS=$(
    grep -v '^#' \
        "${VEP_OUT}" \
    | cut -f2 \
    | sort -u \
    | wc -l
)


echo "Unique VEP locations: ${N_VEP_LOCATIONS}"


# ============================================================
# RAW CONSEQUENCE COUNTS
# ============================================================

echo
echo "============================================================"
echo "RAW VEP CONSEQUENCES"
echo "============================================================"


grep -v '^#' \
    "${VEP_OUT}" \
    | cut -f7 \
    | tr '&' '\n' \
    | sort \
    | uniq -c \
    | sort -nr \
    > "${RAW_COUNTS}"


cat "${RAW_COUNTS}"


# ============================================================
# BASIC VEP FORMAT QC
# ============================================================

echo
echo "============================================================"
echo "FIRST VEP RESULTS"
echo "============================================================"


grep -v '^##' \
    "${VEP_OUT}" \
    | head -15 \
    || true


# ============================================================
# WRITE FINAL QC SUMMARY
# ============================================================

{
    echo "Cactoblastis cactorum VEP candidate-SNP QC"
    echo "==========================================="
    echo

    echo "Date"
    echo "----"
    date
    echo

    echo "Reference FASTA"
    echo "---------------"
    echo "${FASTA}"
    echo

    echo "Original corrected GFF"
    echo "----------------------"
    echo "${GFF_ORIGINAL}"
    echo

    echo "VEP-specific FASTA-compatible GFF"
    echo "---------------------------------"
    echo "${GFF_VEP_GZ}"
    echo

    echo "Final all-SNP VCF"
    echo "-----------------"
    echo "${VCF_ALL}"
    echo

    echo "Candidate SNP master table"
    echo "--------------------------"
    echo "${MASTER}"
    echo

    echo "Candidate positions requested"
    echo "-----------------------------"
    echo "${N_REQUESTED}"
    echo

    echo "Candidate positions recovered"
    echo "-----------------------------"
    echo "${N_EXTRACTED_POS}"
    echo

    echo "Candidate positions missing"
    echo "---------------------------"
    echo "${N_MISSING}"
    echo

    echo "Candidate VCF records"
    echo "---------------------"
    echo "${N_RECORDS}"
    echo

    echo "Sorted candidate VCF records"
    echo "----------------------------"
    echo "${N_SORTED}"
    echo

    echo "Indexed candidate VCF records"
    echo "-----------------------------"
    echo "${N_INDEXED}"
    echo

    echo "Smoke-test VCF records"
    echo "----------------------"
    echo "${N_TEST_VCF}"
    echo

    echo "Smoke-test consequence rows"
    echo "---------------------------"
    echo "${N_TEST_ROWS}"
    echo

    echo "Full VEP consequence rows"
    echo "-------------------------"
    echo "${N_VEP_ROWS}"
    echo

    echo "Unique VEP locations"
    echo "--------------------"
    echo "${N_VEP_LOCATIONS}"
    echo

    echo "Reference contigs"
    echo "-----------------"
    echo "${N_FASTA_CONTIGS}"
    echo

    echo "GFF compatibility report"
    echo "------------------------"
    echo "${GFF_QC}"
    echo

    echo "VEP annotation mode"
    echo "-------------------"
    echo "Custom GFF + custom FASTA"
    echo "No Ensembl cache"
    echo "No --offline flag"
    echo

    echo "Software"
    echo "--------"
    echo "VEP: $(command -v vep)"
    echo "samtools: $(command -v samtools)"
    echo "bcftools: $(command -v bcftools)"
    echo "bgzip: $(command -v bgzip)"
    echo "tabix: $(command -v tabix)"

} > "${QC_SUMMARY}"


# ============================================================
# SUCCESS
# ============================================================

echo
echo "============================================================"
echo "VEP ANALYSIS COMPLETED SUCCESSFULLY"
echo "============================================================"

echo
echo "Final candidate VCF:"
echo "${VCF_SORTED}"

echo
echo "VEP output:"
echo "${VEP_OUT}"

echo
echo "Raw VEP consequence counts:"
echo "${RAW_COUNTS}"

echo
echo "GFF compatibility QC:"
echo "${GFF_QC}"

echo
echo "Final QC summary:"
echo "${QC_SUMMARY}"

echo
echo "Completed:"
date

echo
echo "DONE"

