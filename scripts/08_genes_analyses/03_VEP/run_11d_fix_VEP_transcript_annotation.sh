#!/bin/bash
#SBATCH --job-name=vep_ccac_fix
#SBATCH --output=logs/vep_ccac_fix_%j.out
#SBATCH --error=logs/vep_ccac_fix_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

set -euo pipefail

############################################################
## FINAL CORRECTED VEP
## Cactoblastis cactorum
##
## Purpose:
## Fix transcript recognition in custom AGAT GFF3.
##
## Previous QC established:
##   - 668 unique candidate positions
##   - 672 biallelic SNP records
##   - 668/668 positions recovered
##   - VCF/FASTA contigs compatible
##   - gene/exon/CDS/Parent structure present
##
## Problem:
##   transcript features lack biotype attributes and VEP
##   returned only intergenic_variant.
##
## Fix:
##   - transcript -> mRNA
##   - add biotype=protein_coding
##   - preserve ID/Parent hierarchy
##   - retain only FASTA-compatible contigs
##   - validate REF alleles
##   - smoke-test a SNP inside an annotated transcript
##   - require non-intergenic annotation before full run
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


for prog in vep samtools bcftools bgzip tabix python
do
    if ! command -v "${prog}" >/dev/null 2>&1
    then
        echo "ERROR: executable not found: ${prog}"
        exit 1
    fi
done


echo "============================================================"
echo "SOFTWARE"
echo "============================================================"

echo "VEP:      $(command -v vep)"
echo "samtools: $(command -v samtools)"
echo "bcftools: $(command -v bcftools)"
echo "bgzip:    $(command -v bgzip)"
echo "tabix:    $(command -v tabix)"
echo "python:   $(command -v python)"

echo
samtools --version | head -1
bcftools --version | head -1


# ============================================================
# PATHS
# ============================================================

ROOT="/home/dpoveda/micetral/Population_snps"

BASE="${ROOT}/baypass"

OUTDIR="${BASE}/selection_scan/VEP_ALL_selection_candidates"

mkdir -p "${OUTDIR}"


FASTA="${ROOT}/Ccactorum.fcs_gx_2026_FINAL.fa"

FAI="${FASTA}.fai"


GFF_ORIGINAL="${ROOT}/Ccac_functional-annotation/02_agat/Ccac/out_abril4/Ccac.annotations.corrected.gff3"


VCF="${OUTDIR}/ALL_selection_candidate_SNPs_for_VEP.sorted.vcf.gz"


# ============================================================
# NEW VEP GFF
# ============================================================

GFF_FIXED="${OUTDIR}/Ccac.VEP.fixed.biotype.gff3"

GFF_FIXED_GZ="${GFF_FIXED}.gz"


# ============================================================
# OUTPUTS
# ============================================================

VEP_OUT="${OUTDIR}/ALL_selection_candidate_SNPs.VEP_FINAL.tsv"

RAW_COUNTS="${OUTDIR}/VEP_FINAL_raw_consequence_counts.tsv"

QC="${OUTDIR}/VEP_FINAL_QC.txt"

GFF_QC="${OUTDIR}/VEP_FIXED_GFF_QC.txt"

REF_QC="${OUTDIR}/VCF_reference_allele_QC.txt"


# ============================================================
# SMOKE TEST FILES
# ============================================================

TEST_REGION="${OUTDIR}/VEP_test_selected_variant.tsv"

TEST_VCF="${OUTDIR}/VEP_test_annotated_variant.vcf.gz"

TEST_OUT="${OUTDIR}/VEP_test_annotated_variant.tsv"

TEST_ERR="${OUTDIR}/VEP_test_annotated_variant.err"


# ============================================================
# INPUT QC
# ============================================================

echo
echo "============================================================"
echo "INPUT QC"
echo "============================================================"


for f in \
    "${FASTA}" \
    "${FAI}" \
    "${GFF_ORIGINAL}" \
    "${VCF}"
do

    if [[ ! -s "${f}" ]]
    then
        echo "ERROR: missing/empty input:"
        echo "${f}"
        exit 1
    fi

done


if [[ ! -s "${VCF}.csi" && ! -s "${VCF}.tbi" ]]
then
    echo "Indexing candidate VCF..."

    bcftools index \
        --threads "${SLURM_CPUS_PER_TASK}" \
        -f \
        "${VCF}"
fi


N_VCF=$(bcftools view -H "${VCF}" | wc -l)

echo "Candidate VCF records: ${N_VCF}"


if [[ "${N_VCF}" -ne 672 ]]
then

    echo "WARNING:"
    echo "Expected 672 records from previous validated run."
    echo "Observed: ${N_VCF}"

fi


# ============================================================
# VALIDATE REFERENCE ALLELES
#
# Do this BEFORE VEP.
# If any REF allele disagrees with FASTA, abort.
# ============================================================

echo
echo "============================================================"
echo "REFERENCE ALLELE QC"
echo "============================================================"


rm -f "${REF_QC}"


set +e

bcftools norm \
    -f "${FASTA}" \
    -c e \
    -Ou \
    "${VCF}" \
    >/dev/null \
    2> "${REF_QC}"

REF_STATUS=$?

set -e


if [[ "${REF_STATUS}" -ne 0 ]]
then

    echo
    echo "ERROR: candidate VCF contains REF/FASTA mismatches."
    echo
    cat "${REF_QC}"

    exit 1

fi


echo "All VCF REF alleles match the reference FASTA."


# ============================================================
# BUILD VEP-COMPATIBLE GFF
# ============================================================

echo
echo "============================================================"
echo "BUILDING VEP-COMPATIBLE GFF"
echo "============================================================"


rm -f \
    "${GFF_FIXED}" \
    "${GFF_FIXED_GZ}" \
    "${GFF_FIXED_GZ}.tbi"


python - \
    "${GFF_ORIGINAL}" \
    "${FAI}" \
    "${GFF_FIXED}" \
    "${GFF_QC}" <<'PY'

import sys
from collections import Counter


gff_in = sys.argv[1]
fai = sys.argv[2]
gff_out = sys.argv[3]
qc_out = sys.argv[4]


# ============================================================
# Reference sequence order
# ============================================================

contig_order = {}

with open(fai) as fh:

    for i, line in enumerate(fh):

        chrom = line.rstrip("\n").split("\t")[0]

        contig_order[chrom] = i


# ============================================================
# Helper
# ============================================================

def add_attribute(attr, key, value):

    pieces = [
        x
        for x in attr.split(";")
        if x
    ]

    keys = {
        x.split("=", 1)[0]
        for x in pieces
        if "=" in x
    }

    if key not in keys:

        pieces.append(
            f"{key}={value}"
        )

    return ";".join(pieces)


# ============================================================
# Read + modify
# ============================================================

records = []

types_before = Counter()

types_after = Counter()

removed_contigs = Counter()

n_biotype_added = 0

n_transcript_converted = 0

idx = 0

in_fasta = False


with open(gff_in) as fh:

    for raw in fh:


        if raw.startswith("##FASTA"):

            in_fasta = True
            continue


        if in_fasta:

            continue


        if raw.startswith("#"):

            continue


        line = raw.rstrip("\n")


        if not line:

            continue


        f = line.split("\t")


        if len(f) != 9:

            continue


        chrom = f[0]


        if chrom not in contig_order:

            removed_contigs[chrom] += 1

            continue


        old_type = f[2]

        types_before[old_type] += 1


        # ----------------------------------------------------
        # Convert generic transcript -> standard GFF3 mRNA
        # ----------------------------------------------------

        if f[2] == "transcript":

            f[2] = "mRNA"

            n_transcript_converted += 1


        # ----------------------------------------------------
        # Add protein_coding biotype to genes and transcripts
        #
        # This annotation is AUGUSTUS protein-coding gene
        # prediction and all retained transcript models have
        # CDS features.
        # ----------------------------------------------------

        if f[2] in {
            "gene",
            "mRNA"
        }:

            old_attr = f[8]

            f[8] = add_attribute(
                f[8],
                "biotype",
                "protein_coding"
            )


            if f[8] != old_attr:

                n_biotype_added += 1


        types_after[f[2]] += 1


        start = int(f[3])

        end = int(f[4])


        records.append(
            (
                contig_order[chrom],
                start,
                end,
                idx,
                "\t".join(f)
            )
        )


        idx += 1


# ============================================================
# Sort
# ============================================================

records.sort(
    key=lambda x: (
        x[0],
        x[1],
        x[2],
        x[3]
    )
)


# ============================================================
# Write
# ============================================================

with open(gff_out, "w") as out:

    out.write(
        "##gff-version 3\n"
    )

    for _, _, _, _, line in records:

        out.write(
            line + "\n"
        )


# ============================================================
# Validate hierarchy
# ============================================================

genes = set()

transcripts = set()

transcript_parents = {}

child_parents = []

n_gene = 0
n_mrna = 0
n_exon = 0
n_cds = 0


def parse_attrs(s):

    result = {}

    for item in s.split(";"):

        if "=" not in item:
            continue

        k, v = item.split(
            "=",
            1
        )

        result[k] = v

    return result


with open(gff_out) as fh:

    for line in fh:

        if line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")

        attrs = parse_attrs(
            f[8]
        )


        if f[2] == "gene":

            n_gene += 1

            if "ID" in attrs:

                genes.add(
                    attrs["ID"]
                )


        elif f[2] == "mRNA":

            n_mrna += 1

            tid = attrs.get("ID")

            parent = attrs.get("Parent")


            if tid:

                transcripts.add(
                    tid
                )


            if tid and parent:

                transcript_parents[
                    tid
                ] = parent


        elif f[2] == "exon":

            n_exon += 1

            parent = attrs.get(
                "Parent"
            )

            if parent:

                child_parents.extend(
                    parent.split(",")
                )


        elif f[2] == "CDS":

            n_cds += 1

            parent = attrs.get(
                "Parent"
            )

            if parent:

                child_parents.extend(
                    parent.split(",")
                )


bad_transcript_parent = [
    (t, p)
    for t, p in transcript_parents.items()
    if p not in genes
]


bad_child_parent = [
    p
    for p in child_parents
    if p not in transcripts
]


# ============================================================
# Report
# ============================================================

with open(qc_out, "w") as out:

    out.write(
        f"FASTA_contigs\t{len(contig_order)}\n"
    )

    out.write(
        f"Retained_GFF_records\t{len(records)}\n"
    )

    out.write(
        f"Removed_GFF_contigs\t{len(removed_contigs)}\n"
    )

    out.write(
        f"Transcript_to_mRNA_conversions\t{n_transcript_converted}\n"
    )

    out.write(
        f"Biotype_attributes_added\t{n_biotype_added}\n"
    )

    out.write(
        f"Genes\t{n_gene}\n"
    )

    out.write(
        f"mRNA\t{n_mrna}\n"
    )

    out.write(
        f"Exons\t{n_exon}\n"
    )

    out.write(
        f"CDS\t{n_cds}\n"
    )

    out.write(
        f"Invalid_mRNA_gene_parent_links\t{len(bad_transcript_parent)}\n"
    )

    out.write(
        f"Invalid_exon_CDS_mRNA_parent_links\t{len(bad_child_parent)}\n"
    )

    out.write(
        "\nFeature_types_after_conversion\n"
    )


    for feature, n in types_after.most_common():

        out.write(
            f"{feature}\t{n}\n"
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


print(
    "Retained records:",
    len(records)
)

print(
    "transcript -> mRNA:",
    n_transcript_converted
)

print(
    "biotype attributes added:",
    n_biotype_added
)

print(
    "Genes:",
    n_gene
)

print(
    "mRNA:",
    n_mrna
)

print(
    "Exons:",
    n_exon
)

print(
    "CDS:",
    n_cds
)

print(
    "Invalid mRNA->gene links:",
    len(bad_transcript_parent)
)

print(
    "Invalid exon/CDS->mRNA links:",
    len(bad_child_parent)
)


if n_mrna == 0:

    raise RuntimeError(
        "No mRNA features created"
    )


if n_cds == 0:

    raise RuntimeError(
        "No CDS features retained"
    )


if bad_transcript_parent:

    raise RuntimeError(
        "Invalid mRNA -> gene parent relationships"
    )


if bad_child_parent:

    raise RuntimeError(
        "Invalid exon/CDS -> mRNA parent relationships"
    )
PY


if [[ ! -s "${GFF_FIXED}" ]]
then

    echo "ERROR: corrected VEP GFF was not produced."
    exit 1

fi


echo
cat "${GFF_QC}"


# ============================================================
# CRITICAL GFF CHECKS
# ============================================================

N_MRNA=$(
    grep -c $'\tmRNA\t' \
        "${GFF_FIXED}" \
    || true
)


N_BIOTYPE=$(
    grep -c 'biotype=protein_coding' \
        "${GFF_FIXED}" \
    || true
)


echo
echo "mRNA records:                  ${N_MRNA}"
echo "protein_coding biotype attrs: ${N_BIOTYPE}"


if [[ "${N_MRNA}" -eq 0 ]]
then

    echo "ERROR: no mRNA records."
    exit 1

fi


if [[ "${N_BIOTYPE}" -eq 0 ]]
then

    echo "ERROR: no protein_coding biotypes."
    exit 1

fi


# ============================================================
# BGZIP + TABIX
# ============================================================

echo
echo "============================================================"
echo "COMPRESSING / INDEXING FIXED GFF"
echo "============================================================"


bgzip \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -c "${GFF_FIXED}" \
    > "${GFF_FIXED_GZ}"


tabix \
    -f \
    -p gff \
    "${GFF_FIXED_GZ}"


if [[ ! -s "${GFF_FIXED_GZ}.tbi" ]]
then

    echo "ERROR: corrected GFF index not created."
    exit 1

fi


# ============================================================
# SELECT A REAL ANNOTATED SNP FOR SMOKE TEST
#
# Prefer:
#   1 CDS overlap
#   2 exon overlap
#   3 mRNA overlap
#
# This prevents testing an arbitrary intergenic SNP.
# ============================================================

echo
echo "============================================================"
echo "SELECTING ANNOTATED SMOKE-TEST SNP"
echo "============================================================"


python - \
    "${VCF}" \
    "${GFF_FIXED}" \
    "${TEST_REGION}" <<'PY'

import gzip
import sys


vcf_file = sys.argv[1]
gff_file = sys.argv[2]
outfile = sys.argv[3]


variants = []


with gzip.open(vcf_file, "rt") as fh:

    for line in fh:

        if line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")

        variants.append(
            (
                f[0],
                int(f[1]),
                f[3],
                f[4]
            )
        )


variant_by_contig = {}


for variant in variants:

    variant_by_contig.setdefault(
        variant[0],
        []
    ).append(
        variant
    )


priority = {
    "CDS": 1,
    "exon": 2,
    "mRNA": 3
}


best = None


with open(gff_file) as fh:

    for line in fh:

        if line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")


        feature = f[2]


        if feature not in priority:

            continue


        chrom = f[0]


        if chrom not in variant_by_contig:

            continue


        start = int(f[3])
        end = int(f[4])


        for variant in variant_by_contig[chrom]:

            _, pos, ref, alt = variant


            if start <= pos <= end:

                candidate = (
                    priority[feature],
                    chrom,
                    pos,
                    ref,
                    alt,
                    feature,
                    start,
                    end
                )


                if (
                    best is None
                    or candidate[0] < best[0]
                ):

                    best = candidate


                    if candidate[0] == 1:
                        break


        if best is not None and best[0] == 1:
            break


if best is None:

    raise RuntimeError(
        "No candidate SNP overlaps CDS/exon/mRNA annotation"
    )


_, chrom, pos, ref, alt, feature, start, end = best


with open(outfile, "w") as out:

    out.write(
        "CHROM\tPOS\tREF\tALT\tFEATURE\tFEATURE_START\tFEATURE_END\n"
    )

    out.write(
        f"{chrom}\t{pos}\t{ref}\t{alt}\t"
        f"{feature}\t{start}\t{end}\n"
    )


print(
    "Selected smoke-test SNP:",
    chrom,
    pos,
    ref,
    alt,
    "inside",
    feature
)
PY


cat "${TEST_REGION}"


TEST_CHROM=$(
    awk 'NR==2 {print $1}' \
        "${TEST_REGION}"
)


TEST_POS=$(
    awk 'NR==2 {print $2}' \
        "${TEST_REGION}"
)


# ============================================================
# BUILD TEST VCF
# ============================================================

rm -f \
    "${TEST_VCF}" \
    "${TEST_VCF}.csi" \
    "${TEST_VCF}.tbi" \
    "${TEST_OUT}" \
    "${TEST_ERR}"


bcftools view \
    -r "${TEST_CHROM}:${TEST_POS}-${TEST_POS}" \
    -Oz \
    -o "${TEST_VCF}" \
    "${VCF}"


bcftools index \
    -f \
    "${TEST_VCF}"


N_TEST_VARIANTS=$(
    bcftools view \
        -H \
        "${TEST_VCF}" \
    | wc -l
)


echo "Smoke-test VCF records: ${N_TEST_VARIANTS}"


if [[ "${N_TEST_VARIANTS}" -eq 0 ]]
then

    echo "ERROR: selected smoke-test SNP not recovered."
    exit 1

fi


# ============================================================
# VEP SMOKE TEST
# ============================================================

echo
echo "============================================================"
echo "VEP ANNOTATED-SNP SMOKE TEST"
echo "============================================================"


set +e


vep \
    --input_file "${TEST_VCF}" \
    --output_file "${TEST_OUT}" \
    --gff "${GFF_FIXED_GZ}" \
    --fasta "${FASTA}" \
    --tab \
    --symbol \
    --biotype \
    --numbers \
    --protein \
    --force_overwrite \
    --fork 1 \
    2> "${TEST_ERR}"


TEST_STATUS=$?


set -e


if [[ "${TEST_STATUS}" -ne 0 ]]
then

    echo
    echo "ERROR: VEP smoke test crashed."

    cat "${TEST_ERR}"

    exit 1

fi


if [[ ! -s "${TEST_OUT}" ]]
then

    echo "ERROR: smoke-test output is empty."
    exit 1

fi


echo
echo "Smoke-test VEP output:"

grep -v '^##' \
    "${TEST_OUT}" \
    || true


N_TEST_ROWS=$(
    grep -vc '^#' \
        "${TEST_OUT}" \
    || true
)


N_TEST_NON_INTERGENIC=$(
    grep -v '^#' \
        "${TEST_OUT}" \
    | awk -F'\t' '$7 != "intergenic_variant" {n++} END {print n+0}'
)


echo
echo "Smoke-test consequence rows:      ${N_TEST_ROWS}"
echo "Smoke-test non-intergenic rows:  ${N_TEST_NON_INTERGENIC}"


if [[ "${N_TEST_ROWS}" -eq 0 ]]
then

    echo "ERROR: VEP smoke test produced zero rows."

    exit 1

fi


if [[ "${N_TEST_NON_INTERGENIC}" -eq 0 ]]
then

    echo
    echo "ERROR:"
    echo "The selected SNP is physically inside an annotated"
    echo "CDS/exon/transcript, but VEP still classified it only"
    echo "as intergenic."
    echo
    echo "Do NOT use these annotations."

    exit 1

fi


echo
echo "SMOKE TEST PASSED:"
echo "VEP is recognizing transcript annotation."


# ============================================================
# RUN FULL VEP
# ============================================================

echo
echo "============================================================"
echo "RUNNING FULL VEP"
echo "============================================================"


rm -f "${VEP_OUT}"


vep \
    --input_file "${VCF}" \
    --output_file "${VEP_OUT}" \
    --gff "${GFF_FIXED_GZ}" \
    --fasta "${FASTA}" \
    --tab \
    --symbol \
    --biotype \
    --numbers \
    --protein \
    --force_overwrite \
    --fork "${SLURM_CPUS_PER_TASK}"


# ============================================================
# FULL OUTPUT QC
# ============================================================

echo
echo "============================================================"
echo "FULL VEP OUTPUT QC"
echo "============================================================"


if [[ ! -s "${VEP_OUT}" ]]
then

    echo "ERROR: full VEP output missing."
    exit 1

fi


HEADER_COUNT=$(
    grep -c '^#Uploaded_variation' \
        "${VEP_OUT}" \
    || true
)


if [[ "${HEADER_COUNT}" -ne 1 ]]
then

    echo "ERROR: invalid VEP header."
    exit 1

fi


N_ROWS=$(
    grep -vc '^#' \
        "${VEP_OUT}" \
    || true
)


N_INTERGENIC=$(
    grep -v '^#' \
        "${VEP_OUT}" \
    | awk -F'\t' '$7=="intergenic_variant" {n++} END {print n+0}'
)


N_NON_INTERGENIC=$(
    grep -v '^#' \
        "${VEP_OUT}" \
    | awk -F'\t' '$7!="intergenic_variant" {n++} END {print n+0}'
)


echo "Candidate VCF records:    ${N_VCF}"
echo "VEP consequence rows:     ${N_ROWS}"
echo "Intergenic rows:          ${N_INTERGENIC}"
echo "Non-intergenic rows:      ${N_NON_INTERGENIC}"


if [[ "${N_ROWS}" -eq 0 ]]
then

    echo "ERROR: VEP produced zero rows."
    exit 1

fi


if [[ "${N_NON_INTERGENIC}" -eq 0 ]]
then

    echo
    echo "ERROR:"
    echo "Full VEP output is still entirely intergenic."
    echo "Annotation must not be used."

    exit 1

fi


# ============================================================
# CONSEQUENCE COUNTS
#
# Split compound consequences joined by &
# ============================================================

grep -v '^#' \
    "${VEP_OUT}" \
    | cut -f7 \
    | tr '&' '\n' \
    | sort \
    | uniq -c \
    | sort -nr \
    > "${RAW_COUNTS}"


echo
echo "============================================================"
echo "RAW CONSEQUENCE COUNTS"
echo "============================================================"

cat "${RAW_COUNTS}"


# ============================================================
# UNIQUE LOCATIONS
# ============================================================

N_LOCATIONS=$(
    grep -v '^#' \
        "${VEP_OUT}" \
    | cut -f2 \
    | sort -u \
    | wc -l
)


# ============================================================
# IMPACT COUNTS
# ============================================================

echo
echo "============================================================"
echo "VEP IMPACT COUNTS"
echo "============================================================"


grep -v '^#' \
    "${VEP_OUT}" \
    | cut -f14 \
    | sort \
    | uniq -c \
    | sort -nr


# ============================================================
# FINAL QC REPORT
# ============================================================

{
    echo "FINAL VEP QC - Cactoblastis cactorum"
    echo "===================================="
    echo

    echo "Date"
    echo "----"
    date
    echo

    echo "Candidate VCF"
    echo "-------------"
    echo "${VCF}"
    echo

    echo "Candidate VCF records"
    echo "---------------------"
    echo "${N_VCF}"
    echo

    echo "Reference allele validation"
    echo "---------------------------"
    echo "PASS"
    echo

    echo "Original annotation"
    echo "-------------------"
    echo "${GFF_ORIGINAL}"
    echo

    echo "VEP-fixed GFF"
    echo "-------------"
    echo "${GFF_FIXED_GZ}"
    echo

    echo "mRNA features"
    echo "-------------"
    echo "${N_MRNA}"
    echo

    echo "protein_coding biotype attributes"
    echo "---------------------------------"
    echo "${N_BIOTYPE}"
    echo

    echo "Annotated SNP smoke test"
    echo "------------------------"
    cat "${TEST_REGION}"
    echo

    echo "Smoke-test VEP rows"
    echo "-------------------"
    echo "${N_TEST_ROWS}"
    echo

    echo "Smoke-test non-intergenic rows"
    echo "------------------------------"
    echo "${N_TEST_NON_INTERGENIC}"
    echo

    echo "Full VEP rows"
    echo "-------------"
    echo "${N_ROWS}"
    echo

    echo "Unique VEP locations"
    echo "--------------------"
    echo "${N_LOCATIONS}"
    echo

    echo "Intergenic consequence rows"
    echo "---------------------------"
    echo "${N_INTERGENIC}"
    echo

    echo "Non-intergenic consequence rows"
    echo "-------------------------------"
    echo "${N_NON_INTERGENIC}"
    echo

    echo "Raw consequences"
    echo "----------------"
    cat "${RAW_COUNTS}"

} > "${QC}"


# ============================================================
# FIRST RESULTS
# ============================================================

echo
echo "============================================================"
echo "FIRST FINAL VEP RESULTS"
echo "============================================================"


grep -v '^##' \
    "${VEP_OUT}" \
    | head -20 \
    || true


# ============================================================
# SUCCESS
# ============================================================

echo
echo "============================================================"
echo "FINAL VEP ANALYSIS PASSED"
echo "============================================================"

echo
echo "VEP output:"
echo "${VEP_OUT}"

echo
echo "Consequence counts:"
echo "${RAW_COUNTS}"

echo
echo "QC:"
echo "${QC}"

echo
echo "Fixed VEP GFF:"
echo "${GFF_FIXED_GZ}"

echo
echo "DONE"
echo "$(date)"

