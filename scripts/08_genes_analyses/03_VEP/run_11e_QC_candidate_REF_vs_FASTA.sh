#!/bin/bash
#SBATCH --job-name=qc_ref_ccac
#SBATCH --output=logs/qc_ref_ccac_%j.out
#SBATCH --error=logs/qc_ref_ccac_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

module purge
module load devel/Miniconda/Miniconda3
module load bioinfo/samtools/1.23
module load bioinfo/Bcftools/1.21

source activate bio-env

export LC_ALL=C
export LANG=C


# ============================================================
# PATHS
# ============================================================

ROOT="/home/dpoveda/micetral/Population_snps"

BASE="${ROOT}/baypass"

OUTDIR="${BASE}/selection_scan/VEP_ALL_selection_candidates"

FASTA="${ROOT}/Ccactorum.fcs_gx_2026_FINAL.fa"

VCF="${OUTDIR}/ALL_selection_candidate_SNPs_for_VEP.sorted.vcf.gz"

VARIANTS="${OUTDIR}/QC_candidate_variants_REF_ALT.tsv"

REGIONS="${OUTDIR}/QC_candidate_regions_for_faidx.txt"

FASTA_BASES_RAW="${OUTDIR}/QC_candidate_FASTA_bases.fa"

FINAL="${OUTDIR}/QC_candidate_REF_FASTA_classification.tsv"

SUMMARY="${OUTDIR}/QC_candidate_REF_FASTA_summary.tsv"

MISMATCHES="${OUTDIR}/QC_candidate_REF_FASTA_mismatches.tsv"


# ============================================================
# INPUT CHECKS
# ============================================================

mkdir -p logs

[[ -s "${FASTA}" ]] || {
    echo "ERROR: FASTA not found"
    exit 1
}

[[ -s "${VCF}" ]] || {
    echo "ERROR: candidate VCF not found"
    exit 1
}


if [[ ! -s "${FASTA}.fai" ]]
then
    samtools faidx "${FASTA}"
fi


echo "============================================================"
echo "CANDIDATE VCF REF vs FASTA QC"
echo "============================================================"

echo "FASTA:"
echo "${FASTA}"

echo
echo "VCF:"
echo "${VCF}"


# ============================================================
# EXTRACT VCF ALLELES
# ============================================================

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "${VCF}" \
    > "${VARIANTS}"


N_VARIANTS=$(wc -l < "${VARIANTS}")

echo
echo "VCF records: ${N_VARIANTS}"


# ============================================================
# BUILD FASTA REGIONS
# ============================================================

awk 'BEGIN{OFS=""} {print $1,":",$2,"-",$2}' \
    "${VARIANTS}" \
    > "${REGIONS}"


# ============================================================
# EXTRACT REFERENCE BASES
#
# samtools faidx -r safely handles all requested positions.
# ============================================================

samtools faidx \
    -r "${REGIONS}" \
    "${FASTA}" \
    > "${FASTA_BASES_RAW}"


# ============================================================
# CLASSIFY EVERY SNP
# ============================================================

python - \
    "${VARIANTS}" \
    "${FASTA_BASES_RAW}" \
    "${FINAL}" \
    "${SUMMARY}" \
    "${MISMATCHES}" <<'PY'

import sys
from collections import Counter


variants_file = sys.argv[1]
fasta_file = sys.argv[2]
final_file = sys.argv[3]
summary_file = sys.argv[4]
mismatch_file = sys.argv[5]


# ============================================================
# Complement helper
# ============================================================

COMP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N"
}


def complement(base):

    return "".join(
        COMP.get(
            x.upper(),
            "?"
        )
        for x in base
    )


# ============================================================
# Read VCF variants
# ============================================================

variants = []


with open(variants_file) as fh:

    for line in fh:

        if not line.strip():
            continue

        chrom, pos, ref, alt = line.rstrip("\n").split("\t")[:4]

        variants.append(
            {
                "CHROM": chrom,
                "POS": int(pos),
                "REF": ref.upper(),
                "ALT": alt.upper()
            }
        )


# ============================================================
# Read faidx output
# ============================================================

fasta_bases = []

current = []


with open(fasta_file) as fh:

    for line in fh:

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if current:

                fasta_bases.append(
                    "".join(current).upper()
                )

                current = []

        else:

            current.append(line)


if current:

    fasta_bases.append(
        "".join(current).upper()
    )


if len(variants) != len(fasta_bases):

    raise RuntimeError(
        f"VCF variants ({len(variants)}) != "
        f"FASTA bases ({len(fasta_bases)})"
    )


# ============================================================
# Classification
# ============================================================

results = []

counts = Counter()


for var, fasta_seq in zip(
    variants,
    fasta_bases
):

    ref = var["REF"]
    alt = var["ALT"]

    fasta = fasta_seq.upper()


    # Candidate analysis is SNP-only, so each should be length 1.
    if (
        len(ref) != 1
        or len(alt) != 1
        or len(fasta) != 1
    ):

        category = "NON_SINGLE_BASE"

    elif ref == fasta:

        category = "EXACT_REF_MATCH"

    elif alt == fasta:

        category = "REF_ALT_SWAP"

    elif complement(ref) == fasta:

        category = "STRAND_COMPLEMENT"

    elif complement(alt) == fasta:

        category = "STRAND_COMPLEMENT_PLUS_SWAP"

    else:

        category = "UNRESOLVED_MISMATCH"


    counts[category] += 1


    results.append(
        (
            var["CHROM"],
            var["POS"],
            ref,
            alt,
            fasta,
            complement(ref),
            complement(alt),
            category
        )
    )


# ============================================================
# Write complete table
# ============================================================

header = [
    "CHROM",
    "POS",
    "VCF_REF",
    "VCF_ALT",
    "FASTA_REF",
    "Complement_REF",
    "Complement_ALT",
    "Classification"
]


with open(final_file, "w") as out:

    out.write(
        "\t".join(header) + "\n"
    )

    for row in results:

        out.write(
            "\t".join(
                map(
                    str,
                    row
                )
            )
            + "\n"
        )


# ============================================================
# Write mismatch-only table
# ============================================================

with open(mismatch_file, "w") as out:

    out.write(
        "\t".join(header) + "\n"
    )

    for row in results:

        if row[-1] != "EXACT_REF_MATCH":

            out.write(
                "\t".join(
                    map(
                        str,
                        row
                    )
                )
                + "\n"
            )


# ============================================================
# Summary
# ============================================================

order = [
    "EXACT_REF_MATCH",
    "REF_ALT_SWAP",
    "STRAND_COMPLEMENT",
    "STRAND_COMPLEMENT_PLUS_SWAP",
    "UNRESOLVED_MISMATCH",
    "NON_SINGLE_BASE"
]


with open(summary_file, "w") as out:

    out.write(
        "Classification\tN\tPercent\n"
    )

    total = len(results)


    for category in order:

        n = counts.get(
            category,
            0
        )

        pct = (
            100.0 * n / total
            if total
            else 0
        )

        out.write(
            f"{category}\t{n}\t{pct:.4f}\n"
        )


print()
print("============================================")
print("REF / FASTA CLASSIFICATION")
print("============================================")

print(
    "Total VCF records:",
    len(results)
)

print()


for category in order:

    n = counts.get(
        category,
        0
    )

    pct = (
        100.0 * n / len(results)
        if results
        else 0
    )

    print(
        f"{category:32s} "
        f"{n:5d} "
        f"({pct:7.3f}%)"
    )


print()
print(
    "Full table:",
    final_file
)

print(
    "Mismatch table:",
    mismatch_file
)

print(
    "Summary:",
    summary_file
)
PY


# ============================================================
# SHOW MISMATCH EXAMPLES
# ============================================================

echo
echo "============================================================"
echo "FIRST 30 NON-EXACT RECORDS"
echo "============================================================"

head -31 "${MISMATCHES}"


echo
echo "============================================================"
echo "QC COMPLETE"
echo "============================================================"

cat "${SUMMARY}"

