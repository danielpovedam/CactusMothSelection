#!/bin/bash
#SBATCH --job-name=vep_ccac_final2
#SBATCH --output=logs/vep_ccac_final2_%j.out
#SBATCH --error=logs/vep_ccac_final2_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

set -euo pipefail

############################################################
## FINAL NORMALIZATION + VEP
## Cactoblastis cactorum
##
## Validated candidate set:
##   672 VCF records
##   668 unique positions
##
## REF/FASTA QC:
##   EXACT_REF_MATCH              501
##   REF_ALT_SWAP                 170
##   STRAND_COMPLEMENT              1
##   STRAND_COMPLEMENT_PLUS_SWAP    0
##   UNRESOLVED_MISMATCH            0
##
## This script:
##   1. Creates a VEP-only sites VCF
##   2. Normalizes REF/ALT to FASTA
##   3. Validates all REF alleles
##   4. Builds VEP-compatible GFF
##   5. Validates GFF hierarchy
##   6. Smoke-tests a SNP inside CDS/exon/mRNA
##   7. Runs full VEP
##   8. Rejects an all-intergenic output
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


VCF_ORIGINAL="${OUTDIR}/ALL_selection_candidate_SNPs_for_VEP.sorted.vcf.gz"


REF_CLASSIFICATION="${OUTDIR}/QC_candidate_REF_FASTA_classification.tsv"


# ============================================================
# NORMALIZED VEP VCF
# ============================================================

VCF_NORMALIZED="${OUTDIR}/ALL_selection_candidate_SNPs.VEP_normalized.vcf.gz"

NORMALIZATION_TABLE="${OUTDIR}/VEP_variant_normalization_table.tsv"

NORMALIZATION_QC="${OUTDIR}/VEP_variant_normalization_QC.tsv"


# ============================================================
# VEP GFF
# ============================================================

GFF_FIXED="${OUTDIR}/Ccac.VEP.fixed.biotype.gff3"

GFF_FIXED_GZ="${GFF_FIXED}.gz"

GFF_QC="${OUTDIR}/VEP_FIXED_GFF_QC.txt"


# ============================================================
# VEP OUTPUTS
# ============================================================

VEP_OUT="${OUTDIR}/ALL_selection_candidate_SNPs.VEP_FINAL.tsv"

RAW_COUNTS="${OUTDIR}/VEP_FINAL_raw_consequence_counts.tsv"

FINAL_QC="${OUTDIR}/VEP_FINAL_QC.txt"


# ============================================================
# SMOKE TEST
# ============================================================

TEST_INFO="${OUTDIR}/VEP_test_selected_variant.tsv"

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
    "${VCF_ORIGINAL}" \
    "${REF_CLASSIFICATION}"
do

    if [[ ! -s "${f}" ]]
    then
        echo "ERROR: missing/empty input:"
        echo "${f}"
        exit 1
    fi

done


N_ORIGINAL=$(
    bcftools view \
        -H \
        "${VCF_ORIGINAL}" \
    | wc -l
)


echo "Original candidate VCF records: ${N_ORIGINAL}"


if [[ "${N_ORIGINAL}" -ne 672 ]]
then
    echo "WARNING: expected 672 candidate records."
fi


# ============================================================
# NORMALIZE VARIANT ALLELES FOR VEP
#
# IMPORTANT:
# This creates a SITES-ONLY VCF.
#
# We do NOT alter the population-genomic VCF because swapping
# REF/ALT in a genotype VCF would require genotype recoding.
# ============================================================

echo
echo "============================================================"
echo "NORMALIZING REF / ALT FOR VEP"
echo "============================================================"


rm -f \
    "${VCF_NORMALIZED}" \
    "${VCF_NORMALIZED}.csi" \
    "${VCF_NORMALIZED}.tbi" \
    "${NORMALIZATION_TABLE}" \
    "${NORMALIZATION_QC}"


python - \
    "${VCF_ORIGINAL}" \
    "${REF_CLASSIFICATION}" \
    "${FASTA}" \
    "${OUTDIR}/VEP_normalized.unsorted.vcf" \
    "${NORMALIZATION_TABLE}" \
    "${NORMALIZATION_QC}" <<'PY'

import gzip
import sys
from collections import Counter


vcf_file = sys.argv[1]
classification_file = sys.argv[2]
fasta_file = sys.argv[3]
out_vcf = sys.argv[4]
out_table = sys.argv[5]
out_qc = sys.argv[6]


COMP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}


def comp(base):

    return "".join(
        COMP[x]
        for x in base.upper()
    )


# ============================================================
# Read classifications
#
# Multiple records may exist at the same CHROM/POS,
# therefore key using:
#   CHROM + POS + VCF_REF + VCF_ALT
# ============================================================

class_map = {}


with open(classification_file) as fh:

    header = fh.readline().rstrip("\n").split("\t")

    index = {
        name: i
        for i, name in enumerate(header)
    }


    required = [
        "CHROM",
        "POS",
        "VCF_REF",
        "VCF_ALT",
        "FASTA_REF",
        "Classification"
    ]


    missing = [
        x
        for x in required
        if x not in index
    ]


    if missing:

        raise RuntimeError(
            "Missing classification columns: "
            + ",".join(missing)
        )


    for line in fh:

        if not line.strip():
            continue

        f = line.rstrip("\n").split("\t")


        key = (
            f[index["CHROM"]],
            int(f[index["POS"]]),
            f[index["VCF_REF"]].upper(),
            f[index["VCF_ALT"]].upper()
        )


        class_map[key] = {
            "FASTA_REF":
                f[index["FASTA_REF"]].upper(),

            "Classification":
                f[index["Classification"]]
        }


# ============================================================
# Read original VCF
# ============================================================

records = []

counts = Counter()


with gzip.open(vcf_file, "rt") as fh:

    for line in fh:


        if line.startswith("#"):

            continue


        f = line.rstrip("\n").split("\t")


        chrom = f[0]
        pos = int(f[1])
        vid = f[2]
        ref = f[3].upper()
        alt = f[4].upper()


        if "," in alt:

            raise RuntimeError(
                f"Multiallelic record unexpectedly found: "
                f"{chrom}:{pos}"
            )


        key = (
            chrom,
            pos,
            ref,
            alt
        )


        if key not in class_map:

            raise RuntimeError(
                "Variant absent from classification table: "
                + ":".join(map(str, key))
            )


        info = class_map[key]

        fasta_ref = info["FASTA_REF"]

        category = info["Classification"]


        # ----------------------------------------------------
        # Normalize alleles
        # ----------------------------------------------------

        if category == "EXACT_REF_MATCH":

            new_ref = ref
            new_alt = alt


        elif category == "REF_ALT_SWAP":

            new_ref = alt
            new_alt = ref


        elif category == "STRAND_COMPLEMENT":

            new_ref = comp(ref)
            new_alt = comp(alt)


        elif category == "STRAND_COMPLEMENT_PLUS_SWAP":

            new_ref = comp(alt)
            new_alt = comp(ref)


        else:

            raise RuntimeError(
                f"Unsafe classification {category} "
                f"at {chrom}:{pos}"
            )


        # ----------------------------------------------------
        # Critical validation
        # ----------------------------------------------------

        if new_ref != fasta_ref:

            raise RuntimeError(
                f"Normalized REF still disagrees with FASTA at "
                f"{chrom}:{pos}: "
                f"new REF={new_ref}, FASTA={fasta_ref}"
            )


        if new_ref == new_alt:

            raise RuntimeError(
                f"REF == ALT after normalization at "
                f"{chrom}:{pos}"
            )


        counts[category] += 1


        records.append(
            (
                chrom,
                pos,
                vid,
                ref,
                alt,
                new_ref,
                new_alt,
                fasta_ref,
                category
            )
        )


# ============================================================
# Write minimal sites-only VCF
# ============================================================

with open(out_vcf, "w") as out:

    out.write(
        "##fileformat=VCFv4.2\n"
    )

    out.write(
        "##source=Ccactorum_VEP_normalization\n"
    )

    out.write(
        "##INFO=<ID=ORIG_REF,Number=1,Type=String,"
        "Description=\"Original REF allele\">\n"
    )

    out.write(
        "##INFO=<ID=ORIG_ALT,Number=1,Type=String,"
        "Description=\"Original ALT allele\">\n"
    )

    out.write(
        "##INFO=<ID=NORM_CLASS,Number=1,Type=String,"
        "Description=\"Normalization applied for VEP\">\n"
    )


    # --------------------------------------------------------
    # Copy contig headers from reference FASTA
    # --------------------------------------------------------

    fai = fasta_file + ".fai"


    with open(fai) as fh:

        for line in fh:

            f = line.rstrip("\n").split("\t")

            out.write(
                f"##contig=<ID={f[0]},length={f[1]}>\n"
            )


    out.write(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )


    for (
        chrom,
        pos,
        vid,
        old_ref,
        old_alt,
        new_ref,
        new_alt,
        fasta_ref,
        category
    ) in records:

        info = (
            f"ORIG_REF={old_ref};"
            f"ORIG_ALT={old_alt};"
            f"NORM_CLASS={category}"
        )


        out.write(
            f"{chrom}\t{pos}\t{vid}\t"
            f"{new_ref}\t{new_alt}\t.\tPASS\t{info}\n"
        )


# ============================================================
# Write normalization table
# ============================================================

with open(out_table, "w") as out:

    out.write(
        "CHROM\tPOS\tID\t"
        "Original_REF\tOriginal_ALT\t"
        "Normalized_REF\tNormalized_ALT\t"
        "FASTA_REF\tClassification\n"
    )


    for row in records:

        (
            chrom,
            pos,
            vid,
            old_ref,
            old_alt,
            new_ref,
            new_alt,
            fasta_ref,
            category
        ) = row


        out.write(
            f"{chrom}\t{pos}\t{vid}\t"
            f"{old_ref}\t{old_alt}\t"
            f"{new_ref}\t{new_alt}\t"
            f"{fasta_ref}\t{category}\n"
        )


# ============================================================
# QC
# ============================================================

with open(out_qc, "w") as out:

    out.write(
        "Classification\tN\n"
    )

    for category in [
        "EXACT_REF_MATCH",
        "REF_ALT_SWAP",
        "STRAND_COMPLEMENT",
        "STRAND_COMPLEMENT_PLUS_SWAP"
    ]:

        out.write(
            f"{category}\t{counts.get(category,0)}\n"
        )


    out.write(
        f"TOTAL\t{len(records)}\n"
    )


print(
    "Normalized records:",
    len(records)
)

for k, v in counts.items():

    print(
        k,
        v
    )
PY


# ============================================================
# SORT + BGZIP NORMALIZED VCF
# ============================================================

bcftools sort \
    "${OUTDIR}/VEP_normalized.unsorted.vcf" \
    -Oz \
    -o "${VCF_NORMALIZED}"


rm -f \
    "${OUTDIR}/VEP_normalized.unsorted.vcf"


bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    "${VCF_NORMALIZED}"


N_NORMALIZED=$(
    bcftools view \
        -H \
        "${VCF_NORMALIZED}" \
    | wc -l
)


echo
echo "Normalized VEP records: ${N_NORMALIZED}"


if [[ "${N_NORMALIZED}" -ne "${N_ORIGINAL}" ]]
then

    echo "ERROR: record count changed during normalization."
    exit 1

fi


echo
cat "${NORMALIZATION_QC}"


# ============================================================
# FINAL REF / FASTA VALIDATION
# ============================================================

echo
echo "============================================================"
echo "VALIDATING NORMALIZED VCF AGAINST FASTA"
echo "============================================================"


REF_CHECK_ERR="${OUTDIR}/VEP_normalized_REF_check.err"


set +e


bcftools norm \
    -f "${FASTA}" \
    -c e \
    -Ou \
    "${VCF_NORMALIZED}" \
    >/dev/null \
    2> "${REF_CHECK_ERR}"


REF_STATUS=$?


set -e


if [[ "${REF_STATUS}" -ne 0 ]]
then

    echo "ERROR: normalized VCF still has REF mismatches."

    cat "${REF_CHECK_ERR}"

    exit 1

fi


echo "PASS: all normalized REF alleles match FASTA."


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
# FASTA contigs/order
# ============================================================

order = {}

with open(fai) as fh:

    for i, line in enumerate(fh):

        chrom = line.rstrip("\n").split("\t")[0]

        order[chrom] = i


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


def add_attr(s, key, value):

    attrs = [
        x
        for x in s.split(";")
        if x
    ]

    existing = {
        x.split("=",1)[0]
        for x in attrs
        if "=" in x
    }

    if key not in existing:

        attrs.append(
            f"{key}={value}"
        )

    return ";".join(attrs)


# ============================================================
# Process
# ============================================================

records = []

removed_contigs = Counter()

type_counts = Counter()

n_converted = 0
n_biotype = 0

idx = 0
inside_fasta = False


with open(gff_in) as fh:

    for raw in fh:


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


        f = line.split("\t")


        if len(f) != 9:

            continue


        chrom = f[0]


        if chrom not in order:

            removed_contigs[chrom] += 1
            continue


        if f[2] == "transcript":

            f[2] = "mRNA"

            n_converted += 1


        if f[2] in {
            "gene",
            "mRNA"
        }:

            old = f[8]

            f[8] = add_attr(
                f[8],
                "biotype",
                "protein_coding"
            )

            if old != f[8]:

                n_biotype += 1


        type_counts[
            f[2]
        ] += 1


        start = int(f[3])
        end = int(f[4])


        records.append(
            (
                order[chrom],
                start,
                end,
                idx,
                "\t".join(f)
            )
        )


        idx += 1


records.sort(
    key=lambda x: (
        x[0],
        x[1],
        x[2],
        x[3]
    )
)


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

mrnas = set()

mrna_parents = []

child_parents = []


with open(gff_out) as fh:

    for line in fh:

        if line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")

        attrs = parse_attrs(
            f[8]
        )


        if f[2] == "gene":

            if "ID" in attrs:

                genes.add(
                    attrs["ID"]
                )


        elif f[2] == "mRNA":

            tid = attrs.get(
                "ID"
            )

            parent = attrs.get(
                "Parent"
            )


            if tid:

                mrnas.add(
                    tid
                )


            if tid and parent:

                mrna_parents.append(
                    (
                        tid,
                        parent
                    )
                )


        elif f[2] in {
            "exon",
            "CDS"
        }:

            parent = attrs.get(
                "Parent"
            )


            if parent:

                child_parents.extend(
                    parent.split(",")
                )


bad_mrna = [
    x
    for x in mrna_parents
    if x[1] not in genes
]


bad_child = [
    x
    for x in child_parents
    if x not in mrnas
]


with open(qc_out, "w") as out:

    out.write(
        f"Retained_records\t{len(records)}\n"
    )

    out.write(
        f"Removed_contigs\t{len(removed_contigs)}\n"
    )

    out.write(
        f"Transcript_to_mRNA\t{n_converted}\n"
    )

    out.write(
        f"Biotype_added\t{n_biotype}\n"
    )

    out.write(
        f"Genes\t{len(genes)}\n"
    )

    out.write(
        f"mRNA\t{len(mrnas)}\n"
    )

    out.write(
        f"Invalid_mRNA_gene_links\t{len(bad_mrna)}\n"
    )

    out.write(
        f"Invalid_child_mRNA_links\t{len(bad_child)}\n"
    )


    out.write(
        "\nFeature_types\n"
    )


    for k, v in type_counts.most_common():

        out.write(
            f"{k}\t{v}\n"
        )


print(
    "Retained GFF records:",
    len(records)
)

print(
    "transcript -> mRNA:",
    n_converted
)

print(
    "biotype added:",
    n_biotype
)

print(
    "genes:",
    len(genes)
)

print(
    "mRNAs:",
    len(mrnas)
)

print(
    "invalid mRNA->gene:",
    len(bad_mrna)
)

print(
    "invalid child->mRNA:",
    len(bad_child)
)


if not mrnas:

    raise RuntimeError(
        "No mRNA models generated"
    )


if bad_mrna:

    raise RuntimeError(
        "Broken mRNA -> gene links"
    )


if bad_child:

    raise RuntimeError(
        "Broken exon/CDS -> mRNA links"
    )
PY


cat "${GFF_QC}"


# ============================================================
# COMPRESS / INDEX GFF
# ============================================================

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

    echo "ERROR: fixed GFF index not produced."
    exit 1

fi


# ============================================================
# SELECT AN ANNOTATED SNP FOR SMOKE TEST
#
# Prefer:
#   CDS > exon > mRNA
# ============================================================

echo
echo "============================================================"
echo "SELECTING ANNOTATED SMOKE-TEST SNP"
echo "============================================================"


python - \
    "${VCF_NORMALIZED}" \
    "${GFF_FIXED}" \
    "${TEST_INFO}" <<'PY'

import gzip
import sys


vcf_file = sys.argv[1]
gff_file = sys.argv[2]
outfile = sys.argv[3]


variants = {}


with gzip.open(vcf_file, "rt") as fh:

    for line in fh:

        if line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")

        chrom = f[0]

        variants.setdefault(
            chrom,
            []
        ).append(
            (
                int(f[1]),
                f[3],
                f[4]
            )
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


        if chrom not in variants:

            continue


        start = int(f[3])

        end = int(f[4])


        for pos, ref, alt in variants[chrom]:

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
        "No candidate SNP overlaps CDS/exon/mRNA"
    )


_, chrom, pos, ref, alt, feature, start, end = best


with open(outfile, "w") as out:

    out.write(
        "CHROM\tPOS\tREF\tALT\tFEATURE\tSTART\tEND\n"
    )

    out.write(
        f"{chrom}\t{pos}\t{ref}\t{alt}\t"
        f"{feature}\t{start}\t{end}\n"
    )


print(
    "Selected:",
    chrom,
    pos,
    ref,
    alt,
    "inside",
    feature
)
PY


cat "${TEST_INFO}"


TEST_CHROM=$(
    awk 'NR==2 {print $1}' \
        "${TEST_INFO}"
)


TEST_POS=$(
    awk 'NR==2 {print $2}' \
        "${TEST_INFO}"
)


# ============================================================
# CREATE SMOKE-TEST VCF
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
    "${VCF_NORMALIZED}"


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


if [[ "${N_TEST_VCF}" -eq 0 ]]
then

    echo "ERROR: smoke-test SNP not recovered."
    exit 1

fi


# ============================================================
# RUN SMOKE TEST
# ============================================================

echo
echo "============================================================"
echo "RUNNING VEP SMOKE TEST"
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

    echo "ERROR: smoke-test VEP crashed."

    cat "${TEST_ERR}"

    exit 1

fi


N_TEST_ROWS=$(
    grep -vc '^#' \
        "${TEST_OUT}" \
    || true
)


N_TEST_NON_INTERGENIC=$(
    grep -v '^#' \
        "${TEST_OUT}" \
    | awk -F'\t' '
        $7 != "intergenic_variant" {
            n++
        }
        END {
            print n+0
        }
    '
)


echo
echo "Smoke-test rows:            ${N_TEST_ROWS}"
echo "Smoke-test non-intergenic: ${N_TEST_NON_INTERGENIC}"


grep -v '^##' \
    "${TEST_OUT}" \
    || true


if [[ "${N_TEST_ROWS}" -eq 0 ]]
then

    echo "ERROR: smoke test produced zero VEP rows."
    exit 1

fi


if [[ "${N_TEST_NON_INTERGENIC}" -eq 0 ]]
then

    echo
    echo "ERROR:"
    echo "SNP physically overlaps CDS/exon/mRNA but VEP"
    echo "still reports only intergenic_variant."

    exit 1

fi


echo
echo "SMOKE TEST PASSED."


# ============================================================
# FULL VEP
# ============================================================

echo
echo "============================================================"
echo "RUNNING FULL VEP"
echo "============================================================"


rm -f "${VEP_OUT}"


vep \
    --input_file "${VCF_NORMALIZED}" \
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

    echo "ERROR: VEP output missing."
    exit 1

fi


N_VEP_ROWS=$(
    grep -vc '^#' \
        "${VEP_OUT}" \
    || true
)


N_INTERGENIC=$(
    grep -v '^#' \
        "${VEP_OUT}" \
    | awk -F'\t' '
        $7 == "intergenic_variant" {
            n++
        }
        END {
            print n+0
        }
    '
)


N_NON_INTERGENIC=$(
    grep -v '^#' \
        "${VEP_OUT}" \
    | awk -F'\t' '
        $7 != "intergenic_variant" {
            n++
        }
        END {
            print n+0
        }
    '
)


echo "Normalized VCF records: ${N_NORMALIZED}"
echo "VEP consequence rows:    ${N_VEP_ROWS}"
echo "Intergenic rows:         ${N_INTERGENIC}"
echo "Non-intergenic rows:     ${N_NON_INTERGENIC}"


if [[ "${N_VEP_ROWS}" -eq 0 ]]
then

    echo "ERROR: full VEP produced zero rows."
    exit 1

fi


if [[ "${N_NON_INTERGENIC}" -eq 0 ]]
then

    echo
    echo "ERROR:"
    echo "Full VEP output is entirely intergenic."
    echo "Do not use this annotation."

    exit 1

fi


# ============================================================
# RAW CONSEQUENCE COUNTS
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
echo "RAW VEP CONSEQUENCES"
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
# FINAL QC REPORT
# ============================================================

{
    echo "FINAL VEP QC - Cactoblastis cactorum"
    echo "===================================="
    echo

    echo "Original candidate records"
    echo "--------------------------"
    echo "${N_ORIGINAL}"
    echo

    echo "Normalized candidate records"
    echo "----------------------------"
    echo "${N_NORMALIZED}"
    echo

    echo "REF normalization"
    echo "-----------------"
    cat "${NORMALIZATION_QC}"
    echo

    echo "Normalized REF vs FASTA"
    echo "-----------------------"
    echo "PASS"
    echo

    echo "Fixed GFF"
    echo "---------"
    echo "${GFF_FIXED_GZ}"
    echo

    echo "GFF QC"
    echo "------"
    cat "${GFF_QC}"
    echo

    echo "Smoke-test candidate"
    echo "--------------------"
    cat "${TEST_INFO}"
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
    echo "${N_VEP_ROWS}"
    echo

    echo "Unique VEP locations"
    echo "--------------------"
    echo "${N_LOCATIONS}"
    echo

    echo "Intergenic rows"
    echo "---------------"
    echo "${N_INTERGENIC}"
    echo

    echo "Non-intergenic rows"
    echo "-------------------"
    echo "${N_NON_INTERGENIC}"
    echo

    echo "Raw consequences"
    echo "----------------"
    cat "${RAW_COUNTS}"

} > "${FINAL_QC}"


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
echo "FINAL NORMALIZATION + VEP PASSED"
echo "============================================================"

echo
echo "Normalized VEP VCF:"
echo "${VCF_NORMALIZED}"

echo
echo "Normalization table:"
echo "${NORMALIZATION_TABLE}"

echo
echo "Fixed GFF:"
echo "${GFF_FIXED_GZ}"

echo
echo "Final VEP output:"
echo "${VEP_OUT}"

echo
echo "Raw consequences:"
echo "${RAW_COUNTS}"

echo
echo "Final QC:"
echo "${FINAL_QC}"

echo
echo "DONE"
date

