#!/usr/bin/env bash
#SBATCH --job-name=blast_g11528
#SBATCH --output=logs/blast_g11528_%j.out
#SBATCH --error=logs/blast_g11528_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00

set -euo pipefail


#######################################################################
# Functional annotation of Ccac_g11528
#
# Runs:
#   1. Full protein vs nr
#   2. Full protein vs RefSeq protein
#   3. Conserved/core region vs nr
#   4. Conserved/core region vs RefSeq protein
#
# Cactoblastis cactorum
#######################################################################


# =====================================================================
# MODULES
# =====================================================================

module purge
module load bioinfo/NCBI_Blast+/2.15.0+
module load bioinfo/samtools/1.23


# =====================================================================
# PATHS
# =====================================================================

WORK="/home/dpoveda/micetral/Population_snps/baypass"

QUERY="${WORK}/Ccac_g11528.protein.fa"

OUTDIR="${WORK}/selection_scan/Ccac_g11528_BLAST"

mkdir -p "${OUTDIR}"
mkdir -p "${WORK}/logs"


FULL_NR="${OUTDIR}/Ccac_g11528.full.blastp_nr.tsv"
FULL_REFSEQ="${OUTDIR}/Ccac_g11528.full.blastp_refseq.tsv"

CORE="${OUTDIR}/Ccac_g11528.core_170_290.fa"

CORE_NR="${OUTDIR}/Ccac_g11528.core.blastp_nr.tsv"
CORE_REFSEQ="${OUTDIR}/Ccac_g11528.core.blastp_refseq.tsv"

SUMMARY_FULL_NR="${OUTDIR}/Ccac_g11528.full.blastp_nr.TOP50.tsv"
SUMMARY_FULL_REFSEQ="${OUTDIR}/Ccac_g11528.full.blastp_refseq.TOP50.tsv"

SUMMARY_CORE_NR="${OUTDIR}/Ccac_g11528.core.blastp_nr.TOP50.tsv"
SUMMARY_CORE_REFSEQ="${OUTDIR}/Ccac_g11528.core.blastp_refseq.TOP50.tsv"


# =====================================================================
# CHECK SOFTWARE
# =====================================================================

for prog in blastp blastdbcmd samtools python
do
    if ! command -v "${prog}" >/dev/null 2>&1
    then
        echo "ERROR: ${prog} not found."
        exit 1
    fi
done


echo
echo "============================================================"
echo "Ccac_g11528 FUNCTIONAL ANNOTATION BY BLAST"
echo "============================================================"
echo

echo "BLAST:"
blastp -version | head -2

echo
echo "BLASTDB=${BLASTDB:-NOT_SET}"
echo


# =====================================================================
# CHECK QUERY
# =====================================================================

if [[ ! -s "${QUERY}" ]]
then
    echo "ERROR: query protein not found:"
    echo "${QUERY}"
    exit 1
fi


echo "Query:"
echo "${QUERY}"

echo
cat "${QUERY}"

echo


# =====================================================================
# PROTEIN LENGTH
# =====================================================================

PROTEIN_LENGTH=$(
    grep -v '^>' "${QUERY}" |
    tr -d '\n\r ' |
    awk '{print length($0)}'
)

echo "Protein length: ${PROTEIN_LENGTH} aa"


if [[ "${PROTEIN_LENGTH}" -ne 408 ]]
then
    echo "WARNING: expected 408 aa, observed ${PROTEIN_LENGTH} aa."
fi


# =====================================================================
# CHECK DATABASES
# =====================================================================

echo
echo "Checking nr..."

if ! blastdbcmd -db nr -info >/dev/null 2>&1
then
    echo "ERROR: nr database unavailable."
    exit 1
fi

echo "PASS: nr available."


echo
echo "Checking refseq_protein..."

if ! blastdbcmd -db refseq_protein -info >/dev/null 2>&1
then
    echo "ERROR: refseq_protein database unavailable."
    exit 1
fi

echo "PASS: refseq_protein available."


# =====================================================================
# CREATE CENTRAL / LOWER-COMPLEXITY CORE
#
# Residues 170-290 inclusive.
# This is NOT treated as an inferred domain.
# It is simply a complementary search of the less repetitive central
# region of the predicted protein.
# =====================================================================

echo
echo "Creating residues 170-290 query..."


python - "${QUERY}" "${CORE}" <<'PY'

import sys

query, output = sys.argv[1:]

seq = ""

with open(query) as fh:
    for line in fh:
        line = line.strip()

        if not line or line.startswith(">"):
            continue

        seq += line


start = 170
end = 290

if len(seq) < end:
    raise SystemExit(
        f"Protein only has {len(seq)} aa; cannot extract {start}-{end}."
    )


fragment = seq[start - 1:end]


with open(output, "w") as out:
    out.write(">Ccac_g11528_core_170_290\n")

    for i in range(0, len(fragment), 60):
        out.write(fragment[i:i+60] + "\n")


print(
    f"Full protein length = {len(seq)} aa"
)

print(
    f"Core coordinates = {start}-{end}"
)

print(
    f"Core length = {len(fragment)} aa"
)

print(
    f"Core sequence = {fragment}"
)

PY


# =====================================================================
# BLAST OUTPUT FORMAT
# =====================================================================

OUTFMT="6 qseqid saccver pident length qlen slen qstart qend sstart send evalue bitscore stitle staxids"


# =====================================================================
# 1. FULL PROTEIN vs nr
# =====================================================================

echo
echo "============================================================"
echo "1/4 FULL PROTEIN vs nr"
echo "============================================================"
echo


blastp \
    -query "${QUERY}" \
    -db nr \
    -evalue 1e-5 \
    -max_target_seqs 200 \
    -num_threads "${SLURM_CPUS_PER_TASK:-8}" \
    -seg yes \
    -comp_based_stats 2 \
    -outfmt "${OUTFMT}" \
    -out "${FULL_NR}"


echo "Finished full protein vs nr."

echo


# =====================================================================
# 2. FULL PROTEIN vs RefSeq protein
# =====================================================================

echo
echo "============================================================"
echo "2/4 FULL PROTEIN vs refseq_protein"
echo "============================================================"
echo


blastp \
    -query "${QUERY}" \
    -db refseq_protein \
    -evalue 1e-5 \
    -max_target_seqs 200 \
    -num_threads "${SLURM_CPUS_PER_TASK:-8}" \
    -seg yes \
    -comp_based_stats 2 \
    -outfmt "${OUTFMT}" \
    -out "${FULL_REFSEQ}"


echo "Finished full protein vs RefSeq."

echo


# =====================================================================
# 3. CORE vs nr
# =====================================================================

echo
echo "============================================================"
echo "3/4 CORE 170-290 vs nr"
echo "============================================================"
echo


blastp \
    -query "${CORE}" \
    -db nr \
    -evalue 1e-3 \
    -max_target_seqs 200 \
    -num_threads "${SLURM_CPUS_PER_TASK:-8}" \
    -seg yes \
    -comp_based_stats 2 \
    -outfmt "${OUTFMT}" \
    -out "${CORE_NR}"


echo "Finished core vs nr."

echo


# =====================================================================
# 4. CORE vs RefSeq protein
# =====================================================================

echo
echo "============================================================"
echo "4/4 CORE 170-290 vs refseq_protein"
echo "============================================================"
echo


blastp \
    -query "${CORE}" \
    -db refseq_protein \
    -evalue 1e-3 \
    -max_target_seqs 200 \
    -num_threads "${SLURM_CPUS_PER_TASK:-8}" \
    -seg yes \
    -comp_based_stats 2 \
    -outfmt "${OUTFMT}" \
    -out "${CORE_REFSEQ}"


echo "Finished core vs RefSeq."

echo


# =====================================================================
# MAKE HUMAN-READABLE SUMMARY TABLES
#
# Input columns:
#
# 1 qseqid
# 2 saccver
# 3 pident
# 4 alignment length
# 5 qlen
# 6 slen
# 7 qstart
# 8 qend
# 9 sstart
# 10 send
# 11 evalue
# 12 bitscore
# 13 stitle
# 14 staxids
#
# Approximate query coverage:
# 100 * aligned length / query length
# =====================================================================

make_summary () {

    INPUT="$1"
    OUTPUT="$2"

    {
        printf \
"Accession\tPercent_identity\tAlignment_length\tQuery_length\tQuery_coverage_percent\tSubject_length\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tEvalue\tBitscore\tDescription\tTaxIDs\n"

        awk -F'\t' '
        BEGIN {
            OFS="\t"
        }
        {
            coverage = 100 * $4 / $5

            print \
                $2, \
                $3, \
                $4, \
                $5, \
                sprintf("%.2f", coverage), \
                $6, \
                $7, \
                $8, \
                $9, \
                $10, \
                $11, \
                $12, \
                $13, \
                $14
        }
        ' "${INPUT}" \
        | sort -t $'\t' -k12,12gr \
        | head -50

    } > "${OUTPUT}"
}


make_summary \
    "${FULL_NR}" \
    "${SUMMARY_FULL_NR}"


make_summary \
    "${FULL_REFSEQ}" \
    "${SUMMARY_FULL_REFSEQ}"


make_summary \
    "${CORE_NR}" \
    "${SUMMARY_CORE_NR}"


make_summary \
    "${CORE_REFSEQ}" \
    "${SUMMARY_CORE_REFSEQ}"


# =====================================================================
# QUICK CONSENSUS DESCRIPTION COUNTS
# =====================================================================

echo
echo "============================================================"
echo "TOP HIT DESCRIPTIONS"
echo "============================================================"
echo


for FILE in \
    "${SUMMARY_FULL_NR}" \
    "${SUMMARY_FULL_REFSEQ}" \
    "${SUMMARY_CORE_NR}" \
    "${SUMMARY_CORE_REFSEQ}"
do

    echo
    echo "------------------------------------------------------------"
    echo "$(basename "${FILE}")"
    echo "------------------------------------------------------------"

    head -11 "${FILE}"

done


# =====================================================================
# FINAL QC
# =====================================================================

echo
echo "============================================================"
echo "FINAL QC"
echo "============================================================"
echo


for FILE in \
    "${FULL_NR}" \
    "${FULL_REFSEQ}" \
    "${CORE_NR}" \
    "${CORE_REFSEQ}"
do

    N=$(wc -l < "${FILE}")

    echo "$(basename "${FILE}") : ${N} hits"

done


echo
echo "============================================================"
echo "SUCCESS"
echo "============================================================"
echo

echo "Results directory:"
echo "${OUTDIR}"

echo
echo "Most useful files:"
echo "${SUMMARY_FULL_NR}"
echo "${SUMMARY_FULL_REFSEQ}"
echo "${SUMMARY_CORE_NR}"
echo "${SUMMARY_CORE_REFSEQ}"
