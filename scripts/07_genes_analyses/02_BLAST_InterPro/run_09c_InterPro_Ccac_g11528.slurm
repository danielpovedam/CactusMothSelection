#!/usr/bin/env bash
#SBATCH --job-name=ipr_g11528
#SBATCH --output=logs/ipr_g11528_%j.out
#SBATCH --error=logs/ipr_g11528_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00

set -euo pipefail

module purge
#module load bioinfo/InterProScan/5.73-104.0
module load devel/java/17.0.6 
module load devel/python/Python-3.11.1
module load bioinfo/InterProScan/5.73-104.0

WORK="/home/dpoveda/micetral/Population_snps/baypass"
QUERY="${WORK}/Ccac_g11528.protein.fa"

OUTDIR="${WORK}/selection_scan/Ccac_g11528_InterPro"

mkdir -p "${OUTDIR}"
mkdir -p "${WORK}/logs"

PREFIX="${OUTDIR}/Ccac_g11528.interpro"

echo
echo "============================================================"
echo "INTERPROSCAN: Ccac_g11528"
echo "============================================================"
echo

if [[ ! -s "${QUERY}" ]]
then
    echo "ERROR: query not found:"
    echo "${QUERY}"
    exit 1
fi

echo "Query:"
echo "${QUERY}"

echo
echo "Protein length:"
grep -v '^>' "${QUERY}" \
| tr -d '\n\r ' \
| awk '{print length($0) " aa"}'

echo

if command -v interproscan.sh >/dev/null 2>&1
then
    IPR="interproscan.sh"
elif command -v interproscan >/dev/null 2>&1
then
    IPR="interproscan"
else
    echo "ERROR: InterProScan executable not found."
    exit 1
fi

echo "InterProScan:"
echo "${IPR}"
echo

rm -f \
    "${PREFIX}.tsv" \
    "${PREFIX}.gff3" \
    "${PREFIX}.json" \
    "${PREFIX}.xml"

"${IPR}" \
    -i "${QUERY}" \
    -b "${PREFIX}" \
    -f TSV,GFF3,JSON,XML \
    -cpu "${SLURM_CPUS_PER_TASK:-8}" \
    -goterms \
    -iprlookup \
    -pa \
    -dp

echo
echo "============================================================"
echo "RESULTS"
echo "============================================================"

TSV="${PREFIX}.tsv"

if [[ ! -s "${TSV}" ]]
then
    echo "ERROR: TSV output not found:"
    echo "${TSV}"
    exit 1
fi

echo
echo "Annotation records:"
wc -l "${TSV}"

echo
echo "Top InterProScan records:"
head -30 "${TSV}"

echo
echo "Distinct databases:"
cut -f4 "${TSV}" | sort | uniq -c | sort -nr

echo
echo "Distinct InterPro annotations:"
awk -F'\t' '
    NF >= 12 && $12 != "-" && $12 != "" {
        print $12 "\t" $13
    }
' "${TSV}" \
| sort -u

echo
echo "GO terms:"
awk -F'\t' '
    NF >= 14 && $14 != "-" && $14 != "" {
        print $14
    }
' "${TSV}" \
| sort -u

echo
echo "============================================================"
echo "SUCCESS"
echo "============================================================"

echo
echo "Output prefix:"
echo "${PREFIX}"
