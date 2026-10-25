#!/bin/bash
#SBATCH --job-name=Ccac_BLAST_candidates
#SBATCH --output=logs/Ccac_BLAST_candidates_%j.out
#SBATCH --error=logs/Ccac_BLAST_candidates_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail

###############################################################################
# Cactoblastis cactorum
# BLAST refinement of incompletely annotated XtX, C2 and GEA candidate genes
#
# Author: Daniel Poveda-Martinez
#
# PURPOSE
# -------
# 1. Identify candidate genes whose reference-genome annotation is incomplete:
#      - hypothetical/uncharacterized description
#      - empty gene symbol
#      - gene symbol identical to Gene ID
#
# 2. Explicitly include Ccac_g2721, associated with the TempPC1 locus
#    ptg000019:203701-203719 shown in Figure 4D.
#
# 3. Extract the longest predicted protein isoform per candidate gene.
#
# 4. Search predicted proteins against:
#      - NCBI nr
#      - NCBI RefSeq protein
#
# 5. Select the best-supported BLASTP hit per gene.
#
# 6. Append BLAST evidence to the final XtX, C2 and GEA candidate tables.
#
# BLAST+ version: 2.15.0+
# samtools version: 1.23
###############################################################################

###############################################################################
# 0. PATHS
###############################################################################

BASE="/home/dpoveda/micetral/Population_snps/baypass"

WORKDIR="${BASE}/selection_scan/CANDIDATE_BLAST_REANNOTATION"

PROTEINS="${BASE}/Ccactorum_predicted_proteins.fa"

XTX_TABLE="${BASE}/selection_scan/FINAL_gt1Mb/Table_XtX_candidate_regions_allSNP_gt1Mb_FINAL.tsv"

C2_TABLE="${BASE}/selection_scan/FINAL_gt1Mb/Table_C2_Host_candidate_regions_allSNP_gt1Mb_FINAL.tsv"

GEA_TABLE="${BASE}/selection_scan/GEA_SUPPLEMENTARY_PUBLICATION_FINAL/Table_Sx_ALL_167_GEA_candidate_loci_FINAL.tsv"

# BLAST databases are exposed through BLASTDB=/bank/blastdb
NR_DB="nr"
REFSEQ_DB="refseq_protein"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

###############################################################################
# 1. MODULES
###############################################################################

module purge
module load bioinfo/NCBI_Blast+/2.15.0+
module load bioinfo/samtools/1.23

echo "=============================================================="
echo "STARTING CANDIDATE BLAST REANNOTATION"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "=============================================================="
echo

echo "BLAST version:"
blastp -version
echo

echo "samtools version:"
samtools --version | head -1
echo

echo "BLASTDB=${BLASTDB:-NOT_SET}"
echo

###############################################################################
# 2. CREATE OUTPUT DIRECTORIES
###############################################################################

mkdir -p "${WORKDIR}"
mkdir -p "${WORKDIR}/01_gene_lists"
mkdir -p "${WORKDIR}/02_proteins"
mkdir -p "${WORKDIR}/03_blast"
mkdir -p "${WORKDIR}/04_summary"
mkdir -p "${WORKDIR}/05_tables"
mkdir -p "${WORKDIR}/06_QC"
mkdir -p "${BASE}/logs"

cd "${WORKDIR}"

###############################################################################
# 3. CHECK INPUT FILES
###############################################################################

echo "Checking input files..."
echo

for f in \
    "${PROTEINS}" \
    "${XTX_TABLE}" \
    "${C2_TABLE}" \
    "${GEA_TABLE}"
do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty input:"
        echo "${f}"
        exit 1
    fi
    echo "PASS: ${f}"
done

echo

###############################################################################
# 4. CHECK BLAST DATABASES
###############################################################################

echo "Checking BLAST databases..."
echo

if blastdbcmd -db "${NR_DB}" -info >/dev/null 2>&1; then
    echo "PASS: nr available."
else
    echo "ERROR: nr database unavailable."
    exit 1
fi

if blastdbcmd -db "${REFSEQ_DB}" -info >/dev/null 2>&1; then
    echo "PASS: refseq_protein available."
else
    echo "ERROR: refseq_protein database unavailable."
    exit 1
fi

echo

###############################################################################
# 5. IDENTIFY INCOMPLETELY ANNOTATED CANDIDATE GENES
###############################################################################

echo "Identifying incompletely annotated candidate genes..."
echo

export XTX_TABLE
export C2_TABLE
export GEA_TABLE

python3 <<'PY'
import os
import pandas as pd
from pathlib import Path

tables = {
    "XtX": Path(os.environ["XTX_TABLE"]),
    "C2": Path(os.environ["C2_TABLE"]),
    "GEA": Path(os.environ["GEA_TABLE"]),
}

OUTDIR = Path("01_gene_lists")
OUTDIR.mkdir(exist_ok=True)

def normalise(x):
    return (
        str(x)
        .strip()
        .replace("\\_", "_")
        .replace("*", "")
    )

def find_col(df, candidates):
    norm = {
        str(c).strip().lower()
        .replace(" ", "_")
        .replace("\\_", "_"): c
        for c in df.columns
    }

    for candidate in candidates:
        key = candidate.lower().replace(" ", "_")
        if key in norm:
            return norm[key]

    return None

records = []

for analysis, path in tables.items():

    print(f"\nReading {analysis}: {path}")

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )

    gene_col = find_col(
        df,
        ["Gene ID", "Gene_ID", "gene_id", "GeneID"]
    )

    symbol_col = find_col(
        df,
        ["Gene symbol", "Gene_symbol", "gene_symbol", "Symbol"]
    )

    desc_col = find_col(
        df,
        ["Description", "description"]
    )

    if gene_col is None:
        raise RuntimeError(
            f"Could not identify Gene ID column in {path}\n"
            f"Columns found: {list(df.columns)}"
        )

    for _, row in df.iterrows():

        gene = normalise(row.get(gene_col, ""))

        if gene in {"", "-", "NA", "na", "None", "nan"}:
            continue

        symbol = normalise(row.get(symbol_col, "")) if symbol_col else ""
        desc = str(row.get(desc_col, "")).strip() if desc_col else ""

        desc_lower = desc.lower()

        incomplete = (
            symbol in {"", "-"} or
            symbol == gene or
            "hypothetical" in desc_lower or
            "uncharacterized" in desc_lower or
            "uncharacterised" in desc_lower or
            "unknown protein" in desc_lower
        )

        if incomplete:
            records.append({
                "Analysis": analysis,
                "Gene_ID": gene,
                "Current_symbol": symbol,
                "Current_description": desc
            })

result = pd.DataFrame(records)

if result.empty:
    raise RuntimeError(
        "No incompletely annotated candidate genes were detected."
    )

result = result.drop_duplicates()

# Explicitly include Figure 4D candidate
if "Ccac_g2721" not in set(result["Gene_ID"]):

    result = pd.concat([
        result,
        pd.DataFrame([{
            "Analysis": "GEA",
            "Gene_ID": "Ccac_g2721",
            "Current_symbol": "",
            "Current_description": "Alpha-mannosidase 2"
        }])
    ], ignore_index=True)

result = result.drop_duplicates()

result.to_csv(
    OUTDIR / "candidate_genes_requiring_BLAST_by_analysis.tsv",
    sep="\t",
    index=False
)

genes = sorted(result["Gene_ID"].unique())

with open(
    OUTDIR / "candidate_genes_requiring_BLAST.txt",
    "w"
) as f:
    for gene in genes:
        f.write(gene + "\n")

print("\n========================================")
print("BLAST candidate summary")
print("========================================")
print(f"Unique genes requiring BLAST: {len(genes)}")
print()

print("Unique genes by analysis:")
print(
    result.groupby("Analysis")["Gene_ID"]
    .nunique()
    .sort_index()
)
print()

print("Genes:")
for gene in genes:
    print(gene)

print()

if "Ccac_g2721" in genes:
    print("PASS: Ccac_g2721 included.")
else:
    raise RuntimeError("Ccac_g2721 was not included.")

PY

###############################################################################
# 6. SHOW GENE LIST
###############################################################################

echo
echo "=============================================================="
echo "GENES TO BLAST"
echo "=============================================================="

cat 01_gene_lists/candidate_genes_requiring_BLAST.txt

echo
echo "Total:"
wc -l 01_gene_lists/candidate_genes_requiring_BLAST.txt

###############################################################################
# 7. EXTRACT LONGEST PROTEIN ISOFORM FOR EACH GENE
###############################################################################

echo
echo "Extracting longest predicted protein isoform per gene..."
echo

export PROTEINS

python3 <<'PY'
import os
from pathlib import Path

protein_fasta = Path(os.environ["PROTEINS"])
gene_file = Path(
    "01_gene_lists/candidate_genes_requiring_BLAST.txt"
)

genes = {
    line.strip()
    for line in gene_file.read_text().splitlines()
    if line.strip()
}

gene_records = {gene: [] for gene in genes}

current_header = None
current_seq = []

def save_record(header, seq_lines):

    if header is None:
        return

    seq = "".join(seq_lines).replace("*", "")
    token = header.split()[0]

    for gene in genes:

        if (
            token == gene or
            token.startswith(gene + ".") or
            token.startswith(gene + "-") or
            token.startswith(gene + "|")
        ):
            gene_records[gene].append((token, seq))
            break

with open(protein_fasta) as f:

    for line in f:

        line = line.rstrip()

        if line.startswith(">"):

            save_record(
                current_header,
                current_seq
            )

            current_header = line[1:]
            current_seq = []

        else:
            current_seq.append(line)

save_record(
    current_header,
    current_seq
)

out_fasta = Path(
    "02_proteins/candidate_longest_isoform.fa"
)

out_map = Path(
    "02_proteins/candidate_longest_isoform_mapping.tsv"
)

missing_file = Path(
    "06_QC/missing_candidate_proteins.txt"
)

missing = []

with open(out_fasta, "w") as out, \
     open(out_map, "w") as mapping:

    mapping.write(
        "Gene_ID\tTranscript_ID\tProtein_length_aa\n"
    )

    for gene in sorted(genes):

        records = gene_records.get(gene, [])

        if not records:
            missing.append(gene)
            continue

        transcript, protein = max(
            records,
            key=lambda x: len(x[1])
        )

        out.write(
            f">{gene}|{transcript}\n"
        )

        for i in range(0, len(protein), 60):
            out.write(
                protein[i:i+60] + "\n"
            )

        mapping.write(
            f"{gene}\t{transcript}\t{len(protein)}\n"
        )

with open(missing_file, "w") as f:
    for gene in missing:
        f.write(gene + "\n")

print(f"Genes requested: {len(genes)}")
print(
    "Genes with extracted protein:",
    len(genes) - len(missing)
)
print(
    "Genes without protein:",
    len(missing)
)

if missing:
    print("\nWARNING: missing proteins:")
    for gene in missing:
        print(gene)

if gene_records.get("Ccac_g2721"):
    print("\nPASS: Ccac_g2721 protein found.")
else:
    print("\nWARNING: Ccac_g2721 protein not found.")

PY

echo
echo "Protein mapping:"
column -t -s $'\t' \
    02_proteins/candidate_longest_isoform_mapping.tsv \
    | head -100

echo

###############################################################################
# 8. CHECK EXTRACTED FASTA
###############################################################################

NSEQ=$(grep -c '^>' \
    02_proteins/candidate_longest_isoform.fa)

echo "Number of proteins in BLAST FASTA: ${NSEQ}"

if [[ "${NSEQ}" -eq 0 ]]; then
    echo "ERROR: no protein sequences extracted."
    exit 1
fi

###############################################################################
# 9. BLASTP AGAINST nr
###############################################################################

echo
echo "=============================================================="
echo "BLASTP AGAINST nr"
echo "=============================================================="
echo

blastp \
    -query 02_proteins/candidate_longest_isoform.fa \
    -db "${NR_DB}" \
    -evalue 1e-5 \
    -max_target_seqs 100 \
    -num_threads "${THREADS}" \
    -seg yes \
    -comp_based_stats 2 \
    -outfmt "6 qseqid saccver pident length qlen slen qstart qend sstart send evalue bitscore stitle staxids" \
    -out 03_blast/candidate_proteins.blastp_nr.tsv

echo "Finished nr BLAST."
echo

###############################################################################
# 10. BLASTP AGAINST refseq_protein
###############################################################################

echo
echo "=============================================================="
echo "BLASTP AGAINST refseq_protein"
echo "=============================================================="
echo

blastp \
    -query 02_proteins/candidate_longest_isoform.fa \
    -db "${REFSEQ_DB}" \
    -evalue 1e-5 \
    -max_target_seqs 100 \
    -num_threads "${THREADS}" \
    -seg yes \
    -comp_based_stats 2 \
    -outfmt "6 qseqid saccver pident length qlen slen qstart qend sstart send evalue bitscore stitle staxids" \
    -out 03_blast/candidate_proteins.blastp_refseq.tsv

echo "Finished RefSeq BLAST."
echo

###############################################################################
# 11. PARSE BLAST RESULTS AND CHOOSE BEST HIT PER GENE
###############################################################################

echo "Parsing BLAST results..."
echo

python3 <<'PY'
import pandas as pd
from pathlib import Path
import numpy as np

columns = [
    "qseqid",
    "saccver",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "stitle",
    "staxids"
]

files = {
    "nr":
        Path(
            "03_blast/"
            "candidate_proteins.blastp_nr.tsv"
        ),

    "RefSeq":
        Path(
            "03_blast/"
            "candidate_proteins.blastp_refseq.tsv"
        )
}

all_hits = []

for database, path in files.items():

    if not path.exists() or path.stat().st_size == 0:
        print(
            f"WARNING: no hits in {database}"
        )
        continue

    df = pd.read_csv(
        path,
        sep="\t",
        names=columns,
        header=None
    )

    df["BLAST_database"] = database

    df["Gene_ID"] = (
        df["qseqid"]
        .astype(str)
        .str.split("|", regex=False)
        .str[0]
    )

    for c in [
        "pident",
        "length",
        "qlen",
        "evalue",
        "bitscore"
    ]:
        df[c] = pd.to_numeric(
            df[c],
            errors="coerce"
        )

    # Query coverage based on aligned residues / query length
    df["query_coverage_pct"] = (
        df["length"] /
        df["qlen"] *
        100
    )

    # Alignment length can include gaps.
    # Cap at 100 for reporting.
    df["query_coverage_pct"] = (
        df["query_coverage_pct"]
        .clip(upper=100)
    )

    all_hits.append(df)

if not all_hits:
    raise RuntimeError(
        "No BLAST hits were obtained from either database."
    )

blast = pd.concat(
    all_hits,
    ignore_index=True
)

###############################################################################
# Save all hits
###############################################################################

blast.to_csv(
    "04_summary/ALL_BLAST_hits.tsv",
    sep="\t",
    index=False
)

###############################################################################
# Best hit separately for each database
###############################################################################

best_db = (
    blast
    .sort_values(
        [
            "Gene_ID",
            "BLAST_database",
            "evalue",
            "bitscore",
            "query_coverage_pct",
            "pident"
        ],
        ascending=[
            True,
            True,
            True,
            False,
            False,
            False
        ]
    )
    .groupby(
        ["Gene_ID", "BLAST_database"],
        as_index=False
    )
    .first()
)

best_db.to_csv(
    "04_summary/"
    "BEST_hit_per_gene_per_database.tsv",
    sep="\t",
    index=False
)

###############################################################################
# Best overall hit
###############################################################################

best = (
    blast
    .sort_values(
        [
            "Gene_ID",
            "evalue",
            "bitscore",
            "query_coverage_pct",
            "pident"
        ],
        ascending=[
            True,
            True,
            False,
            False,
            False
        ]
    )
    .groupby(
        "Gene_ID",
        as_index=False
    )
    .first()
)

final = best[[
    "Gene_ID",
    "BLAST_database",
    "saccver",
    "stitle",
    "pident",
    "query_coverage_pct",
    "evalue",
    "bitscore",
    "staxids"
]].copy()

final.columns = [
    "Gene_ID",
    "BLASTP_database",
    "BLASTP_accession",
    "BLASTP_best_hit",
    "BLASTP_identity_percent",
    "BLASTP_query_coverage_percent",
    "BLASTP_Evalue",
    "BLASTP_bitscore",
    "BLASTP_taxid"
]

###############################################################################
# Format numeric columns
###############################################################################

final[
    "BLASTP_identity_percent"
] = final[
    "BLASTP_identity_percent"
].map(
    lambda x: f"{x:.2f}"
)

final[
    "BLASTP_query_coverage_percent"
] = final[
    "BLASTP_query_coverage_percent"
].map(
    lambda x: f"{x:.2f}"
)

final[
    "BLASTP_Evalue"
] = final[
    "BLASTP_Evalue"
].map(
    lambda x: f"{x:.3e}"
)

final[
    "BLASTP_bitscore"
] = final[
    "BLASTP_bitscore"
].map(
    lambda x: f"{x:.1f}"
)

final.to_csv(
    "04_summary/"
    "BEST_BLAST_hit_per_candidate_gene.tsv",
    sep="\t",
    index=False
)

print(
    "\nBest overall hit per candidate gene:\n"
)

print(
    final.to_string(index=False)
)

###############################################################################
# Figure 4D candidate
###############################################################################

print(
    "\n========================================"
)
print(
    "FIGURE 4D: Ccac_g2721"
)
print(
    "========================================"
)

x = final[
    final["Gene_ID"] == "Ccac_g2721"
]

if len(x) == 0:
    print(
        "WARNING: no BLAST hit for Ccac_g2721"
    )
else:
    print(
        x.to_string(index=False)
    )

PY

###############################################################################
# 12. CREATE A COMPACT BLAST ANNOTATION TABLE
###############################################################################

python3 <<'PY'
import pandas as pd

df = pd.read_csv(
    "04_summary/"
    "BEST_BLAST_hit_per_candidate_gene.tsv",
    sep="\t",
    dtype=str,
    keep_default_na=False
)

compact = df[[
    "Gene_ID",
    "BLASTP_best_hit",
    "BLASTP_identity_percent",
    "BLASTP_query_coverage_percent",
    "BLASTP_Evalue",
    "BLASTP_database"
]]

compact.to_csv(
    "04_summary/"
    "BLAST_annotation_compact.tsv",
    sep="\t",
    index=False
)

PY

###############################################################################
# 13. APPEND BLAST COLUMNS TO XtX, C2 AND GEA TABLES
###############################################################################

echo
echo "Adding BLAST evidence to candidate tables..."
echo

export XTX_TABLE
export C2_TABLE
export GEA_TABLE

python3 <<'PY'
import os
from pathlib import Path
import pandas as pd

tables = {
    "XtX": Path(os.environ["XTX_TABLE"]),
    "C2": Path(os.environ["C2_TABLE"]),
    "GEA": Path(os.environ["GEA_TABLE"])
}

blast = pd.read_csv(
    "04_summary/"
    "BEST_BLAST_hit_per_candidate_gene.tsv",
    sep="\t",
    dtype=str,
    keep_default_na=False
)

def find_col(df, candidates):

    mapping = {
        str(c).strip().lower()
        .replace(" ", "_")
        .replace("\\_", "_"): c
        for c in df.columns
    }

    for candidate in candidates:

        key = (
            candidate
            .lower()
            .replace(" ", "_")
        )

        if key in mapping:
            return mapping[key]

    return None

blast_columns = [
    "BLASTP_best_hit",
    "BLASTP_identity_percent",
    "BLASTP_query_coverage_percent",
    "BLASTP_Evalue",
    "BLASTP_database"
]

for analysis, path in tables.items():

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )

    gene_col = find_col(
        df,
        [
            "Gene ID",
            "Gene_ID",
            "gene_id",
            "GeneID"
        ]
    )

    if gene_col is None:
        raise RuntimeError(
            f"No Gene ID column in {path}"
        )

    df["_Gene_ID_merge"] = (
        df[gene_col]
        .astype(str)
        .str.replace(
            "\\_",
            "_",
            regex=False
        )
        .str.replace(
            "*",
            "",
            regex=False
        )
        .str.strip()
    )

    b = blast.rename(
        columns={
            "Gene_ID":
            "_Gene_ID_merge"
        }
    )

    # Keep only publication-oriented compact columns
    b = b[[
        "_Gene_ID_merge",
        "BLASTP_best_hit",
        "BLASTP_identity_percent",
        "BLASTP_query_coverage_percent",
        "BLASTP_Evalue",
        "BLASTP_database"
    ]]

    merged = df.merge(
        b,
        on="_Gene_ID_merge",
        how="left"
    )

    merged = merged.drop(
        columns="_Gene_ID_merge"
    )

    for col in blast_columns:
        merged[col] = (
            merged[col]
            .fillna("")
        )

    outfile = (
        Path("05_tables") /
        f"{analysis}_candidate_table_with_BLAST.tsv"
    )

    merged.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    nannot = (
        merged["BLASTP_Evalue"] != ""
    ).sum()

    print(
        f"{analysis}: "
        f"{len(merged)} rows; "
        f"{nannot} rows with BLAST evidence"
    )

PY

###############################################################################
# 14. DEDICATED FIGURE 4D REPORT
###############################################################################

echo
echo "=============================================================="
echo "FIGURE 4D LOCUS REPORT"
echo "TempPC1 -- ptg000019 -- Ccac_g2721"
echo "=============================================================="
echo

python3 <<'PY'
import pandas as pd

gea = pd.read_csv(
    "05_tables/"
    "GEA_candidate_table_with_BLAST.tsv",
    sep="\t",
    dtype=str,
    keep_default_na=False
)

def get_col(df, possibilities):

    mapping = {
        str(c).strip().lower()
        .replace(" ", "_")
        .replace("\\_", "_"): c
        for c in df.columns
    }

    for p in possibilities:
        key = (
            p.lower()
            .replace(" ", "_")
        )
        if key in mapping:
            return mapping[key]

    return None

cov_col = get_col(
    gea,
    ["Covariate"]
)

scaf_col = get_col(
    gea,
    ["Scaffold"]
)

gene_col = get_col(
    gea,
    ["Gene ID", "Gene_ID"]
)

if cov_col is None or \
   scaf_col is None or \
   gene_col is None:

    raise RuntimeError(
        "Could not identify required GEA columns."
    )

gene_clean = (
    gea[gene_col]
    .str.replace(
        "\\_",
        "_",
        regex=False
    )
    .str.replace(
        "*",
        "",
        regex=False
    )
)

subset = gea[
    (gea[cov_col] == "TempPC1") &
    (gea[scaf_col] == "ptg000019") &
    (gene_clean == "Ccac_g2721")
]

subset.to_csv(
    "06_QC/"
    "Figure4D_Ccac_g2721_annotation.tsv",
    sep="\t",
    index=False
)

print(
    subset.to_string(index=False)
)

PY

###############################################################################
# 15. QC SUMMARY
###############################################################################

echo
echo "=============================================================="
echo "QC SUMMARY"
echo "=============================================================="
echo

python3 <<'PY'
import pandas as pd

genes = pd.read_csv(
    "01_gene_lists/"
    "candidate_genes_requiring_BLAST_by_analysis.tsv",
    sep="\t"
)

mapping = pd.read_csv(
    "02_proteins/"
    "candidate_longest_isoform_mapping.tsv",
    sep="\t"
)

best = pd.read_csv(
    "04_summary/"
    "BEST_BLAST_hit_per_candidate_gene.tsv",
    sep="\t"
)

print(
    "Candidate genes selected:",
    genes["Gene_ID"].nunique()
)

print(
    "Proteins extracted:",
    mapping["Gene_ID"].nunique()
)

print(
    "Genes with >=1 BLAST hit:",
    best["Gene_ID"].nunique()
)

selected = set(
    genes["Gene_ID"]
)

extracted = set(
    mapping["Gene_ID"]
)

hit = set(
    best["Gene_ID"]
)

print(
    "\nSelected but protein missing:"
)

for x in sorted(
    selected - extracted
):
    print(x)

print(
    "\nProtein extracted but no BLAST hit:"
)

for x in sorted(
    extracted - hit
):
    print(x)

print(
    "\nCcac_g2721 status:"
)

print(
    "Selected:",
    "Ccac_g2721" in selected
)

print(
    "Protein extracted:",
    "Ccac_g2721" in extracted
)

print(
    "BLAST hit:",
    "Ccac_g2721" in hit
)

PY

###############################################################################
# 16. OUTPUT LIST
###############################################################################

echo
echo "=============================================================="
echo "FINAL OUTPUT FILES"
echo "=============================================================="

echo
echo "Gene lists:"
ls -lh 01_gene_lists/

echo
echo "Proteins:"
ls -lh 02_proteins/

echo
echo "BLAST:"
ls -lh 03_blast/

echo
echo "Summary:"
ls -lh 04_summary/

echo
echo "Updated candidate tables:"
ls -lh 05_tables/

echo
echo "QC:"
ls -lh 06_QC/

echo
echo "=============================================================="
echo "FINISHED"
echo "Date: $(date)"
echo "=============================================================="
