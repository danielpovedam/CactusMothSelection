#!/bin/bash
#SBATCH --job-name=Ccac_candidate_BLAST
#SBATCH --output=Ccac_candidate_BLAST_%j.out
#SBATCH --error=Ccac_candidate_BLAST_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail

###############################################################################
# Cactoblastis cactorum
# Refine annotations of incompletely annotated XtX, C2 and GEA candidate genes
#
# Strategy:
#   1. Read final XtX, C2 and GEA candidate tables
#   2. Select genes with:
#        - blank gene symbol
#        - gene symbol equal to Gene ID
#        - hypothetical/uncharacterized description
#   3. Extract longest predicted protein isoform per gene
#   4. BLASTP against RefSeq protein and nr
#   5. Select best-supported hit
#   6. Add BLAST evidence columns back to all candidate tables
###############################################################################

THREADS=${SLURM_CPUS_PER_TASK:-8}

###############################################################################
# EDIT ONLY THESE PATHS
###############################################################################

WORKDIR="$HOME/micetral/Population_snps/baypass/selection_scan/functional_annotation"

PROTEINS="$HOME/micetral/Population_snps/Ccactorum_predicted_proteins.fa"

# FINAL candidate tables.
# Replace these paths only if your actual filenames differ.
XTX_TABLE="$HOME/micetral/Population_snps/baypass/selection_scan/FINAL_gt1Mb/Table_XtX_candidate_regions_allSNP_gt1Mb_FINAL.tsv"

C2_TABLE="$HOME/micetral/Population_snps/baypass/selection_scan/FINAL_gt1Mb/Table_C2_Host_candidate_regions_allSNP_gt1Mb_FINAL.tsv"

GEA_TABLE="$HOME/micetral/Population_snps/baypass/selection_scan/Table_GEA_candidate_loci_FINAL.tsv"

###############################################################################
# BLAST databases
#
# IMPORTANT:
# Put here the SAME nr and RefSeq BLAST database paths that you used for Masc.
###############################################################################

NR_DB="/PATH/TO/BLAST_DB/nr"
REFSEQ_DB="/PATH/TO/BLAST_DB/refseq_protein"

###############################################################################

mkdir -p "${WORKDIR}"/{01_gene_lists,02_proteins,03_blast,04_summary,05_tables}
cd "${WORKDIR}"

echo "============================================================"
echo "Candidate annotation refinement"
echo "Started: $(date)"
echo "============================================================"

###############################################################################
# MODULES
###############################################################################

module purge

# Adapt module names if your cluster uses slightly different names.
module load bioinfo/BLAST/2.15.0
module load bioinfo/samtools

blastp -version
samtools --version | head -1

###############################################################################
# CHECK INPUTS
###############################################################################

for f in "$PROTEINS" "$XTX_TABLE" "$C2_TABLE" "$GEA_TABLE"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing input file:"
        echo "$f"
        exit 1
    fi
done

###############################################################################
# STEP 1
# Identify genes requiring annotation refinement
###############################################################################

python3 <<'PY'
import pandas as pd
import re
from pathlib import Path

tables = {
    "XtX": Path("""'"$XTX_TABLE"'""".strip("'")),
    "C2":  Path("""'"$C2_TABLE"'""".strip("'")),
    "GEA": Path("""'"$GEA_TABLE"'""".strip("'"))
}

outdir = Path("01_gene_lists")
outdir.mkdir(exist_ok=True)

all_rows = []

def find_column(df, possibilities):
    norm = {
        str(c).strip().lower().replace(" ", "_"): c
        for c in df.columns
    }
    for x in possibilities:
        key = x.lower().replace(" ", "_")
        if key in norm:
            return norm[key]
    return None

for analysis, path in tables.items():

    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    gene_col = find_column(
        df,
        ["Gene ID", "Gene_ID", "gene_id", "GeneID"]
    )

    symbol_col = find_column(
        df,
        ["Gene symbol", "Gene_symbol", "gene_symbol", "Symbol"]
    )

    desc_col = find_column(
        df,
        ["Description", "description"]
    )

    if gene_col is None:
        raise RuntimeError(
            f"Could not identify Gene ID column in {path}\n"
            f"Columns = {list(df.columns)}"
        )

    for _, r in df.iterrows():

        gene = str(r.get(gene_col, "")).strip()
        symbol = str(r.get(symbol_col, "")).strip() if symbol_col else ""
        desc = str(r.get(desc_col, "")).strip() if desc_col else ""

        if gene in {"", "-", "NA", "na", "None"}:
            continue

        # Standardize possible markdown escapes
        gene = gene.replace("\\_", "_")
        symbol_clean = symbol.replace("*", "").replace("\\_", "_").strip()

        desc_lower = desc.lower()

        incomplete = (
            symbol_clean == "" or
            symbol_clean == "-" or
            symbol_clean == gene or
            "hypothetical" in desc_lower or
            "uncharacterized" in desc_lower or
            "unknown protein" in desc_lower
        )

        if incomplete:
            all_rows.append({
                "Analysis": analysis,
                "Gene_ID": gene,
                "Current_symbol": symbol_clean,
                "Current_description": desc
            })

result = pd.DataFrame(all_rows)

if result.empty:
    raise RuntimeError("No incompletely annotated genes were detected.")

result = result.drop_duplicates()

result.to_csv(
    outdir / "candidate_genes_requiring_BLAST_by_analysis.tsv",
    sep="\t",
    index=False
)

genes = sorted(result["Gene_ID"].unique())

with open(outdir / "candidate_genes_requiring_BLAST.txt", "w") as f:
    for g in genes:
        f.write(g + "\n")

print("\nGenes requiring BLAST refinement:", len(genes))
print("\nBy analysis:")
print(result.groupby("Analysis")["Gene_ID"].nunique())

print("\nGenes:")
for x in genes:
    print(x)

# Explicit check for Figure 4D candidate
if "Ccac_g2721" in genes:
    print("\nSUCCESS: Ccac_g2721 is included for BLAST refinement.")
else:
    print(
        "\nWARNING: Ccac_g2721 was not automatically selected. "
        "Adding it explicitly because it is the Figure 4D candidate."
    )
    genes.append("Ccac_g2721")
    genes = sorted(set(genes))

    with open(outdir / "candidate_genes_requiring_BLAST.txt", "w") as f:
        for g in genes:
            f.write(g + "\n")
PY

###############################################################################
# STEP 2
# Find all protein isoforms belonging to each gene and retain longest isoform
###############################################################################

python3 <<'PY'
from pathlib import Path
import re

protein_fasta = Path("""'"$PROTEINS"'""".strip("'"))
gene_file = Path("01_gene_lists/candidate_genes_requiring_BLAST.txt")

genes = {
    x.strip()
    for x in gene_file.read_text().splitlines()
    if x.strip()
}

seqs = {}

header = None
seq = []

def save_record(header, seq):
    if header is None:
        return

    name = header.split()[0]

    # Supports:
    # Ccac_g2721
    # Ccac_g2721.t1
    # Ccac_g2721.t2
    # etc.
    for gene in genes:
        if (
            name == gene or
            name.startswith(gene + ".") or
            name.startswith(gene + "-")
        ):
            protein = "".join(seq).replace("*", "")
            seqs.setdefault(gene, []).append((name, protein))
            break

with open(protein_fasta) as f:
    for line in f:
        line = line.rstrip()

        if line.startswith(">"):
            save_record(header, seq)
            header = line[1:]
            seq = []
        else:
            seq.append(line)

save_record(header, seq)

out = Path("02_proteins/candidate_longest_isoform.fa")
mapping = Path("02_proteins/candidate_longest_isoform_mapping.tsv")

missing = []

with open(out, "w") as fo, open(mapping, "w") as fm:

    fm.write("Gene_ID\tTranscript_ID\tProtein_length_aa\n")

    for gene in sorted(genes):

        if gene not in seqs:
            missing.append(gene)
            continue

        # Longest protein isoform
        transcript, protein = max(
            seqs[gene],
            key=lambda x: len(x[1])
        )

        fo.write(f">{gene}|{transcript}\n")

        for i in range(0, len(protein), 60):
            fo.write(protein[i:i+60] + "\n")

        fm.write(
            f"{gene}\t{transcript}\t{len(protein)}\n"
        )

print("Genes requested:", len(genes))
print("Genes with protein:", len(seqs))
print("Missing proteins:", len(missing))

if missing:
    print("\nMissing:")
    for x in missing:
        print(x)

# Key check
if "Ccac_g2721" in seqs:
    print("\nCcac_g2721 protein successfully extracted.")
PY

echo
echo "Protein mapping:"
column -t -s $'\t' 02_proteins/candidate_longest_isoform_mapping.tsv | head -50

###############################################################################
# STEP 3
# BLASTP against RefSeq protein
###############################################################################

echo
echo "Running BLASTP against RefSeq protein..."
echo

blastp \
    -query 02_proteins/candidate_longest_isoform.fa \
    -db "${REFSEQ_DB}" \
    -evalue 1e-3 \
    -max_target_seqs 20 \
    -num_threads "${THREADS}" \
    -seg yes \
    -soft_masking true \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs stitle" \
    > 03_blast/refseq_protein.BLASTP.tsv

###############################################################################
# STEP 4
# BLASTP against nr
###############################################################################

echo
echo "Running BLASTP against nr..."
echo

blastp \
    -query 02_proteins/candidate_longest_isoform.fa \
    -db "${NR_DB}" \
    -evalue 1e-3 \
    -max_target_seqs 20 \
    -num_threads "${THREADS}" \
    -seg yes \
    -soft_masking true \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs stitle" \
    > 03_blast/nr.BLASTP.tsv

###############################################################################
# STEP 5
# Parse all BLAST results and choose best-supported hit
#
# Ranking:
#   lower E-value
#   then higher bit score
#   then higher query coverage
#   then higher identity
###############################################################################

python3 <<'PY'
import pandas as pd
from pathlib import Path
import numpy as np

cols = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "qlen",
    "sstart",
    "send",
    "slen",
    "evalue",
    "bitscore",
    "qcovs",
    "stitle"
]

files = {
    "RefSeq": Path("03_blast/refseq_protein.BLASTP.tsv"),
    "nr": Path("03_blast/nr.BLASTP.tsv")
}

dfs = []

for db, f in files.items():

    if not f.exists() or f.stat().st_size == 0:
        continue

    x = pd.read_csv(
        f,
        sep="\t",
        names=cols,
        header=None
    )

    x["BLAST_database"] = db

    # qseqid = GeneID|Transcript
    x["Gene_ID"] = (
        x["qseqid"]
        .astype(str)
        .str.split("|", regex=False)
        .str[0]
    )

    dfs.append(x)

if not dfs:
    raise RuntimeError("No BLAST results found.")

blast = pd.concat(dfs, ignore_index=True)

for col in ["pident", "qcovs", "evalue", "bitscore"]:
    blast[col] = pd.to_numeric(blast[col], errors="coerce")

###############################################################################
# Save all hits
###############################################################################

blast.to_csv(
    "04_summary/ALL_candidate_BLAST_hits.tsv",
    sep="\t",
    index=False
)

###############################################################################
# Choose best hit per gene
###############################################################################

best = (
    blast
    .sort_values(
        ["Gene_ID", "evalue", "bitscore", "qcovs", "pident"],
        ascending=[True, True, False, False, False]
    )
    .groupby("Gene_ID", as_index=False)
    .first()
)

best_out = best[[
    "Gene_ID",
    "BLAST_database",
    "sseqid",
    "stitle",
    "pident",
    "qcovs",
    "evalue",
    "bitscore"
]].copy()

best_out.columns = [
    "Gene_ID",
    "BLAST_database",
    "BLAST_accession",
    "BLAST_hit",
    "BLAST_identity_pct",
    "BLAST_query_coverage_pct",
    "BLAST_Evalue",
    "BLAST_bitscore"
]

best_out.to_csv(
    "04_summary/BEST_BLAST_hit_per_candidate_gene.tsv",
    sep="\t",
    index=False
)

print("\nBest BLAST hit per gene:")
print(best_out.to_string(index=False))

print("\nCcac_g2721:")
print(
    best_out[
        best_out["Gene_ID"] == "Ccac_g2721"
    ].to_string(index=False)
)
PY

###############################################################################
# STEP 6
# Add BLAST columns to XtX, C2 and GEA tables
###############################################################################

python3 <<'PY'
import pandas as pd
from pathlib import Path

tables = {
    "XtX": Path("""'"$XTX_TABLE"'""".strip("'")),
    "C2":  Path("""'"$C2_TABLE"'""".strip("'")),
    "GEA": Path("""'"$GEA_TABLE"'""".strip("'"))
}

outdir = Path("05_tables")
outdir.mkdir(exist_ok=True)

blast = pd.read_csv(
    "04_summary/BEST_BLAST_hit_per_candidate_gene.tsv",
    sep="\t",
    dtype=str,
    keep_default_na=False
)

def find_column(df, possibilities):
    norm = {
        str(c).strip().lower().replace(" ", "_"): c
        for c in df.columns
    }

    for x in possibilities:
        key = x.lower().replace(" ", "_")
        if key in norm:
            return norm[key]

    return None

blast_cols = [
    "BLAST_database",
    "BLAST_accession",
    "BLAST_hit",
    "BLAST_identity_pct",
    "BLAST_query_coverage_pct",
    "BLAST_Evalue",
    "BLAST_bitscore"
]

for analysis, path in tables.items():

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )

    gene_col = find_column(
        df,
        ["Gene ID", "Gene_ID", "gene_id", "GeneID"]
    )

    if gene_col is None:
        raise RuntimeError(
            f"Could not identify Gene ID in {path}"
        )

    df["_Gene_ID_merge"] = (
        df[gene_col]
        .astype(str)
        .str.replace("\\_", "_", regex=False)
        .str.strip()
    )

    b = blast.rename(
        columns={"Gene_ID": "_Gene_ID_merge"}
    )

    df = df.merge(
        b,
        on="_Gene_ID_merge",
        how="left"
    )

    df.drop(
        columns=["_Gene_ID_merge"],
        inplace=True
    )

    # Empty instead of NA for genes not BLASTed
    for c in blast_cols:
        if c in df.columns:
            df[c] = df[c].fillna("")

    outfile = outdir / f"{analysis}_candidate_table_with_BLAST.tsv"

    df.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    print(
        analysis,
        "rows =", len(df),
        "BLAST-annotated rows =",
        (df["BLAST_Evalue"] != "").sum()
    )

print("\nFinished updated tables.")
PY

###############################################################################
# STEP 7
# Specific Figure 4D report
###############################################################################

echo
echo "============================================================"
echo "FIGURE 4D: ptg000019 TemperaturePC1 candidate"
echo "============================================================"

python3 <<'PY'
import pandas as pd

gea = pd.read_csv(
    "05_tables/GEA_candidate_table_with_BLAST.tsv",
    sep="\t",
    dtype=str,
    keep_default_na=False
)

# Flexible column matching
def col(names):
    lookup = {
        x.lower().replace(" ", "_"): x
        for x in gea.columns
    }
    for n in names:
        k = n.lower().replace(" ", "_")
        if k in lookup:
            return lookup[k]
    return None

cov = col(["Covariate"])
scaf = col(["Scaffold"])
gene = col(["Gene ID", "Gene_ID"])

subset = gea[
    (gea[cov] == "TempPC1") &
    (gea[scaf] == "ptg000019")
]

print(subset.to_string(index=False))

print("\nRows involving Ccac_g2721:")
print(
    gea[
        gea[gene]
        .str.replace("\\_", "_", regex=False)
        .eq("Ccac_g2721")
    ].to_string(index=False)
)
PY

echo
echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"

