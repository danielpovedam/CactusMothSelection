#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# PATHS
# ============================================================

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

GENO_FILE <- file.path(
  BASE,
  "../baypass_preparation/baypass/all_snps.baypass.geno"
)

META_FILE <- file.path(
  BASE,
  "../baypass_preparation/metadata/all_snps.baypass.snp_metadata.tsv"
)

OUTDIR <- file.path(
  BASE,
  "sensitivity_C2_geographic"
)

INPUT_DIR <- file.path(
  OUTDIR,
  "input"
)

QC_DIR <- file.path(
  OUTDIR,
  "qc"
)

dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# ORIGINAL BAYPASS POPULATION ORDER
# ============================================================

FULL_POP_ORDER <- c(
  "BAM", "BAP", "BAS", "CAQ", "CAS",
  "CBA", "CBF", "CBR", "CBS", "CBV",
  "CHP", "ERB", "ERC", "ERD", "ERF",
  "ERL", "ERN", "FOG", "JUJ", "JUS",
  "LPA", "LPS", "LPV", "RNR",
  "SEB", "SEL", "SEP", "TUN"
)

# ============================================================
# GEOGRAPHICALLY RESTRICTED POPULATIONS
#
# Keep FULL BayPass order.
# ============================================================

KEEP_POPS <- c(
  "BAM",
  "BAP",
  "BAS",
  "CBA",
  "CBF",
  "CBR",
  "CBV",
  "ERB",
  "ERC",
  "ERD",
  "ERF",
  "ERL",
  "ERN",
  "SEB",
  "SEL",
  "SEP"
)

HOST_C2 <- c(
  BAM =  1,
  BAP = -1,
  BAS =  1,
  CBA = -1,
  CBF = -1,
  CBR =  1,
  CBV =  1,
  ERB = -1,
  ERC = -1,
  ERD =  1,
  ERF =  1,
  ERL = -1,
  ERN = -1,
  SEB =  1,
  SEL = -1,
  SEP =  1
)

EXPECTED_SNPS <- 132283L
EXPECTED_FULL_POPS <- 28L
EXPECTED_KEEP_POPS <- 16L

# ============================================================
# BASIC CHECKS
# ============================================================

if (!file.exists(GENO_FILE)) {
  stop("Missing genotype file:\n", GENO_FILE)
}

if (!file.exists(META_FILE)) {
  stop("Missing metadata file:\n", META_FILE)
}

if (length(FULL_POP_ORDER) != EXPECTED_FULL_POPS) {
  stop("FULL_POP_ORDER does not contain 28 populations.")
}

if (length(KEEP_POPS) != EXPECTED_KEEP_POPS) {
  stop("KEEP_POPS does not contain 16 populations.")
}

if (!all(KEEP_POPS %in% FULL_POP_ORDER)) {
  stop(
    "These restricted populations are absent from FULL_POP_ORDER:\n",
    paste(
      setdiff(KEEP_POPS, FULL_POP_ORDER),
      collapse = ", "
    )
  )
}

if (!identical(
  KEEP_POPS,
  FULL_POP_ORDER[FULL_POP_ORDER %in% KEEP_POPS]
)) {
  stop(
    "KEEP_POPS is not in the same order as the original BayPass genotype file."
  )
}

if (!identical(
  names(HOST_C2),
  KEEP_POPS
)) {
  stop("HOST_C2 names/order do not match KEEP_POPS.")
}

if (sum(HOST_C2 == 1) != 8L ||
    sum(HOST_C2 == -1) != 8L) {
  stop("Host contrast is not balanced 8 vs 8.")
}

# ============================================================
# READ GENOTYPE MATRIX
# ============================================================

cat("[41A] Reading full BayPass genotype matrix...\n")

geno <- fread(
  GENO_FILE,
  header = FALSE,
  showProgress = TRUE
)

if (nrow(geno) != EXPECTED_SNPS) {
  stop(
    "Expected ",
    EXPECTED_SNPS,
    " SNP rows but found ",
    nrow(geno)
  )
}

if (ncol(geno) != 2L * EXPECTED_FULL_POPS) {
  stop(
    "Expected 56 columns for 28 populations but found ",
    ncol(geno)
  )
}

# ============================================================
# IDENTIFY COLUMN PAIRS
# ============================================================

pop_index <- match(
  KEEP_POPS,
  FULL_POP_ORDER
)

keep_columns <- as.vector(
  rbind(
    2L * pop_index - 1L,
    2L * pop_index
  )
)

if (length(keep_columns) != 32L) {
  stop("Did not identify 32 allele-count columns.")
}

restricted <- geno[
  ,
  ..keep_columns
]

if (ncol(restricted) != 32L) {
  stop("Restricted genotype matrix does not contain 32 columns.")
}

# ============================================================
# QC: TOTAL ALLELE COUNTS AND POLYMORPHISM
#
# BayPass representation:
# columns 1/2 = allele 1/2 counts in population 1,
# columns 3/4 = allele 1/2 counts in population 2, etc.
# ============================================================

allele1_cols <- seq(
  1L,
  ncol(restricted),
  by = 2L
)

allele2_cols <- seq(
  2L,
  ncol(restricted),
  by = 2L
)

A1_total <- rowSums(
  as.matrix(
    restricted[, ..allele1_cols]
  )
)

A2_total <- rowSums(
  as.matrix(
    restricted[, ..allele2_cols]
  )
)

total_called <- A1_total + A2_total

polymorphic <- (
  A1_total > 0 &
  A2_total > 0
)

fixed_allele1 <- (
  A1_total > 0 &
  A2_total == 0
)

fixed_allele2 <- (
  A1_total == 0 &
  A2_total > 0
)

no_data <- (
  total_called == 0
)

n_poly <- sum(polymorphic)
n_fixed <- sum(fixed_allele1 | fixed_allele2)
n_no_data <- sum(no_data)

cat(
  "[41A] SNPs polymorphic in restricted subset: ",
  n_poly,
  "\n",
  sep = ""
)

cat(
  "[41A] SNPs fixed in restricted subset: ",
  n_fixed,
  "\n",
  sep = ""
)

cat(
  "[41A] SNPs with zero allele counts: ",
  n_no_data,
  "\n",
  sep = ""
)

# ============================================================
# SAVE FULL SAME-ROW RESTRICTED MATRIX
#
# This preserves exact correspondence with the 132,283-SNP
# primary analysis.
# ============================================================

OUT_GENO_ALL <- file.path(
  INPUT_DIR,
  "all_snps.C2_geographic_16pop.same132283.geno"
)

fwrite(
  restricted,
  OUT_GENO_ALL,
  sep = " ",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# ALSO SAVE POLYMORPHIC-ONLY VERSION
#
# Do not use this automatically yet.
# It is prepared so we can decide after QC.
# ============================================================

OUT_GENO_POLY <- file.path(
  INPUT_DIR,
  "all_snps.C2_geographic_16pop.polymorphic_only.geno"
)

fwrite(
  restricted[
    polymorphic
  ],
  OUT_GENO_POLY,
  sep = " ",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# SNP METADATA / ROW MAP
# ============================================================

meta <- fread(
  META_FILE,
  showProgress = FALSE
)

if (nrow(meta) != EXPECTED_SNPS) {
  stop(
    "Metadata does not contain ",
    EXPECTED_SNPS,
    " rows."
  )
}

meta[
  ,
  `:=`(
    Restricted_A1_total = A1_total,
    Restricted_A2_total = A2_total,
    Restricted_total_called = total_called,
    Restricted_polymorphic = polymorphic
  )
]

fwrite(
  meta,
  file.path(
    QC_DIR,
    "restricted_C2_allSNP_row_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

fwrite(
  meta[
    polymorphic
  ],
  file.path(
    INPUT_DIR,
    "restricted_C2_polymorphic_SNP_metadata.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

# ============================================================
# HOST CONTRAST FILE
# ============================================================

contrast <- data.table(
  Host_C2 = as.numeric(
    HOST_C2[KEEP_POPS]
  )
)

fwrite(
  contrast,
  file.path(
    INPUT_DIR,
    "Host_C2_geographic_16pop.txt"
  ),
  sep = " ",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# POPULATION ORDER FILE
# ============================================================

pop_table <- data.table(
  BayPass_order = seq_along(
    KEEP_POPS
  ),
  Population = KEEP_POPS,
  Host_C2 = as.numeric(
    HOST_C2[KEEP_POPS]
  ),
  Host_group = ifelse(
    HOST_C2[KEEP_POPS] == 1,
    "O_ficus_indica",
    "native_Opuntia"
  )
)

fwrite(
  pop_table,
  file.path(
    INPUT_DIR,
    "population_order_C2_geographic_16pop.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

# ============================================================
# QC SUMMARY
# ============================================================

qc <- data.table(
  Metric = c(
    "Original SNP rows",
    "Original populations",
    "Restricted populations",
    "O_ficus_indica populations",
    "Native Opuntia populations",
    "Restricted genotype columns",
    "Polymorphic SNPs in restricted subset",
    "Fixed SNPs in restricted subset",
    "SNPs with zero counts in restricted subset"
  ),
  Value = c(
    EXPECTED_SNPS,
    EXPECTED_FULL_POPS,
    EXPECTED_KEEP_POPS,
    sum(HOST_C2 == 1),
    sum(HOST_C2 == -1),
    ncol(restricted),
    n_poly,
    n_fixed,
    n_no_data
  )
)

fwrite(
  qc,
  file.path(
    QC_DIR,
    "QC_C2_geographic_subset_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

cat(
  "\n[41A] OK\n",
  "  Restricted populations: 16 (8 exotic / 8 native)\n",
  "  Full-row genotype file: ", OUT_GENO_ALL, "\n",
  "  Polymorphic-only genotype file: ", OUT_GENO_POLY, "\n",
  "  Host contrast: ",
  file.path(INPUT_DIR, "Host_C2_geographic_16pop.txt"),
  "\n",
  "  QC: ",
  file.path(QC_DIR, "QC_C2_geographic_subset_summary.tsv"),
  "\n",
  sep = ""
)
