#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =========================================================
# SETTINGS
# =========================================================

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

RUNBASE <- file.path(
  BASE,
  "runs_10x",
  "all_snps"
)

META <- paste0(
  "/home/dpoveda/micetral/Population_snps/",
  "baypass_preparation/metadata/",
  "all_snps.baypass.snp_metadata.tsv"
)

OUTDIR <- file.path(
  BASE,
  "selection_scan",
  "GEA_BF_median3seeds_allSNP"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

OUTFILE <- file.path(
  OUTDIR,
  "merged_sorted_BF_median3seeds.tsv"
)

SEEDS <- c(
  "12345",
  "67890",
  "13579"
)

COV_NAMES <- c(
  "TempPC1",
  "PrecipPC1",
  "AridityIndex",
  "Elevation",
  "Latitude",
  "Longitude",
  "BIO1",
  "BIO12"
)

# =========================================================
# HELPERS
# =========================================================

read_seed <- function(seed) {

  f <- file.path(
    RUNBASE,
    paste0("seed_", seed),
    paste0(
      "BP_all_snps_seed",
      seed,
      "_summary_betai_reg.out"
    )
  )

  stopifnot(file.exists(f))

  x <- fread(f)

  stopifnot(
    all(
      c(
        "COVARIABLE",
        "MRK",
        "BF(dB)",
        "Beta_is",
        "SD_Beta_is",
        "eBPis"
      ) %in% names(x)
    )
  )

  x <- x[, .(
    COVARIABLE = as.integer(COVARIABLE),
    MRK = as.integer(MRK),
    BFdb = as.numeric(get("BF(dB)")),
    Beta_is = as.numeric(Beta_is),
    eBPis = as.numeric(eBPis)
  )]

  setnames(
    x,
    c("BFdb", "Beta_is", "eBPis"),
    c(
      paste0("BFdb_seed", seed),
      paste0("Beta_seed", seed),
      paste0("eBPis_seed", seed)
    )
  )

  setkey(
    x,
    COVARIABLE,
    MRK
  )

  x
}

# =========================================================
# READ THREE BAYPASS RUNS
# =========================================================

cat("Reading BayPass runs...\n")

x1 <- read_seed("12345")
x2 <- read_seed("67890")
x3 <- read_seed("13579")

stopifnot(
  nrow(x1) == nrow(x2),
  nrow(x1) == nrow(x3)
)

# =========================================================
# MERGE USING KEYS
# =========================================================

x <- merge(
  x1,
  x2,
  by = c("COVARIABLE", "MRK"),
  all = TRUE
)

x <- merge(
  x,
  x3,
  by = c("COVARIABLE", "MRK"),
  all = TRUE
)

stopifnot(
  nrow(x) == nrow(x1)
)

# =========================================================
# MEDIAN BF ACROSS THREE RUNS
# =========================================================

x[, BFdb := apply(
  .SD,
  1,
  median,
  na.rm = TRUE
),
.SDcols = c(
  "BFdb_seed12345",
  "BFdb_seed67890",
  "BFdb_seed13579"
)]

# Keep summaries for QC only.

x[, BFdb_min := pmin(
  BFdb_seed12345,
  BFdb_seed67890,
  BFdb_seed13579,
  na.rm = TRUE
)]

x[, BFdb_max := pmax(
  BFdb_seed12345,
  BFdb_seed67890,
  BFdb_seed13579,
  na.rm = TRUE
)]

# =========================================================
# ADD COVARIATE NAMES
# =========================================================

stopifnot(
  all(
    sort(unique(x$COVARIABLE)) ==
      seq_along(COV_NAMES)
  )
)

x[, Covariate := COV_NAMES[COVARIABLE]]

# =========================================================
# ADD SNP METADATA
# =========================================================

meta <- fread(META)

stopifnot(
  all(
    c(
      "BayPass_row",
      "CHROM",
      "POS",
      "SNP",
      "REF",
      "ALT"
    ) %in% names(meta)
  )
)

meta <- meta[, .(
  MRK = as.integer(BayPass_row),
  chr = as.character(CHROM),
  pos = as.integer(POS),
  SNP,
  REF,
  ALT
)]

stopifnot(!anyDuplicated(meta$MRK))

x <- merge(
  x,
  meta,
  by = "MRK",
  all.x = TRUE
)

stopifnot(
  !anyNA(x$chr),
  !anyNA(x$pos)
)

# =========================================================
# ADD M_P
# =========================================================

PI_FILE <- file.path(
  RUNBASE,
  "seed_12345",
  "BP_all_snps_seed12345_summary_pi_xtx.out"
)

stopifnot(file.exists(PI_FILE))

pi <- fread(PI_FILE)

stopifnot(
  all(
    c("MRK", "M_P") %in% names(pi)
  )
)

pi <- pi[, .(
  MRK = as.integer(MRK),
  M_P = as.numeric(M_P)
)]

stopifnot(!anyDuplicated(pi$MRK))

x <- merge(
  x,
  pi,
  by = "MRK",
  all.x = TRUE
)

stopifnot(!anyNA(x$M_P))

# =========================================================
# DATASET LABEL
# =========================================================

x[, subset := "all_snps"]

# =========================================================
# ORDER EXACTLY LIKE MOUSE INPUT
# =========================================================

setorder(
  x,
  COVARIABLE,
  chr,
  pos,
  MRK
)

# =========================================================
# QC
# =========================================================

cat("\n=========================================\n")
cat("THREE-SEED BF MEDIAN QC\n")
cat("=========================================\n\n")

for (cov in sort(unique(x$COVARIABLE))) {

  z <- x[COVARIABLE == cov]

  cat(
    "\n",
    cov,
    " ",
    COV_NAMES[cov],
    "\n",
    sep = ""
  )

  cmat <- cor(
    z[, .(
      BFdb_seed12345,
      BFdb_seed67890,
      BFdb_seed13579
    )],
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  print(cmat)

  cat(
    "Median BF max:",
    max(z$BFdb),
    "\n"
  )
}

# =========================================================
# WRITE
# =========================================================

fwrite(
  x,
  OUTFILE,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n=========================================\n")
cat("BF MEDIAN FILE COMPLETED\n")
cat("=========================================\n")

cat(
  "Rows:",
  nrow(x),
  "\n"
)

cat(
  "Output:\n",
  OUTFILE,
  "\n",
  sep = ""
)
