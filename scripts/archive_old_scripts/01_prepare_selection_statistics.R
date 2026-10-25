#!/usr/bin/env Rscript

############################################################
## PREPARE BAYPASS SELECTION STATISTICS
## Cactoblastis cactorum
##
## XtX:
##   representative converged seed = 12345
##
## C2:
##   representative converged seed = 12345
##
## GEA:
##   mean BF(dB) across seeds:
##   12345, 67890, 13579
##
## Datasets:
##   random_snp
##   all_snps
############################################################


suppressPackageStartupMessages({
  library(data.table)
})


############################################################
## 0. SETTINGS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

RUNBASE <- file.path(
  BASE,
  "runs_10x"
)

METABASE <- paste0(
  "/home/dpoveda/micetral/Population_snps/",
  "baypass_preparation/metadata"
)

OUTBASE <- file.path(
  BASE,
  "selection_scan",
  "prepared_statistics"
)

dir.create(
  OUTBASE,
  recursive = TRUE,
  showWarnings = FALSE
)


DATASETS <- c(
  "random_snp",
  "all_snps"
)


SEEDS <- c(
  "12345",
  "67890",
  "13579"
)


REP_SEED <- "12345"


############################################################
## 1. COVARIATE MAP
############################################################

COV_MAP <- data.table(
  COVARIABLE = 1:8,
  Covariate = c(
    "TempPC1",
    "PrecipPC1",
    "AridityIndex",
    "Elevation",
    "Latitude",
    "Longitude",
    "BIO1",
    "BIO12"
  )
)


############################################################
## 2. FUNCTIONS
############################################################

get_prefix <- function(dataset, seed) {

  file.path(
    RUNBASE,
    dataset,
    paste0("seed_", seed),
    paste0(
      "BP_",
      dataset,
      "_seed",
      seed
    )
  )
}


read_metadata <- function(dataset) {

  file <- file.path(
    METABASE,
    paste0(
      dataset,
      ".baypass.snp_metadata.tsv"
    )
  )

  if (!file.exists(file)) {
    stop(
      "Missing metadata file: ",
      file
    )
  }


  dt <- fread(
    file,
    showProgress = FALSE
  )


  required <- c(
    "BayPass_row",
    "CHROM",
    "SNP",
    "POS",
    "ALT",
    "REF"
  )


  if (!all(required %in% names(dt))) {
    stop(
      "Metadata missing required columns: ",
      file
    )
  }


  dt[, BayPass_row := as.integer(BayPass_row)]

  dt[, POS := as.integer(POS)]


  if (anyDuplicated(dt$BayPass_row)) {
    stop(
      "Duplicated BayPass_row values in metadata: ",
      dataset
    )
  }


  return(dt)
}


read_xtx <- function(dataset, seed) {

  file <- paste0(
    get_prefix(
      dataset,
      seed
    ),
    "_summary_pi_xtx.out"
  )


  if (!file.exists(file)) {
    stop(
      "Missing XtX file: ",
      file
    )
  }


  dt <- fread(
    file,
    showProgress = FALSE
  )


  required <- c(
    "MRK",
    "M_P",
    "SD_P",
    "M_XtX",
    "SD_XtX",
    "XtXst",
    "log10(1/pval)"
  )


  if (!all(required %in% names(dt))) {
    stop(
      "XtX file missing required columns: ",
      file,
      "\nColumns found: ",
      paste(
        names(dt),
        collapse = ", "
      )
    )
  }


  dt[, MRK := as.integer(MRK)]


  return(dt)
}


read_c2 <- function(dataset, seed) {

  file <- paste0(
    get_prefix(
      dataset,
      seed
    ),
    "_summary_contrast.out"
  )


  if (!file.exists(file)) {
    stop(
      "Missing contrast file: ",
      file
    )
  }


  dt <- fread(
    file,
    showProgress = FALSE
  )


  required <- c(
    "CONTRAST",
    "MRK",
    "C2",
    "log10(1/pval)"
  )


  if (!all(required %in% names(dt))) {
    stop(
      "C2 file missing required columns: ",
      file,
      "\nColumns found: ",
      paste(
        names(dt),
        collapse = ", "
      )
    )
  }


  dt[, CONTRAST := as.integer(CONTRAST)]

  dt[, MRK := as.integer(MRK)]


  return(dt)
}


read_bf <- function(dataset, seed) {

  file <- paste0(
    get_prefix(
      dataset,
      seed
    ),
    "_summary_betai_reg.out"
  )


  if (!file.exists(file)) {
    stop(
      "Missing BF file: ",
      file
    )
  }


  dt <- fread(
    file,
    showProgress = FALSE
  )


  required <- c(
    "COVARIABLE",
    "MRK",
    "BF(dB)"
  )


  if (!all(required %in% names(dt))) {
    stop(
      "BF file missing required columns: ",
      file,
      "\nColumns found: ",
      paste(
        names(dt),
        collapse = ", "
      )
    )
  }


  dt[, COVARIABLE := as.integer(COVARIABLE)]

  dt[, MRK := as.integer(MRK)]


  dt <- dt[
    ,
    .(
      MRK,
      COVARIABLE,
      BFdb = as.numeric(`BF(dB)`)
    )
  ]


  return(dt)
}


############################################################
## 3. GLOBAL BF CORRELATION CONTAINER
############################################################

BF_COR_RESULTS <- list()

counter <- 1


############################################################
## 4. PROCESS DATASETS
############################################################

for (dataset in DATASETS) {

  cat(
    "\n========================================\n"
  )

  cat(
    "DATASET: ",
    dataset,
    "\n",
    sep = ""
  )

  cat(
    "========================================\n"
  )


  out_ds <- file.path(
    OUTBASE,
    dataset
  )


  dir.create(
    out_ds,
    recursive = TRUE,
    showWarnings = FALSE
  )


  ##########################################################
  ## 4.1 Metadata
  ##########################################################

  META <- read_metadata(
    dataset
  )


  cat(
    "Metadata SNPs: ",
    nrow(META),
    "\n",
    sep = ""
  )


  ##########################################################
  ## 4.2 XtX from representative seed
  ##########################################################

  XTX_RAW <- read_xtx(
    dataset,
    REP_SEED
  )


  if (nrow(XTX_RAW) != nrow(META)) {

    stop(
      "XtX and metadata SNP counts differ for ",
      dataset,
      ": ",
      nrow(XTX_RAW),
      " vs ",
      nrow(META)
    )
  }


  XTX <- merge(
    META,
    XTX_RAW,
    by.x = "BayPass_row",
    by.y = "MRK",
    all.x = TRUE,
    sort = FALSE
  )


  if (anyNA(XTX$M_XtX)) {

    stop(
      "Missing XtX values after metadata merge: ",
      dataset
    )
  }


  setorder(
    XTX,
    CHROM,
    POS
  )


  XTX_OUT <- XTX[
    ,
    .(
      BayPass_row,
      CHROM,
      SNP,
      POS,
      REF,
      ALT,
      M_P,
      SD_P,
      M_XtX,
      SD_XtX,
      XtXst,
      log10_p = `log10(1/pval)`
    )
  ]


  fwrite(
    XTX_OUT,
    file.path(
      out_ds,
      paste0(
        "XtX_seed",
        REP_SEED,
        ".tsv"
      )
    ),
    sep = "\t"
  )


  cat(
    "XtX prepared: ",
    nrow(XTX_OUT),
    " SNPs\n",
    sep = ""
  )


  ##########################################################
  ## M_P lookup
  ##########################################################

  MP <- XTX[
    ,
    .(
      BayPass_row,
      M_P
    )
  ]


  ##########################################################
  ## 4.3 C2 from representative seed
  ##########################################################

  C2_RAW <- read_c2(
    dataset,
    REP_SEED
  )


  if (length(unique(C2_RAW$CONTRAST)) != 1) {

    stop(
      "Expected exactly one C2 contrast for ",
      dataset
    )
  }


  if (nrow(C2_RAW) != nrow(META)) {

    stop(
      "C2 and metadata SNP counts differ for ",
      dataset
    )
  }


  C2 <- merge(
    META,
    C2_RAW,
    by.x = "BayPass_row",
    by.y = "MRK",
    all.x = TRUE,
    sort = FALSE
  )


  C2 <- merge(
    C2,
    MP,
    by = "BayPass_row",
    all.x = TRUE,
    sort = FALSE
  )


  if (
    anyNA(C2$C2) ||
    anyNA(C2$M_P)
  ) {

    stop(
      "Missing C2 or M_P values after merge: ",
      dataset
    )
  }


  setorder(
    C2,
    CHROM,
    POS
  )


  C2_OUT <- C2[
    ,
    .(
      BayPass_row,
      CHROM,
      SNP,
      POS,
      REF,
      ALT,
      M_P,
      CONTRAST,
      C2,
      log10_p = `log10(1/pval)`
    )
  ]


  fwrite(
    C2_OUT,
    file.path(
      out_ds,
      paste0(
        "C2_Host_seed",
        REP_SEED,
        ".tsv"
      )
    ),
    sep = "\t"
  )


  cat(
    "C2 prepared: ",
    nrow(C2_OUT),
    " SNPs\n",
    sep = ""
  )


  ##########################################################
  ## 4.4 Read GEA BF for three seeds
  ##########################################################

  BF_SEEDS <- list()


  for (seed in SEEDS) {

    cat(
      "Reading BF seed ",
      seed,
      "...\n",
      sep = ""
    )


    tmp <- read_bf(
      dataset,
      seed
    )


    new_name <- paste0(
      "BF_seed",
      seed
    )


    setnames(
      tmp,
      "BFdb",
      new_name
    )


    BF_SEEDS[[seed]] <- tmp
  }


  ##########################################################
  ## 4.5 Merge seed-specific BFs
  ##########################################################

  BF <- Reduce(
    function(x, y) {

      merge(
        x,
        y,
        by = c(
          "MRK",
          "COVARIABLE"
        ),
        all = TRUE,
        sort = FALSE
      )

    },
    BF_SEEDS
  )


  bf_seed_cols <- paste0(
    "BF_seed",
    SEEDS
  )


  if (
    anyNA(
      BF[, ..bf_seed_cols]
    )
  ) {

    stop(
      "Missing BF values among seeds for ",
      dataset
    )
  }


  ##########################################################
  ## 4.6 Mean / SD / min / max across seeds
  ##########################################################

  BF[
    ,
    BF_mean := rowMeans(
      .SD
    ),
    .SDcols = bf_seed_cols
  ]


  BF[
    ,
    BF_sd := apply(
      .SD,
      1,
      sd
    ),
    .SDcols = bf_seed_cols
  ]


  BF[
    ,
    BF_min := apply(
      .SD,
      1,
      min
    ),
    .SDcols = bf_seed_cols
  ]


  BF[
    ,
    BF_max := apply(
      .SD,
      1,
      max
    ),
    .SDcols = bf_seed_cols
  ]


  ##########################################################
  ## 4.7 Add covariate names
  ##########################################################

  BF <- merge(
    BF,
    COV_MAP,
    by = "COVARIABLE",
    all.x = TRUE,
    sort = FALSE
  )


  if (anyNA(BF$Covariate)) {

    stop(
      "Unknown covariate number in ",
      dataset
    )
  }


  ##########################################################
  ## 4.8 Add SNP metadata
  ##########################################################

  META_BF <- copy(
    META
  )


  setnames(
    META_BF,
    "BayPass_row",
    "MRK"
  )


  BF <- merge(
    BF,
    META_BF,
    by = "MRK",
    all.x = TRUE,
    sort = FALSE
  )


  setnames(
    BF,
    "MRK",
    "BayPass_row"
  )


  ##########################################################
  ## 4.9 Add M_P
  ##########################################################

  BF <- merge(
    BF,
    MP,
    by = "BayPass_row",
    all.x = TRUE,
    sort = FALSE
  )


  if (
    anyNA(BF$CHROM) ||
    anyNA(BF$POS) ||
    anyNA(BF$M_P)
  ) {

    stop(
      "Metadata/M_P merge failed for ",
      dataset
    )
  }


  ##########################################################
  ## Expected BF row count
  ##########################################################

  expected_bf_rows <-
    nrow(META) *
    nrow(COV_MAP)


  if (nrow(BF) != expected_bf_rows) {

    stop(
      "Unexpected BF row count for ",
      dataset,
      ": found ",
      nrow(BF),
      ", expected ",
      expected_bf_rows
    )
  }


  ##########################################################
  ## 4.10 BF correlations among seeds
  ##########################################################

  seed_pairs <- combn(
    SEEDS,
    2,
    simplify = FALSE
  )


  for (cov in COV_MAP$COVARIABLE) {

    cov_name <- COV_MAP[
      COVARIABLE == cov,
      Covariate
    ]


    bf_cov <- BF[
      COVARIABLE == cov
    ]


    for (pair in seed_pairs) {

      s1 <- pair[1]

      s2 <- pair[2]


      col1 <- paste0(
        "BF_seed",
        s1
      )


      col2 <- paste0(
        "BF_seed",
        s2
      )


      BF_COR_RESULTS[[counter]] <- data.table(
        Dataset = dataset,
        COVARIABLE = cov,
        Covariate = cov_name,
        Seed1 = s1,
        Seed2 = s2,
        N_SNPs = nrow(bf_cov),
        Spearman_rho = cor(
          bf_cov[[col1]],
          bf_cov[[col2]],
          method = "spearman",
          use = "complete.obs"
        ),
        Pearson_r = cor(
          bf_cov[[col1]],
          bf_cov[[col2]],
          method = "pearson",
          use = "complete.obs"
        )
      )


      counter <- counter + 1
    }
  }


  ##########################################################
  ## 4.11 Sort final BF table
  ##########################################################

  setorder(
    BF,
    COVARIABLE,
    CHROM,
    POS
  )


  ##########################################################
  ## 4.12 Write final mean-BF table
  ##########################################################

  BF_OUT <- BF[
    ,
    .(
      BayPass_row,
      CHROM,
      SNP,
      POS,
      REF,
      ALT,
      M_P,
      COVARIABLE,
      Covariate,
      BF_seed12345,
      BF_seed67890,
      BF_seed13579,
      BF_mean,
      BF_sd,
      BF_min,
      BF_max
    )
  ]


  fwrite(
    BF_OUT,
    file.path(
      out_ds,
      "GEA_BF_mean3seeds.tsv"
    ),
    sep = "\t"
  )


  cat(
    "BF prepared: ",
    nrow(BF_OUT),
    " SNP x covariate rows\n",
    sep = ""
  )


  cat(
    "Finished dataset: ",
    dataset,
    "\n",
    sep = ""
  )
}


############################################################
## 5. COMBINE BF CORRELATIONS
############################################################

BF_COR <- rbindlist(
  BF_COR_RESULTS,
  use.names = TRUE
)


setorder(
  BF_COR,
  Dataset,
  COVARIABLE,
  Seed1,
  Seed2
)


fwrite(
  BF_COR,
  file.path(
    OUTBASE,
    "BF_seed_correlations.tsv"
  ),
  sep = "\t"
)


############################################################
## 6. BF CORRELATION SUMMARY
############################################################

BF_COR_SUMMARY <- BF_COR[
  ,
  .(
    Mean_Spearman = mean(
      Spearman_rho
    ),

    Median_Spearman = median(
      Spearman_rho
    ),

    Min_Spearman = min(
      Spearman_rho
    ),

    Max_Spearman = max(
      Spearman_rho
    ),

    Mean_Pearson = mean(
      Pearson_r
    )
  ),
  by = .(
    Dataset,
    COVARIABLE,
    Covariate
  )
]


setorder(
  BF_COR_SUMMARY,
  Dataset,
  COVARIABLE
)


fwrite(
  BF_COR_SUMMARY,
  file.path(
    OUTBASE,
    "BF_seed_correlations_summary.tsv"
  ),
  sep = "\t"
)


############################################################
## 7. PRINT SUMMARY
############################################################

cat(
  "\n========================================\n"
)

cat(
  "BF SEED CORRELATION SUMMARY\n"
)

cat(
  "========================================\n"
)


print(
  BF_COR_SUMMARY
)


cat(
  "\n========================================\n"
)

cat(
  "OUTPUT DIRECTORY\n"
)

cat(
  "========================================\n"
)

cat(
  OUTBASE,
  "\n"
)


cat(
  "\nDONE\n"
)
