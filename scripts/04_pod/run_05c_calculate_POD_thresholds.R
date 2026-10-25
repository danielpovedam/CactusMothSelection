#!/usr/bin/env Rscript

############################################################
## 49C_calculate_POD_thresholds_and_outliers.R
##
## Cactoblastis cactorum
##
## FINAL POD calibration for BayPass XtX and C2
##
## PURPOSE
## -------
## 1. Read the neutral 100,000-SNP POD simulated in 49A.
##
## 2. Calculate observed pooled MAF for every simulated SNP
##    directly from the BayPass allele-count matrix.
##
## 3. Define two POD calibration sets:
##
##      A. ALL POD loci
##      B. ASCERTAINED POD loci with observed MAF >= 0.05
##
##    The second set is the PRIMARY calibration because the
##    empirical all-SNP BayPass dataset had already been
##    filtered to MAF >= 0.05.
##
## 4. Calculate 99%, 99.5%, and 99.9% POD thresholds for:
##
##      XtXst
##      C2
##
## 5. Freeze the PRIMARY SNP-outlier criterion as:
##
##      XtXst > 99.9th percentile of MAF-filtered POD
##      C2    > 99.9th percentile of MAF-filtered POD
##
## 6. Apply these thresholds to the empirical 132,283 SNPs.
##
## IMPORTANT
## ---------
## These POD-calibrated SNPs are a complementary SNP-level
## layer. They DO NOT replace the local-score candidate
## regions used as the primary region-level inference.
##
## No cross-analysis overlap is calculated here.
############################################################


suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


options(stringsAsFactors = FALSE)


############################################################
## 0. PATHS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"


############################################################
## POD simulation input from 49A
############################################################

POD_GENO_FILE <- file.path(
  BASE,
  "POD_analysis",
  "01_simulation",
  "G.POD_neutral_100k"
)


############################################################
## POD BayPass results from 49B
############################################################

POD_PREFIX <- file.path(
  BASE,
  "POD_analysis",
  "02_baypass_run",
  "BP_POD_neutral100k_seed24680"
)


POD_XTX_FILE <- paste0(
  POD_PREFIX,
  "_summary_pi_xtx.out"
)


POD_C2_FILE <- paste0(
  POD_PREFIX,
  "_summary_contrast.out"
)


############################################################
## Empirical BayPass results
##
## Representative converged all-SNP run:
## seed 12345
############################################################

EMP_PREFIX <- file.path(
  BASE,
  "runs_10x",
  "all_snps",
  "seed_12345",
  "BP_all_snps_seed12345"
)


EMP_XTX_FILE <- paste0(
  EMP_PREFIX,
  "_summary_pi_xtx.out"
)


EMP_C2_FILE <- paste0(
  EMP_PREFIX,
  "_summary_contrast.out"
)


############################################################
## Empirical SNP metadata
############################################################

META_FILE <- file.path(
  "/home/dpoveda/micetral/Population_snps",
  "baypass_preparation",
  "metadata",
  "all_snps.baypass.snp_metadata.tsv"
)


############################################################
## Output directories
############################################################

OUT_THRESHOLD <- file.path(
  BASE,
  "POD_analysis",
  "03_thresholds"
)


OUT_OUTLIER <- file.path(
  BASE,
  "POD_analysis",
  "04_empirical_outliers"
)


OUT_FIG <- file.path(
  BASE,
  "POD_analysis",
  "05_figures"
)


OUT_TABLE <- file.path(
  BASE,
  "POD_analysis",
  "06_tables"
)


OUT_REPRO <- file.path(
  BASE,
  "POD_analysis",
  "reproducibility"
)


for (
  d in c(
    OUT_THRESHOLD,
    OUT_OUTLIER,
    OUT_FIG,
    OUT_TABLE,
    OUT_REPRO
  )
) {

  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


############################################################
## 1. SETTINGS — FROZEN BEFORE OVERLAP ANALYSIS
############################################################

N_EMPIRICAL <- 132283L

N_POD_EXPECTED <- 100000L

NPOP <- 28L

MAF_THRESHOLD <- 0.05


############################################################
## POD thresholds to report
############################################################

QUANTILES <- c(
  0.990,
  0.995,
  0.999
)


############################################################
## PRIMARY FINAL threshold
############################################################

PRIMARY_QUANTILE <- 0.999


############################################################
## PRIMARY calibration set
##
## TRUE =
## thresholds are obtained from POD loci whose OBSERVED
## pooled MAF is >= 0.05.
############################################################

PRIMARY_USE_MAF_FILTER <- TRUE


############################################################
## 2. START REPORT
############################################################

cat(
  "\n============================================================\n"
)

cat(
  "49C — FINAL POD CALIBRATION\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Primary calibration:\n"
)

cat(
  "  POD loci with observed MAF >= ",
  MAF_THRESHOLD,
  "\n",
  sep = ""
)


cat(
  "Primary quantile:\n"
)

cat(
  "  ",
  PRIMARY_QUANTILE,
  "\n\n",
  sep = ""
)


############################################################
## 3. INPUT QC
############################################################

FILES <- c(
  POD_GENO_FILE,
  POD_XTX_FILE,
  POD_C2_FILE,
  EMP_XTX_FILE,
  EMP_C2_FILE,
  META_FILE
)


missing_files <- FILES[
  !file.exists(FILES)
]


if (
  length(missing_files) > 0L
) {

  stop(
    "\nMissing required input file(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    "\n"
  )
}


empty_files <- FILES[
  file.info(FILES)$size == 0
]


if (
  length(empty_files) > 0L
) {

  stop(
    "\nEmpty required input file(s):\n",
    paste(
      empty_files,
      collapse = "\n"
    ),
    "\n"
  )
}


############################################################
## 4. READ POD COUNT MATRIX
############################################################

cat(
  "------------------------------------------------------------\n"
)

cat(
  "Reading POD count matrix...\n"
)


POD_G <- fread(
  POD_GENO_FILE,
  header = FALSE,
  showProgress = TRUE
)


if (
  nrow(POD_G) != N_POD_EXPECTED
) {

  stop(
    "Expected ",
    N_POD_EXPECTED,
    " POD loci; found ",
    nrow(POD_G)
  )
}


if (
  ncol(POD_G) != 2L * NPOP
) {

  stop(
    "Expected ",
    2L * NPOP,
    " BayPass allele-count columns; found ",
    ncol(POD_G)
  )
}


############################################################
## Ensure numeric
############################################################

for (
  j in seq_len(
    ncol(POD_G)
  )
) {

  set(
    POD_G,
    j = j,
    value = as.numeric(
      POD_G[[j]]
    )
  )
}


if (
  anyNA(
    POD_G
  )
) {

  stop(
    "Missing values detected in POD genotype matrix."
  )
}


if (
  any(
    as.matrix(POD_G) < 0
  )
) {

  stop(
    "Negative allele counts detected in POD genotype matrix."
  )
}


############################################################
## 5. CALCULATE OBSERVED MAF FOR EACH POD SNP
##
## BayPass count structure:
##
## population 1:
##   col 1 = allele A
##   col 2 = allele B
##
## population 2:
##   col 3 = allele A
##   col 4 = allele B
##
## ...
############################################################

ODD_COLS <- seq(
  1L,
  2L * NPOP,
  by = 2L
)


EVEN_COLS <- seq(
  2L,
  2L * NPOP,
  by = 2L
)


ALLELE1_TOTAL <- rowSums(
  as.matrix(
    POD_G[
      ,
      ..ODD_COLS
    ]
  )
)


ALLELE2_TOTAL <- rowSums(
  as.matrix(
    POD_G[
      ,
      ..EVEN_COLS
    ]
  )
)


TOTAL_ALLELES <- ALLELE1_TOTAL +
  ALLELE2_TOTAL


if (
  any(
    TOTAL_ALLELES <= 0
  )
) {

  stop(
    "At least one POD locus has no observed allele counts."
  )
}


POD_FREQ <- ALLELE1_TOTAL /
  TOTAL_ALLELES


POD_MAF <- pmin(
  POD_FREQ,
  1 - POD_FREQ
)


POD_FIXED <- POD_MAF == 0


POD_PASS_MAF <- POD_MAF >=
  MAF_THRESHOLD


############################################################
## POD ascertainment table
############################################################

POD_MAF_DT <- data.table(

  MRK =
    seq_len(
      N_POD_EXPECTED
    ),

  Allele1_total =
    ALLELE1_TOTAL,

  Allele2_total =
    ALLELE2_TOTAL,

  Total_alleles =
    TOTAL_ALLELES,

  Observed_frequency =
    POD_FREQ,

  Observed_MAF =
    POD_MAF,

  Fixed =
    POD_FIXED,

  Pass_MAF_0.05 =
    POD_PASS_MAF
)


fwrite(
  POD_MAF_DT,
  file.path(
    OUT_TABLE,
    "Table_POD_observed_MAF_all100k.tsv"
  ),
  sep = "\t"
)


############################################################
## POD ascertainment QC
############################################################

N_FIXED <- sum(
  POD_FIXED
)


N_MAF_LT05 <- sum(
  POD_MAF < MAF_THRESHOLD
)


N_MAF_GE05 <- sum(
  POD_PASS_MAF
)


cat(
  "\nPOD ascertainment:\n"
)


cat(
  "  Total POD SNPs       = ",
  length(POD_MAF),
  "\n",
  sep = ""
)


cat(
  "  Fixed POD SNPs       = ",
  N_FIXED,
  "\n",
  sep = ""
)


cat(
  "  MAF < 0.05           = ",
  N_MAF_LT05,
  "\n",
  sep = ""
)


cat(
  "  MAF >= 0.05          = ",
  N_MAF_GE05,
  "\n\n",
  sep = ""
)


if (
  N_MAF_GE05 < 50000L
) {

  warning(
    "Fewer than 50,000 POD SNPs remain after MAF ascertainment."
  )
}


############################################################
## 6. READ POD XTX RESULTS
############################################################

cat(
  "Reading POD XtX output...\n"
)


POD_XTX <- fread(
  POD_XTX_FILE
)


required_xtx_columns <- c(
  "MRK",
  "M_XtX",
  "XtXst"
)


if (
  !all(
    required_xtx_columns %in%
      names(POD_XTX)
  )
) {

  stop(
    "\nPOD XtX output is missing required columns.\n",
    "Found:\n",
    paste(
      names(POD_XTX),
      collapse = ", "
    )
  )
}


if (
  nrow(POD_XTX) != N_POD_EXPECTED
) {

  stop(
    "Expected ",
    N_POD_EXPECTED,
    " POD XtX rows; found ",
    nrow(POD_XTX)
  )
}


POD_XTX[
  ,
  MRK :=
    as.integer(
      MRK
    )
]


if (
  anyDuplicated(
    POD_XTX$MRK
  )
) {

  stop(
    "Duplicated MRK IDs in POD XtX output."
  )
}


if (
  anyNA(
    POD_XTX$XtXst
  )
) {

  stop(
    "Missing XtXst values in POD XtX output."
  )
}


############################################################
## 7. READ POD C2 RESULTS
############################################################

cat(
  "Reading POD C2 output...\n"
)


POD_C2 <- fread(
  POD_C2_FILE
)


required_c2_columns <- c(
  "CONTRAST",
  "MRK",
  "C2"
)


if (
  !all(
    required_c2_columns %in%
      names(POD_C2)
  )
) {

  stop(
    "\nPOD C2 output is missing required columns.\n",
    "Found:\n",
    paste(
      names(POD_C2),
      collapse = ", "
    )
  )
}


if (
  uniqueN(
    POD_C2$CONTRAST
  ) != 1L
) {

  stop(
    "Expected exactly one host C2 contrast in POD output."
  )
}


if (
  nrow(POD_C2) != N_POD_EXPECTED
) {

  stop(
    "Expected ",
    N_POD_EXPECTED,
    " POD C2 rows; found ",
    nrow(POD_C2)
  )
}


POD_C2[
  ,
  MRK :=
    as.integer(
      MRK
    )
]


if (
  anyDuplicated(
    POD_C2$MRK
  )
) {

  stop(
    "Duplicated MRK IDs in POD C2 output."
  )
}


if (
  anyNA(
    POD_C2$C2
  )
) {

  stop(
    "Missing C2 values in POD C2 output."
  )
}


############################################################
## 8. JOIN POD MAF TO BAYPASS RESULTS
############################################################

POD_XTX <- merge(
  POD_XTX,
  POD_MAF_DT[
    ,
    .(
      MRK,
      Observed_frequency,
      Observed_MAF,
      Fixed,
      Pass_MAF_0.05
    )
  ],
  by = "MRK",
  all.x = TRUE,
  sort = FALSE
)


POD_C2 <- merge(
  POD_C2,
  POD_MAF_DT[
    ,
    .(
      MRK,
      Observed_frequency,
      Observed_MAF,
      Fixed,
      Pass_MAF_0.05
    )
  ],
  by = "MRK",
  all.x = TRUE,
  sort = FALSE
)


if (
  anyNA(
    POD_XTX$Observed_MAF
  ) ||
  anyNA(
    POD_C2$Observed_MAF
  )
) {

  stop(
    "Failed to map POD MAF information to BayPass results."
  )
}


############################################################
## 9. DEFINE POD CALIBRATION SETS
############################################################

POD_XTX_ALL <- POD_XTX


POD_C2_ALL <- POD_C2


POD_XTX_MAF <- POD_XTX[
  Pass_MAF_0.05 == TRUE
]


POD_C2_MAF <- POD_C2[
  Pass_MAF_0.05 == TRUE
]


if (
  nrow(POD_XTX_MAF) !=
    N_MAF_GE05
) {

  stop(
    "XtX MAF-filtered POD count does not match MAF QC count."
  )
}


if (
  nrow(POD_C2_MAF) !=
    N_MAF_GE05
) {

  stop(
    "C2 MAF-filtered POD count does not match MAF QC count."
  )
}


############################################################
## 10. CALCULATE ALL THRESHOLDS
##
## We calculate thresholds for:
##
## A. all POD SNPs
## B. MAF >= 0.05 POD SNPs
##
## Quantile type 8:
## approximately median-unbiased quantile estimator.
############################################################

calculate_thresholds <- function(
  x,
  statistic,
  calibration_set,
  n_calibration
) {

  qs <- as.numeric(
    quantile(
      x,
      probs = QUANTILES,
      na.rm = TRUE,
      names = FALSE,
      type = 8
    )
  )


  data.table(

    Statistic =
      statistic,

    Calibration_set =
      calibration_set,

    N_POD_SNPs =
      n_calibration,

    Quantile =
      QUANTILES,

    Threshold =
      qs
  )
}


THRESHOLDS <- rbindlist(
  list(

    calculate_thresholds(
      POD_XTX_ALL$XtXst,
      "XtXst",
      "All_POD",
      nrow(POD_XTX_ALL)
    ),

    calculate_thresholds(
      POD_XTX_MAF$XtXst,
      "XtXst",
      "MAF_ge_0.05",
      nrow(POD_XTX_MAF)
    ),

    calculate_thresholds(
      POD_C2_ALL$C2,
      "C2",
      "All_POD",
      nrow(POD_C2_ALL)
    ),

    calculate_thresholds(
      POD_C2_MAF$C2,
      "C2",
      "MAF_ge_0.05",
      nrow(POD_C2_MAF)
    )
  )
)


THRESHOLDS[
  ,
  Primary :=
    Calibration_set ==
      "MAF_ge_0.05" &
    Quantile ==
      PRIMARY_QUANTILE
]


setorder(
  THRESHOLDS,
  Statistic,
  Calibration_set,
  Quantile
)


fwrite(
  THRESHOLDS,
  file.path(
    OUT_THRESHOLD,
    "Table_POD_thresholds_XtXst_C2_ALL_vs_MAFfiltered.tsv"
  ),
  sep = "\t"
)


cat(
  "\n------------------------------------------------------------\n"
)

cat(
  "POD THRESHOLDS\n"
)

cat(
  "------------------------------------------------------------\n"
)


print(
  THRESHOLDS
)


############################################################
## 11. PRIMARY THRESHOLDS
############################################################

XTX_PRIMARY <- THRESHOLDS[

  Statistic == "XtXst" &

    Calibration_set ==
      "MAF_ge_0.05" &

    Quantile ==
      PRIMARY_QUANTILE,

  Threshold
]


C2_PRIMARY <- THRESHOLDS[

  Statistic == "C2" &

    Calibration_set ==
      "MAF_ge_0.05" &

    Quantile ==
      PRIMARY_QUANTILE,

  Threshold
]


if (
  length(XTX_PRIMARY) != 1L ||
  length(C2_PRIMARY) != 1L
) {

  stop(
    "Failed to identify unique primary POD thresholds."
  )
}


cat(
  "\nPRIMARY 99.9% MAF-FILTERED POD THRESHOLDS\n"
)


cat(
  "  XtXst = ",
  XTX_PRIMARY,
  "\n",
  sep = ""
)


cat(
  "  C2    = ",
  C2_PRIMARY,
  "\n\n",
  sep = ""
)


############################################################
## 12. COMPARE ALL-POD VS MAF-FILTERED THRESHOLDS
############################################################

THRESHOLD_COMPARE <- dcast(
  THRESHOLDS,
  Statistic + Quantile ~ Calibration_set,
  value.var = "Threshold"
)


THRESHOLD_COMPARE[
  ,
  Absolute_difference :=
    MAF_ge_0.05 -
    All_POD
]


THRESHOLD_COMPARE[
  ,
  Percent_difference :=
    100 *
    (
      MAF_ge_0.05 -
      All_POD
    ) /
    All_POD
]


fwrite(
  THRESHOLD_COMPARE,
  file.path(
    OUT_TABLE,
    "Table_POD_threshold_comparison_ALL_vs_MAFfiltered.tsv"
  ),
  sep = "\t"
)


############################################################
## 13. READ EMPIRICAL XTX
############################################################

cat(
  "Reading empirical XtX output...\n"
)


EMP_XTX <- fread(
  EMP_XTX_FILE
)


if (
  !all(
    c(
      "MRK",
      "M_XtX",
      "XtXst"
    ) %in%
      names(EMP_XTX)
  )
) {

  stop(
    "Unexpected empirical XtX file structure."
  )
}


if (
  nrow(EMP_XTX) != N_EMPIRICAL
) {

  stop(
    "Expected ",
    N_EMPIRICAL,
    " empirical XtX SNPs; found ",
    nrow(EMP_XTX)
  )
}


EMP_XTX[
  ,
  MRK :=
    as.integer(
      MRK
    )
]


if (
  anyDuplicated(
    EMP_XTX$MRK
  )
) {

  stop(
    "Duplicated empirical XtX MRKs."
  )
}


############################################################
## 14. READ EMPIRICAL C2
############################################################

cat(
  "Reading empirical C2 output...\n"
)


EMP_C2 <- fread(
  EMP_C2_FILE
)


if (
  !all(
    c(
      "CONTRAST",
      "MRK",
      "C2"
    ) %in%
      names(EMP_C2)
  )
) {

  stop(
    "Unexpected empirical C2 file structure."
  )
}


if (
  uniqueN(
    EMP_C2$CONTRAST
  ) != 1L
) {

  stop(
    "Expected exactly one empirical C2 contrast."
  )
}


if (
  nrow(EMP_C2) != N_EMPIRICAL
) {

  stop(
    "Expected ",
    N_EMPIRICAL,
    " empirical C2 SNPs; found ",
    nrow(EMP_C2)
  )
}


EMP_C2[
  ,
  MRK :=
    as.integer(
      MRK
    )
]


if (
  anyDuplicated(
    EMP_C2$MRK
  )
) {

  stop(
    "Duplicated empirical C2 MRKs."
  )
}


############################################################
## 15. READ SNP METADATA
############################################################

META <- fread(
  META_FILE
)


required_meta <- c(
  "BayPass_row",
  "CHROM",
  "SNP",
  "POS",
  "ALT",
  "REF"
)


if (
  !all(
    required_meta %in%
      names(META)
  )
) {

  stop(
    "\nMetadata missing required columns.\n",
    "Expected:\n",
    paste(
      required_meta,
      collapse = ", "
    ),
    "\nFound:\n",
    paste(
      names(META),
      collapse = ", "
    )
  )
}


if (
  nrow(META) != N_EMPIRICAL
) {

  stop(
    "Expected ",
    N_EMPIRICAL,
    " metadata rows; found ",
    nrow(META)
  )
}


META[
  ,
  BayPass_row :=
    as.integer(
      BayPass_row
    )
]


if (
  anyDuplicated(
    META$BayPass_row
  )
) {

  stop(
    "BayPass_row is not unique in metadata."
  )
}


############################################################
## 16. FLAG EMPIRICAL OUTLIERS AT ALL THRESHOLDS
##
## We save flags for:
##
## 99.0%
## 99.5%
## 99.9%
##
## using the PRIMARY MAF-filtered POD calibration.
############################################################

for (
  q in QUANTILES
) {

  label <- paste0(
    "POD_MAF_q",
    sprintf(
      "%03d",
      round(
        q * 1000
      )
    )
  )


  xtx_thr <- THRESHOLDS[

    Statistic == "XtXst" &

      Calibration_set ==
        "MAF_ge_0.05" &

      Quantile == q,

    Threshold
  ]


  c2_thr <- THRESHOLDS[

    Statistic == "C2" &

      Calibration_set ==
        "MAF_ge_0.05" &

      Quantile == q,

    Threshold
  ]


  EMP_XTX[
    ,
    (label) :=
      XtXst >
      xtx_thr
  ]


  EMP_C2[
    ,
    (label) :=
      C2 >
      c2_thr
  ]
}


############################################################
## Explicit primary flag
############################################################

EMP_XTX[
  ,
  POD999_primary_outlier :=
    XtXst >
    XTX_PRIMARY
]


EMP_C2[
  ,
  POD999_primary_outlier :=
    C2 >
    C2_PRIMARY
]


############################################################
## 17. OUTLIER COUNTS AT EACH THRESHOLD
############################################################

OUTLIER_SUMMARY <- rbindlist(
  lapply(
    QUANTILES,
    function(q) {

      label <- paste0(
        "POD_MAF_q",
        sprintf(
          "%03d",
          round(
            q * 1000
          )
        )
      )


      xtx_thr <- THRESHOLDS[

        Statistic == "XtXst" &

          Calibration_set ==
            "MAF_ge_0.05" &

          Quantile == q,

        Threshold
      ]


      c2_thr <- THRESHOLDS[

        Statistic == "C2" &

          Calibration_set ==
            "MAF_ge_0.05" &

          Quantile == q,

        Threshold
      ]


      data.table(

        Statistic =
          c(
            "XtXst",
            "C2"
          ),

        POD_calibration =
          "Observed MAF >= 0.05",

        POD_quantile =
          q,

        POD_threshold =
          c(
            xtx_thr,
            c2_thr
          ),

        N_empirical_SNPs =
          N_EMPIRICAL,

        N_empirical_outliers =
          c(
            sum(
              EMP_XTX[[label]],
              na.rm = TRUE
            ),

            sum(
              EMP_C2[[label]],
              na.rm = TRUE
            )
          ),

        Percent_empirical_outliers =
          100 *
          c(
            mean(
              EMP_XTX[[label]],
              na.rm = TRUE
            ),

            mean(
              EMP_C2[[label]],
              na.rm = TRUE
            )
          ),

        Primary =
          q ==
            PRIMARY_QUANTILE
      )
    }
  )
)


fwrite(
  OUTLIER_SUMMARY,
  file.path(
    OUT_TABLE,
    "Table_POD_empirical_outlier_counts_by_threshold.tsv"
  ),
  sep = "\t"
)


############################################################
## 18. EXTRACT PRIMARY 99.9% OUTLIERS
############################################################

XTX_OUT <- EMP_XTX[
  POD999_primary_outlier == TRUE
]


C2_OUT <- EMP_C2[
  POD999_primary_outlier == TRUE
]


############################################################
## 19. ADD GENOMIC METADATA
############################################################

XTX_OUT_META <- merge(
  XTX_OUT,
  META,
  by.x = "MRK",
  by.y = "BayPass_row",
  all.x = TRUE,
  sort = FALSE
)


C2_OUT_META <- merge(
  C2_OUT,
  META,
  by.x = "MRK",
  by.y = "BayPass_row",
  all.x = TRUE,
  sort = FALSE
)


if (
  nrow(XTX_OUT_META) !=
    nrow(XTX_OUT)
) {

  stop(
    "XtX metadata merge changed row number."
  )
}


if (
  nrow(C2_OUT_META) !=
    nrow(C2_OUT)
) {

  stop(
    "C2 metadata merge changed row number."
  )
}


if (
  anyNA(
    XTX_OUT_META$CHROM
  )
) {

  stop(
    "Missing genomic metadata among XtX POD outliers."
  )
}


if (
  anyNA(
    C2_OUT_META$CHROM
  )
) {

  stop(
    "Missing genomic metadata among C2 POD outliers."
  )
}


setorder(
  XTX_OUT_META,
  CHROM,
  POS
)


setorder(
  C2_OUT_META,
  CHROM,
  POS
)


############################################################
## 20. SAVE PRIMARY OUTLIER TABLES
############################################################

fwrite(
  XTX_OUT_META,
  file.path(
    OUT_OUTLIER,
    "XtX_POD999_MAFfiltered_outlier_SNPs.tsv"
  ),
  sep = "\t"
)


fwrite(
  C2_OUT_META,
  file.path(
    OUT_OUTLIER,
    "C2_POD999_MAFfiltered_outlier_SNPs.tsv"
  ),
  sep = "\t"
)


############################################################
## 21. SAVE COMPLETE EMPIRICAL TABLES WITH FLAGS
############################################################

fwrite(
  EMP_XTX,
  file.path(
    OUT_OUTLIER,
    "XtX_all132283_with_POD_MAF_threshold_flags.tsv"
  ),
  sep = "\t"
)


fwrite(
  EMP_C2,
  file.path(
    OUT_OUTLIER,
    "C2_all132283_with_POD_MAF_threshold_flags.tsv"
  ),
  sep = "\t"
)


############################################################
## 22. POD ASCERTAINMENT SUMMARY
############################################################

ASCERTAINMENT_SUMMARY <- data.table(

  Metric = c(
    "Total_POD_SNPs",
    "Fixed_POD_SNPs",
    "POD_SNPs_MAF_lt_0.05",
    "POD_SNPs_MAF_ge_0.05",
    "Percent_POD_retained_MAF_ge_0.05",
    "Observed_MAF_min",
    "Observed_MAF_median",
    "Observed_MAF_mean",
    "Observed_MAF_max"
  ),

  Value = c(
    N_POD_EXPECTED,
    N_FIXED,
    N_MAF_LT05,
    N_MAF_GE05,
    100 *
      N_MAF_GE05 /
      N_POD_EXPECTED,
    min(POD_MAF),
    median(POD_MAF),
    mean(POD_MAF),
    max(POD_MAF)
  )
)


fwrite(
  ASCERTAINMENT_SUMMARY,
  file.path(
    OUT_TABLE,
    "Table_POD_MAF_ascertainment_summary.tsv"
  ),
  sep = "\t"
)


############################################################
## 23. DISTRIBUTION SUMMARY
############################################################

DISTRIBUTION_SUMMARY <- rbindlist(
  list(

    data.table(

      Statistic = "XtXst",

      Calibration_set = "All_POD",

      N = nrow(POD_XTX_ALL),

      Mean = mean(
        POD_XTX_ALL$XtXst
      ),

      SD = sd(
        POD_XTX_ALL$XtXst
      ),

      Median = median(
        POD_XTX_ALL$XtXst
      ),

      Min = min(
        POD_XTX_ALL$XtXst
      ),

      Max = max(
        POD_XTX_ALL$XtXst
      )
    ),


    data.table(

      Statistic = "XtXst",

      Calibration_set = "MAF_ge_0.05",

      N = nrow(POD_XTX_MAF),

      Mean = mean(
        POD_XTX_MAF$XtXst
      ),

      SD = sd(
        POD_XTX_MAF$XtXst
      ),

      Median = median(
        POD_XTX_MAF$XtXst
      ),

      Min = min(
        POD_XTX_MAF$XtXst
      ),

      Max = max(
        POD_XTX_MAF$XtXst
      )
    ),


    data.table(

      Statistic = "C2",

      Calibration_set = "All_POD",

      N = nrow(POD_C2_ALL),

      Mean = mean(
        POD_C2_ALL$C2
      ),

      SD = sd(
        POD_C2_ALL$C2
      ),

      Median = median(
        POD_C2_ALL$C2
      ),

      Min = min(
        POD_C2_ALL$C2
      ),

      Max = max(
        POD_C2_ALL$C2
      )
    ),


    data.table(

      Statistic = "C2",

      Calibration_set = "MAF_ge_0.05",

      N = nrow(POD_C2_MAF),

      Mean = mean(
        POD_C2_MAF$C2
      ),

      SD = sd(
        POD_C2_MAF$C2
      ),

      Median = median(
        POD_C2_MAF$C2
      ),

      Min = min(
        POD_C2_MAF$C2
      ),

      Max = max(
        POD_C2_MAF$C2
      )
    )
  )
)


fwrite(
  DISTRIBUTION_SUMMARY,
  file.path(
    OUT_TABLE,
    "Table_POD_statistic_distribution_summary.tsv"
  ),
  sep = "\t"
)


############################################################
## 24. FIGURE DATA
############################################################

FIG_XTX <- rbindlist(
  list(

    data.table(
      Value =
        POD_XTX_ALL$XtXst,
      Calibration =
        "All POD"
    ),

    data.table(
      Value =
        POD_XTX_MAF$XtXst,
      Calibration =
        "MAF >= 0.05"
    )
  )
)


FIG_C2 <- rbindlist(
  list(

    data.table(
      Value =
        POD_C2_ALL$C2,
      Calibration =
        "All POD"
    ),

    data.table(
      Value =
        POD_C2_MAF$C2,
      Calibration =
        "MAF >= 0.05"
    )
  )
)


############################################################
## 25. SUPPLEMENTARY FIGURE A — XTX
############################################################

PX <- ggplot(
  FIG_XTX,
  aes(
    x = Value,
    linetype = Calibration
  )
) +

  geom_density(
    linewidth = 0.8,
    adjust = 1
  ) +

  geom_vline(
    xintercept =
      XTX_PRIMARY,
    linewidth = 0.7,
    linetype = 2
  ) +

  labs(
    title = "A. XtX POD calibration",
    x = expression(XtX[st]),
    y = "Density",
    linetype = NULL
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    legend.position =
      "top"
  )


############################################################
## 26. SUPPLEMENTARY FIGURE B — C2
############################################################

PC <- ggplot(
  FIG_C2,
  aes(
    x = Value,
    linetype = Calibration
  )
) +

  geom_density(
    linewidth = 0.8,
    adjust = 1
  ) +

  geom_vline(
    xintercept =
      C2_PRIMARY,
    linewidth = 0.7,
    linetype = 2
  ) +

  labs(
    title = "B. C2 POD calibration",
    x = expression(C^2),
    y = "Density",
    linetype = NULL
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    legend.position =
      "top"
  )


############################################################
## Save individual panels
############################################################

ggsave(
  file.path(
    OUT_FIG,
    "Figure_SX_POD_XtX_calibration.png"
  ),
  plot = PX,
  width = 6.5,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  file.path(
    OUT_FIG,
    "Figure_SX_POD_XtX_calibration.pdf"
  ),
  plot = PX,
  width = 6.5,
  height = 4.5,
  units = "in",
  bg = "white"
)


ggsave(
  file.path(
    OUT_FIG,
    "Figure_SX_POD_C2_calibration.png"
  ),
  plot = PC,
  width = 6.5,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  file.path(
    OUT_FIG,
    "Figure_SX_POD_C2_calibration.pdf"
  ),
  plot = PC,
  width = 6.5,
  height = 4.5,
  units = "in",
  bg = "white"
)


############################################################
## 27. PRIMARY DECISION RECORD
##
## This is deliberately written BEFORE subsequent overlap
## analyses are performed.
############################################################

DECISION_RECORD <- c(

  "Cactoblastis cactorum — BayPass POD calibration",

  "",

  paste0(
    "Date: ",
    Sys.Date()
  ),

  "",

  "Empirical BayPass SNP universe: 132283 SNPs",

  paste0(
    "Neutral POD simulated SNPs: ",
    N_POD_EXPECTED
  ),

  paste0(
    "POD SNPs retained after observed MAF >= 0.05: ",
    N_MAF_GE05
  ),

  "",

  "PRIMARY CALIBRATION DECISION:",

  "POD thresholds are calculated from simulated loci with observed pooled MAF >= 0.05, matching the ascertainment of the empirical all-SNP dataset.",

  "",

  "PRIMARY SNP-OUTLIER QUANTILE:",

  paste0(
    PRIMARY_QUANTILE
  ),

  "",

  "Sensitivity quantiles:",

  paste(
    QUANTILES[
      QUANTILES !=
        PRIMARY_QUANTILE
    ],
    collapse = ", "
  ),

  "",

  "XtX statistic used:",

  "XtXst",

  "",

  "C2 statistic used:",

  "C2",

  "",

  paste0(
    "Primary XtXst POD threshold = ",
    signif(
      XTX_PRIMARY,
      10
    )
  ),

  paste0(
    "Primary C2 POD threshold = ",
    signif(
      C2_PRIMARY,
      10
    )
  ),

  "",

  paste0(
    "Primary empirical XtX POD outlier SNPs = ",
    nrow(XTX_OUT_META)
  ),

  paste0(
    "Primary empirical C2 POD outlier SNPs = ",
    nrow(C2_OUT_META)
  ),

  "",

  "These thresholds were defined before cross-analysis SNP-overlap testing.",

  "",

  "POD-calibrated SNP outliers are complementary to, and do not replace, local-score candidate regions."
)


writeLines(
  DECISION_RECORD,
  con = file.path(
    OUT_THRESHOLD,
    "POD_threshold_decision_record_FINAL.txt"
  )
)


############################################################
## 28. SESSION INFO
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_REPRO,
    "sessionInfo_49C.txt"
  )
)


############################################################
## 29. FINAL CONSOLE REPORT
############################################################

cat(
  "\n============================================================\n"
)

cat(
  "49C COMPLETE\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "POD SNPs simulated:\n"
)

cat(
  "  ",
  N_POD_EXPECTED,
  "\n",
  sep = ""
)


cat(
  "\nPOD SNPs retained at observed MAF >= 0.05:\n"
)

cat(
  "  ",
  N_MAF_GE05,
  "\n",
  sep = ""
)


cat(
  "\nPrimary threshold definition:\n"
)

cat(
  "  99.9th percentile of MAF-filtered POD\n"
)


cat(
  "\nXtXst primary threshold:\n"
)

cat(
  "  ",
  signif(
    XTX_PRIMARY,
    10
  ),
  "\n",
  sep = ""
)


cat(
  "\nC2 primary threshold:\n"
)

cat(
  "  ",
  signif(
    C2_PRIMARY,
    10
  ),
  "\n",
  sep = ""
)


cat(
  "\nEmpirical XtX POD99.9 SNP outliers:\n"
)

cat(
  "  ",
  nrow(XTX_OUT_META),
  "\n",
  sep = ""
)


cat(
  "\nEmpirical C2 POD99.9 SNP outliers:\n"
)

cat(
  "  ",
  nrow(C2_OUT_META),
  "\n",
  sep = ""
)


cat(
  "\n------------------------------------------------------------\n"
)

cat(
  "IMPORTANT\n"
)

cat(
  "------------------------------------------------------------\n"
)

cat(
  "No XtX/C2/GEA overlap was calculated by this script.\n"
)

cat(
  "Inspect POD calibration and outlier counts before proceeding.\n"
)


cat(
  "\nMain outputs:\n\n"
)


cat(
  file.path(
    OUT_THRESHOLD,
    "Table_POD_thresholds_XtXst_C2_ALL_vs_MAFfiltered.tsv"
  ),
  "\n"
)


cat(
  file.path(
    OUT_TABLE,
    "Table_POD_empirical_outlier_counts_by_threshold.tsv"
  ),
  "\n"
)


cat(
  file.path(
    OUT_OUTLIER,
    "XtX_POD999_MAFfiltered_outlier_SNPs.tsv"
  ),
  "\n"
)


cat(
  file.path(
    OUT_OUTLIER,
    "C2_POD999_MAFfiltered_outlier_SNPs.tsv"
  ),
  "\n"
)


cat(
  "\n[49C] DONE\n"
)
