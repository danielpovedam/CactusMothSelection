#!/usr/bin/env Rscript

############################################################
## 41C_gt1Mb_prepare_and_local_score_C2_geographic_sensitivity.R
##
## Cactoblastis cactorum
##
## FINAL geographically restricted C2 sensitivity analysis
##
## The BayPass restricted analysis itself is unchanged.
##
## Final local-score eligibility:
##   - scaffold length > 1 Mb
##   - min.maf = 0.05
##   - min.nsnp = 20
##   - xi = 2
##   - alpha = 0.01
##
## This matches the FINAL primary C2 local-score analysis.
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 0. SETTINGS
# ============================================================

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

POPBASE <- "/home/dpoveda/micetral/Population_snps"

ROOT <- file.path(
  BASE,
  "sensitivity_C2_geographic"
)

REP_SEED <- "12345"

RUNDIR <- file.path(
  ROOT,
  "runs",
  paste0("seed_", REP_SEED)
)

META_FILE <- file.path(
  ROOT,
  "input",
  "restricted_C2_polymorphic_SNP_metadata.tsv"
)

FAI_FILE <- file.path(
  POPBASE,
  "Ccactorum.fcs_gx_2026_FINAL.fa.fai"
)

OUTBASE <- file.path(
  ROOT,
  "local_score_gt1Mb",
  "C2_Host"
)

UTILS <- paste0(
  "/home/dpoveda/micetral/selection/baypass/",
  "baypass_public/utils/baypass_utils.R"
)

EXPECTED_SNPS <- 130732L

MIN_SCAFFOLD_BP <- 1000000L
MIN_MAF <- 0.05
MIN_NSNP <- 20L
XI_C2 <- 2
ALPHA_C2 <- 0.01


# ============================================================
# 1. HELPERS
# ============================================================

check_file <- function(f) {

  if (
    !file.exists(f) ||
    is.na(file.info(f)$size) ||
    file.info(f)$size <= 0
  ) {

    stop(
      "Missing or empty file:\n",
      f
    )
  }
}


normalize_scaffold <- function(x) {

  x <- as.character(x)

  x <- sub(
    "^lcl\\|",
    "",
    x
  )

  x <- sub(
    "l_1_1$",
    "",
    x
  )

  x <- sub(
    "l_1$",
    "",
    x
  )

  x
}


find_one <- function(
  pattern,
  label
) {

  ff <- list.files(
    RUNDIR,
    pattern = pattern,
    full.names = TRUE
  )

  ff <- ff[
    file.exists(ff) &
      file.info(ff)$size > 0
  ]

  if (
    length(ff) == 0L
  ) {

    stop(
      "\nCould not find ",
      label,
      " in:\n",
      RUNDIR,
      "\nPattern: ",
      pattern
    )
  }

  if (
    length(ff) > 1L
  ) {

    stop(
      "\nMore than one candidate ",
      label,
      " file found:\n",
      paste(
        ff,
        collapse = "\n"
      )
    )
  }

  ff[1L]
}


pick_column <- function(
  dt,
  preferred,
  regex = NULL,
  label
) {

  nms <- names(dt)

  for (
    x in preferred
  ) {

    if (
      x %in% nms
    ) {

      return(x)
    }
  }

  if (
    !is.null(regex)
  ) {

    hits <- grep(
      regex,
      nms,
      value = TRUE,
      ignore.case = TRUE
    )

    if (
      length(hits) == 1L
    ) {

      return(hits)
    }
  }

  stop(
    "Could not uniquely identify ",
    label,
    " column.\nObserved columns: ",
    paste(
      nms,
      collapse = ", "
    )
  )
}


for (
  f in c(
    UTILS,
    META_FILE,
    FAI_FILE
  )
) {

  check_file(f)
}


if (
  !dir.exists(RUNDIR)
) {

  stop(
    "Representative BayPass run directory not found:\n",
    RUNDIR
  )
}


dir.create(
  OUTBASE,
  recursive = TRUE,
  showWarnings = FALSE
)


source(UTILS)


if (
  !exists("compute.local.scores")
) {

  stop(
    "compute.local.scores() not found."
  )
}


cat(
  "\n============================================\n"
)

cat(
  "FINAL GEOGRAPHIC C2 SENSITIVITY: >1 Mb\n"
)

cat(
  "============================================\n"
)


# ============================================================
# 2. BAYPASS OUTPUT FILES
# ============================================================

C2_FILE <- find_one(
  "summary_contrast.*\\.out$",
  "BayPass summary_contrast"
)

PI_FILE <- find_one(
  "summary_pi_xtx.*\\.out$",
  "BayPass summary_pi_xtx"
)


cat(
  "C2 file: ",
  basename(C2_FILE),
  "\n",
  sep = ""
)

cat(
  "M_P file: ",
  basename(PI_FILE),
  "\n",
  sep = ""
)


# ============================================================
# 3. TRUE SCAFFOLD LENGTHS
# ============================================================

fai <- fread(
  FAI_FILE,
  header = FALSE,
  showProgress = FALSE
)


if (
  ncol(fai) < 2L
) {

  stop(
    "FAI has fewer than two columns."
  )
}


setnames(
  fai,
  1:2,
  c(
    "FASTA_scaffold",
    "Scaffold_length_bp"
  )
)


fai <- fai[
  ,
  .(
    FASTA_scaffold =
      as.character(FASTA_scaffold),

    Scaffold_length_bp =
      as.integer(Scaffold_length_bp)
  )
]


fai[
  ,
  Scaffold_norm :=
    normalize_scaffold(
      FASTA_scaffold
    )
]


fai[
  ,
  Eligible_gt1Mb :=
    Scaffold_length_bp >
    MIN_SCAFFOLD_BP
]


if (
  anyDuplicated(
    fai$Scaffold_norm
  )
) {

  stop(
    "Duplicated normalized scaffold names in FAI."
  )
}


cat(
  "FASTA scaffolds: ",
  nrow(fai),
  "\n",
  sep = ""
)

cat(
  "FASTA scaffolds >1 Mb: ",
  sum(fai$Eligible_gt1Mb),
  "\n",
  sep = ""
)


# ============================================================
# 4. RESTRICTED SNP METADATA
# ============================================================

meta <- fread(
  META_FILE,
  showProgress = FALSE
)


required_meta <- c(
  "BayPass_row",
  "CHROM",
  "SNP",
  "POS",
  "REF",
  "ALT"
)


missing_meta <- setdiff(
  required_meta,
  names(meta)
)


if (
  length(missing_meta) > 0L
) {

  stop(
    "Restricted metadata missing columns: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}


if (
  nrow(meta) != EXPECTED_SNPS
) {

  stop(
    "Expected ",
    EXPECTED_SNPS,
    " metadata rows; found ",
    nrow(meta)
  )
}


meta[
  ,
  Restricted_MRK :=
    seq_len(.N)
]


meta[
  ,
  `:=`(
    CHROM =
      as.character(CHROM),

    POS =
      as.integer(POS),

    SNP =
      as.character(SNP),

    REF =
      as.character(REF),

    ALT =
      as.character(ALT),

    BayPass_row =
      as.integer(BayPass_row),

    Scaffold_norm =
      normalize_scaffold(CHROM)
  )
]


# ============================================================
# 5. READ BAYPASS C2
# ============================================================

c2_raw <- fread(
  C2_FILE,
  header = TRUE,
  showProgress = FALSE
)


if (
  nrow(c2_raw) != EXPECTED_SNPS
) {

  stop(
    "C2 output contains ",
    nrow(c2_raw),
    " rows; expected ",
    EXPECTED_SNPS
  )
}


C2_COL <- pick_column(
  c2_raw,
  preferred = c(
    "C2",
    "C2_1",
    "C2contrast1"
  ),
  regex = "^C2($|[_\\.])",
  label = "C2 statistic"
)


CONTRAST_COL <- NULL


if (
  "CONTRAST" %in%
    names(c2_raw)
) {

  CONTRAST_COL <- "CONTRAST"

} else {

  tmp <- grep(
    "contrast",
    names(c2_raw),
    value = TRUE,
    ignore.case = TRUE
  )

  tmp <- setdiff(
    tmp,
    C2_COL
  )

  if (
    length(tmp) == 1L
  ) {

    CONTRAST_COL <- tmp
  }
}


MRK_C2_COL <- if (
  "MRK" %in%
    names(c2_raw)
) {

  "MRK"

} else {

  NULL
}


c2 <- data.table(

  Restricted_MRK =
    if (
      !is.null(MRK_C2_COL)
    ) {

      as.integer(
        c2_raw[[MRK_C2_COL]]
      )

    } else {

      seq_len(
        nrow(c2_raw)
      )
    },

  C2 =
    as.numeric(
      c2_raw[[C2_COL]]
    )
)


if (
  !is.null(CONTRAST_COL)
) {

  c2[
    ,
    CONTRAST :=
      as.numeric(
        c2_raw[[CONTRAST_COL]]
      )
  ]

} else {

  c2[
    ,
    CONTRAST :=
      1
  ]
}


# ============================================================
# 6. READ M_P
# ============================================================

pi_raw <- fread(
  PI_FILE,
  header = TRUE,
  showProgress = FALSE
)


if (
  nrow(pi_raw) != EXPECTED_SNPS
) {

  stop(
    "summary_pi_xtx contains ",
    nrow(pi_raw),
    " rows; expected ",
    EXPECTED_SNPS
  )
}


MP_COL <- pick_column(
  pi_raw,
  preferred = c(
    "M_P",
    "MP",
    "Mean_P"
  ),
  regex = "^M[_\\.]?P$",
  label = "M_P"
)


MRK_PI_COL <- if (
  "MRK" %in%
    names(pi_raw)
) {

  "MRK"

} else {

  NULL
}


pi_dt <- data.table(

  Restricted_MRK =
    if (
      !is.null(MRK_PI_COL)
    ) {

      as.integer(
        pi_raw[[MRK_PI_COL]]
      )

    } else {

      seq_len(
        nrow(pi_raw)
      )
    },

  M_P =
    as.numeric(
      pi_raw[[MP_COL]]
    )
)


# ============================================================
# 7. BUILD COMPLETE 130,732-SNP PREPARED TABLE
# ============================================================

dt <- merge(
  meta,
  pi_dt,
  by = "Restricted_MRK",
  all.x = TRUE,
  sort = FALSE
)


dt <- merge(
  dt,
  c2,
  by = "Restricted_MRK",
  all.x = TRUE,
  sort = FALSE
)


setorder(
  dt,
  Restricted_MRK
)


if (
  nrow(dt) != EXPECTED_SNPS ||
  anyNA(dt$M_P) ||
  anyNA(dt$C2)
) {

  stop(
    "Failed to join metadata, M_P and C2."
  )
}


# ============================================================
# 8. C2 -> -log10(P)
# ============================================================

dt[
  ,
  p_value :=
    pchisq(
      C2,
      df = 1,
      lower.tail = FALSE
    )
]


dt[
  ,
  log10_p :=
    -log10(
      pmax(
        p_value,
        .Machine$double.xmin
      )
    )
]


# ============================================================
# 9. ATTACH TRUE SCAFFOLD LENGTHS
# ============================================================

dt <- merge(

  dt,

  fai[
    ,
    .(
      Scaffold_norm,
      Scaffold_length_bp,
      Eligible_gt1Mb
    )
  ],

  by = "Scaffold_norm",

  all.x = TRUE,

  sort = FALSE
)


if (
  anyNA(
    dt$Scaffold_length_bp
  )
) {

  bad <- unique(
    dt[
      is.na(Scaffold_length_bp),
      CHROM
    ]
  )

  stop(
    "Restricted SNP scaffolds failed to map to FAI:\n",
    paste(
      bad,
      collapse = "\n"
    )
  )
}


# ============================================================
# 10. SAVE COMPLETE RESTRICTED BAYPASS TABLE
#
# Keep this because the genome-wide primary-vs-restricted
# C2 correlation still legitimately uses all 130,732
# restricted polymorphic SNPs.
# ============================================================

prepared_full <- dt[
  ,
  .(
    Restricted_MRK,
    BayPass_row,
    CHROM,
    SNP,
    POS,
    REF,
    ALT,
    M_P,
    CONTRAST,
    C2,
    p_value,
    log10_p,
    Scaffold_norm,
    Scaffold_length_bp,
    Eligible_gt1Mb
  )
]


PREPARED_FULL_FILE <- file.path(
  OUTBASE,
  "C2_Host_geographic16_seed12345_prepared_ALL130732.tsv"
)


fwrite(
  prepared_full,
  PREPARED_FULL_FILE,
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 11. >1 Mb LOCAL-SCORE INPUT
# ============================================================

prepared_gt1Mb <- prepared_full[
  Eligible_gt1Mb == TRUE
]


setorder(
  prepared_gt1Mb,
  CHROM,
  POS
)


prepared_gt1Mb[
  ,
  MAF_proxy :=
    pmin(
      M_P,
      1 - M_P
    )
]


scaffold_qc <- prepared_gt1Mb[
  MAF_proxy >= MIN_MAF,
  .(
    N_MAF_eligible_SNPs =
      .N
  ),
  by = .(
    CHROM,
    Scaffold_norm,
    Scaffold_length_bp
  )
]


scaffold_qc[
  ,
  Eligible_gt20SNP :=
    N_MAF_eligible_SNPs >
    MIN_NSNP
]


INPUT_GT1MB_FILE <- file.path(
  OUTBASE,
  "C2_Host_geographic16_local_score_input_gt1Mb.tsv"
)


fwrite(
  prepared_gt1Mb,
  INPUT_GT1MB_FILE,
  sep = "\t",
  quote = FALSE
)


fwrite(
  scaffold_qc,
  file.path(
    OUTBASE,
    "C2_Host_geographic16_scaffold_QC_gt1Mb.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 12. QC BEFORE LOCAL SCORE
# ============================================================

cat(
  "\nRestricted BayPass SNPs: ",
  nrow(prepared_full),
  "\n",
  sep = ""
)


cat(
  "SNPs on scaffolds >1 Mb: ",
  nrow(prepared_gt1Mb),
  "\n",
  sep = ""
)


cat(
  "Percent retained by >1 Mb criterion: ",
  round(
    100 *
      nrow(prepared_gt1Mb) /
      nrow(prepared_full),
    4
  ),
  "%\n",
  sep = ""
)


cat(
  "Scaffolds >1 Mb represented: ",
  uniqueN(
    prepared_gt1Mb$Scaffold_norm
  ),
  "\n",
  sep = ""
)


cat(
  "Scaffolds >1 Mb with >20 MAF-eligible SNPs: ",
  sum(
    scaffold_qc$Eligible_gt20SNP
  ),
  "\n",
  sep = ""
)


# ============================================================
# 13. FINAL RESTRICTED LOCAL SCORE
# ============================================================

res <- compute.local.scores(

  snp.position = data.frame(
    chr =
      prepared_gt1Mb$CHROM,
    pos =
      prepared_gt1Mb$POS
  ),

  snp.pi =
    prepared_gt1Mb$M_P,

  snp.pvalue =
    prepared_gt1Mb$log10_p,

  xi =
    XI_C2,

  pval.local.score.thres =
    ALPHA_C2,

  min.maf =
    MIN_MAF,

  min.nsnp =
    MIN_NSNP,

  plot.pvalhist =
    FALSE,

  manplot =
    FALSE
)


# ============================================================
# 14. LOCAL-SCORE SNP OUTPUT
# ============================================================

scores <- as.data.table(
  res$res.local.scores
)


if (
  nrow(scores) > 0L
) {

  scores[
    ,
    `:=`(
      Dataset =
        "geographic_16pop_gt1Mb",

      Analysis =
        "C2_Host_sensitivity"
    )
  ]


  setcolorder(
    scores,
    c(
      "Dataset",
      "Analysis",
      setdiff(
        names(scores),
        c(
          "Dataset",
          "Analysis"
        )
      )
    )
  )
}


SCORE_FILE <- file.path(
  OUTBASE,
  "C2_Host_geographic16_local_scores_gt1Mb.tsv"
)


fwrite(
  scores,
  SCORE_FILE,
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 15. FINAL RESTRICTED CANDIDATE WINDOWS
# ============================================================

windows <- as.data.table(
  res$significant.windows
)


if (
  nrow(windows) > 0L
) {

  windows[
    ,
    `:=`(
      Dataset =
        "geographic_16pop_gt1Mb",

      Analysis =
        "C2_Host_sensitivity"
    )
  ]


  setcolorder(
    windows,
    c(
      "Dataset",
      "Analysis",
      setdiff(
        names(windows),
        c(
          "Dataset",
          "Analysis"
        )
      )
    )
  )
}


WINDOW_FILE <- file.path(
  OUTBASE,
  "C2_Host_geographic16_candidate_windows_gt1Mb.tsv"
)


fwrite(
  windows,
  WINDOW_FILE,
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 16. FINAL QC SUMMARY
# ============================================================

settings <- data.table(

  Parameter = c(
    "Representative_seed",
    "Restricted_BayPass_SNPs",
    "Minimum_scaffold_length_bp",
    "N_SNPs_on_gt1Mb_scaffolds",
    "Percent_SNPs_retained_gt1Mb",
    "min_maf_local_score",
    "min_nsnp",
    "xi_C2",
    "alpha_C2",
    "N_gt1Mb_scaffolds",
    "N_gt1Mb_scaffolds_gt20_MAFeligible_SNPs",
    "N_local_score_SNPs",
    "N_candidate_windows"
  ),

  Value = c(
    REP_SEED,
    nrow(prepared_full),
    MIN_SCAFFOLD_BP,
    nrow(prepared_gt1Mb),
    100 *
      nrow(prepared_gt1Mb) /
      nrow(prepared_full),
    MIN_MAF,
    MIN_NSNP,
    XI_C2,
    ALPHA_C2,
    uniqueN(
      prepared_gt1Mb$Scaffold_norm
    ),
    sum(
      scaffold_qc$Eligible_gt20SNP
    ),
    nrow(scores),
    nrow(windows)
  )
)


SETTINGS_FILE <- file.path(
  OUTBASE,
  "C2_Host_geographic16_local_score_gt1Mb_QC.tsv"
)


fwrite(
  settings,
  SETTINGS_FILE,
  sep = "\t",
  quote = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "FINAL RESTRICTED C2 LOCAL-SCORE SUMMARY\n"
)

cat(
  "============================================\n"
)


print(settings)


cat(
  "\nCandidate windows: ",
  nrow(windows),
  "\n",
  sep = ""
)


cat(
  "\nOutput:\n",
  WINDOW_FILE,
  "\n",
  sep = ""
)


cat(
  "\nDONE\n"
)
