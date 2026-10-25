#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

############################################################
## 21B_build_FINAL_gt1Mb_candidate_SNPs_for_VEP.R
##
## FINAL candidate SNP universe for functional annotation
##
## XtX:
##   all SNPs inside 12 FINAL >1 Mb local-score regions
##
## C2:
##   all SNPs inside 17 FINAL >1 Mb local-score regions
##
## GEA:
##   all BF >20 SNPs for six primary covariates
##
## IMPORTANT:
## POD99.9 is NOT used to redefine candidate membership.
############################################################


# ============================================================
# PATHS
# ============================================================

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"


XTX_SCORE_FILE <- file.path(
  BASE,
  "selection_scan",
  "local_score_gt1Mb",
  "all_snps",
  "XtX",
  "XtX_local_scores_gt1Mb.tsv"
)


XTX_WINDOW_FILE <- file.path(
  BASE,
  "selection_scan",
  "local_score_gt1Mb",
  "all_snps",
  "XtX",
  "XtX_candidate_windows_gt1Mb.tsv"
)


C2_SCORE_FILE <- file.path(
  BASE,
  "selection_scan",
  "local_score_gt1Mb",
  "all_snps",
  "C2_Host",
  "C2_Host_local_scores_gt1Mb.tsv"
)


C2_WINDOW_FILE <- file.path(
  BASE,
  "selection_scan",
  "local_score_gt1Mb",
  "all_snps",
  "C2_Host",
  "C2_Host_candidate_windows_gt1Mb.tsv"
)


GEA_BF_FILE <- file.path(
  BASE,
  "selection_scan",
  "GEA_BF_median3seeds_allSNP",
  "merged_sorted_BF_median3seeds.tsv"
)


OUT_DIR <- file.path(
  BASE,
  "selection_scan",
  "VEP_FINAL_gt1Mb_selection_candidates"
)


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# SETTINGS
# ============================================================

BF_THRESHOLD <- 20


PRIMARY_COVARIATES <- c(
  "TempPC1",
  "PrecipPC1",
  "AridityIndex",
  "Elevation",
  "Latitude",
  "Longitude"
)


# ============================================================
# HELPERS
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


clean_scaffold <- function(x) {

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


read_windows <- function(
  filename,
  expected_n
) {

  x <- fread(
    filename,
    showProgress = FALSE
  )


  required <- c(
    "chr",
    "beg",
    "end"
  )


  miss <- setdiff(
    required,
    names(x)
  )


  if (
    length(miss) > 0L
  ) {

    stop(
      filename,
      " missing columns: ",
      paste(
        miss,
        collapse = ", "
      )
    )
  }


  x[
    ,
    `:=`(
      chr =
        as.character(chr),

      Scaffold_norm =
        clean_scaffold(chr),

      beg =
        as.integer(beg),

      end =
        as.integer(end)
    )
  ]


  if (
    nrow(x) != expected_n
  ) {

    stop(
      "Expected ",
      expected_n,
      " windows in ",
      filename,
      "; found ",
      nrow(x)
    )
  }


  x
}


extract_window_snps <- function(
  score_file,
  windows,
  analysis_name
) {

  s <- fread(
    score_file,
    showProgress = FALSE
  )


  required <- c(
    "chr",
    "pos"
  )


  miss <- setdiff(
    required,
    names(s)
  )


  if (
    length(miss) > 0L
  ) {

    stop(
      score_file,
      " missing columns: ",
      paste(
        miss,
        collapse = ", "
      )
    )
  }


  s[
    ,
    `:=`(
      chr =
        as.character(chr),

      Scaffold_norm =
        clean_scaffold(chr),

      pos =
        as.integer(pos)
    )
  ]


  keep <- vector(
    "list",
    nrow(windows)
  )


  for (
    i in seq_len(
      nrow(windows)
    )
  ) {

    keep[[i]] <- s[

      Scaffold_norm ==
        windows$Scaffold_norm[i] &

      pos >=
        windows$beg[i] &

      pos <=
        windows$end[i],

      .(
        CHROM = chr,
        POS = pos
      )
    ]
  }


  out <- unique(
    rbindlist(
      keep,
      use.names = TRUE
    )
  )


  out[
    ,
    Analysis :=
      analysis_name
  ]


  out
}


# ============================================================
# CHECK FILES
# ============================================================

for (
  f in c(
    XTX_SCORE_FILE,
    XTX_WINDOW_FILE,
    C2_SCORE_FILE,
    C2_WINDOW_FILE,
    GEA_BF_FILE
  )
) {

  check_file(f)
}


# ============================================================
# FINAL XtX / C2 WINDOWS
# ============================================================

xtx_windows <- read_windows(
  XTX_WINDOW_FILE,
  12L
)


c2_windows <- read_windows(
  C2_WINDOW_FILE,
  17L
)


# ============================================================
# ALL SNPs INSIDE FINAL REGIONS
# ============================================================

xtx_snps <- extract_window_snps(
  XTX_SCORE_FILE,
  xtx_windows,
  "XtX"
)


c2_snps <- extract_window_snps(
  C2_SCORE_FILE,
  c2_windows,
  "C2"
)


# ============================================================
# GEA BF >20 SNPs — UNCHANGED
# ============================================================

gea <- fread(
  GEA_BF_FILE,
  showProgress = FALSE
)


required_gea <- c(
  "chr",
  "pos",
  "Covariate",
  "BFdb"
)


miss <- setdiff(
  required_gea,
  names(gea)
)


if (
  length(miss) > 0L
) {

  stop(
    "GEA table missing columns: ",
    paste(
      miss,
      collapse = ", "
    )
  )
}


gea[
  ,
  `:=`(
    chr =
      as.character(chr),

    pos =
      as.integer(pos),

    Covariate =
      as.character(Covariate),

    BFdb =
      as.numeric(BFdb)
  )
]


gea <- gea[
  Covariate %in%
    PRIMARY_COVARIATES &
    BFdb > BF_THRESHOLD
]


gea_snps <- gea[
  ,
  .(
    GEA_covariates =
      paste(
        sort(
          unique(Covariate)
        ),
        collapse = ";"
      ),

    GEA_max_BFdb =
      max(
        BFdb,
        na.rm = TRUE
      )
  ),
  by = .(
    CHROM = chr,
    POS = pos
  )
]


# ============================================================
# MASTER POSITION SET
# ============================================================

xtx_key <- unique(
  xtx_snps[
    ,
    .(
      CHROM,
      POS
    )
  ]
)


xtx_key[
  ,
  In_XtX :=
    TRUE
]


c2_key <- unique(
  c2_snps[
    ,
    .(
      CHROM,
      POS
    )
  ]
)


c2_key[
  ,
  In_C2 :=
    TRUE
]


all_positions <- unique(

  rbindlist(

    list(

      xtx_key[
        ,
        .(
          CHROM,
          POS
        )
      ],

      c2_key[
        ,
        .(
          CHROM,
          POS
        )
      ],

      gea_snps[
        ,
        .(
          CHROM,
          POS
        )
      ]

    ),

    use.names = TRUE

  )
)


# ============================================================
# ADD MEMBERSHIP
# ============================================================

all_positions <- merge(
  all_positions,
  xtx_key,
  by = c(
    "CHROM",
    "POS"
  ),
  all.x = TRUE
)


all_positions <- merge(
  all_positions,
  c2_key,
  by = c(
    "CHROM",
    "POS"
  ),
  all.x = TRUE
)


all_positions <- merge(
  all_positions,
  gea_snps,
  by = c(
    "CHROM",
    "POS"
  ),
  all.x = TRUE
)


all_positions[
  is.na(In_XtX),
  In_XtX :=
    FALSE
]


all_positions[
  is.na(In_C2),
  In_C2 :=
    FALSE
]


all_positions[
  ,
  In_GEA :=
    !is.na(
      GEA_covariates
    )
]


# ============================================================
# ANALYSIS COMBINATIONS
# ============================================================

all_positions[
  ,
  N_analyses :=
    as.integer(In_XtX) +
    as.integer(In_C2) +
    as.integer(In_GEA)
]


all_positions[
  ,
  Analysis_combination := fcase(

    In_XtX &
      !In_C2 &
      !In_GEA,
    "XtX",

    !In_XtX &
      In_C2 &
      !In_GEA,
    "C2",

    !In_XtX &
      !In_C2 &
      In_GEA,
    "GEA",

    In_XtX &
      In_C2 &
      !In_GEA,
    "XtX+C2",

    In_XtX &
      !In_C2 &
      In_GEA,
    "XtX+GEA",

    !In_XtX &
      In_C2 &
      In_GEA,
    "C2+GEA",

    In_XtX &
      In_C2 &
      In_GEA,
    "XtX+C2+GEA",

    default =
      "ERROR"
  )
]


if (
  any(
    all_positions$Analysis_combination ==
      "ERROR"
  )
) {

  stop(
    "Unclassified analysis combination."
  )
}


setorder(
  all_positions,
  CHROM,
  POS
)


# ============================================================
# QC
# ============================================================

cat(
  "\n=============================================\n"
)

cat(
  "FINAL >1 Mb CANDIDATE SNP UNION\n"
)

cat(
  "=============================================\n"
)


cat(
  "XtX candidate-region SNPs: ",
  nrow(xtx_key),
  "\n",
  sep = ""
)


cat(
  "C2 candidate-region SNPs: ",
  nrow(c2_key),
  "\n",
  sep = ""
)


cat(
  "GEA BF>20 unique SNPs: ",
  nrow(gea_snps),
  "\n",
  sep = ""
)


cat(
  "FINAL non-redundant union: ",
  nrow(all_positions),
  "\n",
  sep = ""
)


cat(
  "\nAnalysis combinations:\n"
)


print(
  all_positions[
    ,
    .N,
    by =
      Analysis_combination
  ][
    order(
      -N
    )
  ]
)


# ============================================================
# WRITE
# ============================================================

MASTER_OUT <- file.path(
  OUT_DIR,
  "ALL_selection_candidate_SNPs_master_gt1Mb_FINAL.tsv"
)


fwrite(
  all_positions,
  MASTER_OUT,
  sep = "\t",
  quote = FALSE,
  na = ""
)


REGIONS_OUT <- file.path(
  OUT_DIR,
  "ALL_selection_candidate_SNP_positions_gt1Mb_FINAL.tsv"
)


fwrite(

  all_positions[
    ,
    .(
      CHROM,
      FROM = POS,
      TO = POS
    )
  ],

  REGIONS_OUT,

  sep = "\t",

  col.names = FALSE,

  quote = FALSE

)


cat(
  "\nMaster table:\n",
  MASTER_OUT,
  "\n",
  sep = ""
)


cat(
  "\nDONE\n"
)
