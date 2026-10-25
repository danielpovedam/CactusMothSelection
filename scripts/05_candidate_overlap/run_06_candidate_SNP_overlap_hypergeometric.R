#!/usr/bin/env Rscript

############################################################
## 47_FINAL_gt1Mb_POD_exactSNP_hypergeometric_v2.R
##
## Cactoblastis cactorum
##
## FINAL >1 Mb downstream update
##
## IMPORTANT:
##
## POD calibration is NOT recalculated here.
##
## We use the frozen empirical POD99.9 outlier files from
## the finalized 49C analysis:
##
##   XtXst > POD99.9
##   C2    > POD99.9
##
## We only:
##
## 1. restrict POD99.9 SNPs to the revised local-score
##    candidate regions:
##       XtX = 12 regions
##       C2  = 17 regions
##
## 2. define the revised exact-SNP universe from the
##    actual >1 Mb + >20-SNP local-score analysis
##
## 3. recalculate:
##       XtX vs C2 exact SNP overlap
##       XtX vs six GEA covariates
##       C2  vs six GEA covariates
##       pooled XtX vs GEA
##       pooled C2 vs GEA
##
## 4. perform one-sided hypergeometric tests
##
## 5. BH-adjust the 12 covariate-specific tests
############################################################


suppressPackageStartupMessages({
  library(data.table)
})


############################################################
## 0. PATHS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"


############################################################
## FINAL >1 Mb local-score regions
############################################################

XTX_REGIONS_FILE <- file.path(
  BASE,
  "selection_scan",
  "FINAL_gt1Mb",
  "Table_XtX_candidate_regions_allSNP_gt1Mb_FINAL.tsv"
)

C2_REGIONS_FILE <- file.path(
  BASE,
  "selection_scan",
  "FINAL_gt1Mb",
  "Table_C2_Host_candidate_regions_allSNP_gt1Mb_FINAL.tsv"
)


############################################################
## Frozen POD99.9 empirical outliers from 49C
############################################################

XTX_POD_FILE <- file.path(
  BASE,
  "POD_analysis",
  "04_empirical_outliers",
  "XtX_POD999_MAFfiltered_outlier_SNPs.tsv"
)

C2_POD_FILE <- file.path(
  BASE,
  "POD_analysis",
  "04_empirical_outliers",
  "C2_POD999_MAFfiltered_outlier_SNPs.tsv"
)


############################################################
## Revised local-score universe
##
## This contains the 123,428 SNP records actually retained
## after:
##
## scaffold >1 Mb
## AND
## min.nsnp >20
############################################################

UNIVERSE_FILE <- file.path(
  BASE,
  "selection_scan",
  "local_score_gt1Mb",
  "all_snps",
  "XtX",
  "XtX_local_scores_gt1Mb.tsv"
)


############################################################
## Authoritative 207 BF>20 GEA SNP-covariate associations
############################################################

GEA_FILE <- file.path(
  BASE,
  "selection_scan",
  "GEA_SUPPLEMENTARY_PUBLICATION_FINAL",
  "Table_Sy_ALL_207_BF20_SNP_covariate_associations_FINAL.tsv"
)


############################################################
## Output
############################################################

OUTDIR <- file.path(
  BASE,
  "selection_scan",
  "FINAL_gt1Mb",
  "POD_exactSNP_hypergeometric_v2"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 1. SETTINGS
############################################################

COVARIATE_ORDER <- c(
  "TempPC1",
  "PrecipPC1",
  "AridityIndex",
  "Elevation",
  "Latitude",
  "Longitude"
)


############################################################
## 2. HELPERS
############################################################

check_file <- function(f) {

  if (
    !file.exists(f) ||
    is.na(file.info(f)$size) ||
    file.info(f)$size <= 0
  ) {

    stop(
      "Missing or empty input file:\n",
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


make_snp_id <- function(
  scaffold,
  position
) {

  paste0(
    normalize_scaffold(scaffold),
    ":",
    as.integer(position)
  )
}


find_col <- function(
  dt,
  candidates,
  required = TRUE
) {

  hit <- candidates[
    candidates %in%
      names(dt)
  ]

  if (
    length(hit) > 0L
  ) {

    return(
      hit[1]
    )
  }


  if (
    required
  ) {

    stop(
      "\nCould not find any expected column:\n",
      paste(
        candidates,
        collapse = ", "
      ),
      "\n\nAvailable columns:\n",
      paste(
        names(dt),
        collapse = ", "
      )
    )
  }


  NA_character_
}


############################################################
## One-sided enrichment test
############################################################

hyper_test <- function(
  setA,
  setB,
  universe,
  comparison
) {

  universe <- unique(
    universe
  )

  setA <- unique(
    intersect(
      setA,
      universe
    )
  )

  setB <- unique(
    intersect(
      setB,
      universe
    )
  )


  N <- length(
    universe
  )

  A <- length(
    setA
  )

  B <- length(
    setB
  )

  O <- length(
    intersect(
      setA,
      setB
    )
  )


  E <- A * B / N


  fold <- ifelse(
    E > 0,
    O / E,
    NA_real_
  )


  p <- phyper(
    q = O - 1,
    m = A,
    n = N - A,
    k = B,
    lower.tail = FALSE
  )


  data.table(
    Comparison = comparison,
    SNP_universe = N,
    Candidate_SNPs_A = A,
    Candidate_SNPs_B = B,
    Observed_shared_SNPs = O,
    Expected_shared_SNPs = E,
    Fold_enrichment = fold,
    P_hypergeometric = p
  )
}


############################################################
## 3. CHECK FILES
############################################################

for (
  f in c(
    XTX_REGIONS_FILE,
    C2_REGIONS_FILE,
    XTX_POD_FILE,
    C2_POD_FILE,
    UNIVERSE_FILE,
    GEA_FILE
  )
) {

  check_file(f)
}


############################################################
## 4. READ FINAL >1 Mb REGIONS
############################################################

xtx_regions <- fread(
  XTX_REGIONS_FILE,
  showProgress = FALSE
)

c2_regions <- fread(
  C2_REGIONS_FILE,
  showProgress = FALSE
)


if (
  nrow(xtx_regions) != 12L
) {

  stop(
    "Expected 12 XtX regions; found ",
    nrow(xtx_regions)
  )
}


if (
  nrow(c2_regions) != 17L
) {

  stop(
    "Expected 17 C2 regions; found ",
    nrow(c2_regions)
  )
}


xtx_regions[
  ,
  `:=`(
    Scaffold_norm =
      normalize_scaffold(
        Scaffold_accession
      ),

    Window_start =
      as.integer(
        Window_start
      ),

    Window_end =
      as.integer(
        Window_end
      )
  )
]


c2_regions[
  ,
  `:=`(
    Scaffold_norm =
      normalize_scaffold(
        Scaffold_accession
      ),

    Window_start =
      as.integer(
        Window_start
      ),

    Window_end =
      as.integer(
        Window_end
      )
  )
]


############################################################
## 5. BUILD EXACT REVISED SNP UNIVERSE
############################################################

cat(
  "\n============================================\n"
)

cat(
  "FINAL >1 Mb SNP UNIVERSE\n"
)

cat(
  "============================================\n"
)


univ <- fread(
  UNIVERSE_FILE,
  showProgress = FALSE
)


UNIV_CHR <- find_col(
  univ,
  c(
    "chr",
    "CHROM"
  )
)


UNIV_POS <- find_col(
  univ,
  c(
    "pos",
    "POS"
  )
)


univ[
  ,
  SNP_ID :=
    make_snp_id(
      get(UNIV_CHR),
      get(UNIV_POS)
    )
]


UNIVERSE <- unique(
  univ$SNP_ID
)


cat(
  "SNP records retained by local score: ",
  nrow(univ),
  "\n",
  sep = ""
)


cat(
  "Unique eligible SNP positions: ",
  length(UNIVERSE),
  "\n",
  sep = ""
)


cat(
  "Duplicate positions: ",
  nrow(univ) -
    length(UNIVERSE),
  "\n",
  sep = ""
)


############################################################
## 6. READ FROZEN POD99.9 OUTLIERS
############################################################

xtx_pod <- fread(
  XTX_POD_FILE,
  showProgress = FALSE
)

c2_pod <- fread(
  C2_POD_FILE,
  showProgress = FALSE
)


############################################################
## Mandatory columns from authoritative 49C outputs
############################################################

required_xtx <- c(
  "MRK",
  "XtXst",
  "CHROM",
  "POS"
)


required_c2 <- c(
  "MRK",
  "C2",
  "CHROM",
  "POS"
)


if (
  !all(
    required_xtx %in%
      names(xtx_pod)
  )
) {

  stop(
    "XtX POD outlier file missing expected columns.\n",
    paste(
      names(xtx_pod),
      collapse = ", "
    )
  )
}


if (
  !all(
    required_c2 %in%
      names(c2_pod)
  )
) {

  stop(
    "C2 POD outlier file missing expected columns.\n",
    paste(
      names(c2_pod),
      collapse = ", "
    )
  )
}


xtx_pod[
  ,
  `:=`(
    MRK =
      as.integer(MRK),

    POS =
      as.integer(POS),

    Scaffold_norm =
      normalize_scaffold(CHROM),

    SNP_ID =
      make_snp_id(
        CHROM,
        POS
      )
  )
]


c2_pod[
  ,
  `:=`(
    MRK =
      as.integer(MRK),

    POS =
      as.integer(POS),

    Scaffold_norm =
      normalize_scaffold(CHROM),

    SNP_ID =
      make_snp_id(
        CHROM,
        POS
      )
  )
]


cat(
  "\nFrozen genome-wide POD99.9 outliers:\n"
)

cat(
  "XtX = ",
  nrow(xtx_pod),
  "\n",
  sep = ""
)

cat(
  "C2  = ",
  nrow(c2_pod),
  "\n",
  sep = ""
)


############################################################
## 7. ASSIGN POD OUTLIERS TO FINAL REGIONS
############################################################

assign_pod_to_regions <- function(
  pod,
  regions,
  analysis
) {

  result <- list()

  counter <- 1L


  for (
    i in seq_len(
      nrow(regions)
    )
  ) {

    hits <- pod[

      Scaffold_norm ==
        regions$Scaffold_norm[i] &

      POS >=
        regions$Window_start[i] &

      POS <=
        regions$Window_end[i]

    ]


    if (
      nrow(hits) > 0L
    ) {

      x <- copy(
        hits
      )


      x[
        ,
        `:=`(
          Analysis =
            analysis,

          Region_ID =
            regions$Region_ID[i],

          Region_start =
            regions$Window_start[i],

          Region_end =
            regions$Window_end[i],

          Gene_ID =
            regions$Gene_ID[i],

          Gene_symbol =
            regions$Gene_symbol[i]
        )
      ]


      result[[counter]] <-
        x


      counter <- counter + 1L
    }
  }


  if (
    length(result) == 0L
  ) {

    return(
      data.table()
    )
  }


  rbindlist(
    result,
    fill = TRUE
  )
}


xtx_inside <- assign_pod_to_regions(
  xtx_pod,
  xtx_regions,
  "XtX"
)


c2_inside <- assign_pod_to_regions(
  c2_pod,
  c2_regions,
  "C2"
)


############################################################
## 8. UNIQUE POD99.9 SNP SETS INSIDE FINAL REGIONS
############################################################

xtx_unique <- unique(
  xtx_inside[
    ,
    .(
      SNP_ID,
      MRK,
      CHROM,
      Scaffold_norm,
      POS,
      SNP,
      XtXst
    )
  ],
  by = "SNP_ID"
)


c2_unique <- unique(
  c2_inside[
    ,
    .(
      SNP_ID,
      MRK,
      CHROM,
      Scaffold_norm,
      POS,
      SNP,
      C2
    )
  ],
  by = "SNP_ID"
)


############################################################
## Restrict to revised exact-SNP universe
############################################################

xtx_unique <- xtx_unique[
  SNP_ID %in%
    UNIVERSE
]


c2_unique <- c2_unique[
  SNP_ID %in%
    UNIVERSE
]


XTX_SET <- xtx_unique$SNP_ID

C2_SET <- c2_unique$SNP_ID


cat(
  "\n============================================\n"
)

cat(
  "REVISED POD99.9 SETS INSIDE LOCAL-SCORE REGIONS\n"
)

cat(
  "============================================\n"
)


cat(
  "XtX unique POD99.9 SNPs inside 12 regions: ",
  length(XTX_SET),
  "\n",
  sep = ""
)


cat(
  "C2 unique POD99.9 SNPs inside 17 regions: ",
  length(C2_SET),
  "\n",
  sep = ""
)


############################################################
## 9. REGION-LEVEL POD SUPPORT
############################################################

region_support <- function(
  inside,
  regions,
  analysis
) {

  out <- copy(
    regions
  )


  out[
    ,
    N_POD99.9_SNPs :=
      0L
  ]


  for (
    i in seq_len(
      nrow(out)
    )
  ) {

    ids <- inside[

      Region_ID ==
        out$Region_ID[i],

      unique(
        SNP_ID
      )

    ]


    out$N_POD99.9_SNPs[i] <-
      length(ids)
  }


  out[
    ,
    Contains_POD99.9 :=
      N_POD99.9_SNPs >
      0L
  ]


  out[
    ,
    Analysis :=
      analysis
  ]


  out[
    ,
    .(
      Analysis,
      Region_ID,
      Scaffold_accession,
      Window_start,
      Window_end,
      Gene_ID,
      Gene_symbol,
      N_POD99.9_SNPs,
      Contains_POD99.9
    )
  ]
}


xtx_region_support <- region_support(
  xtx_inside,
  xtx_regions,
  "XtX"
)


c2_region_support <- region_support(
  c2_inside,
  c2_regions,
  "C2"
)


############################################################
## 10. READ AUTHORITATIVE GEA BF>20 ASSOCIATIONS
############################################################

gea <- fread(
  GEA_FILE,
  showProgress = FALSE
)


############################################################
## Inspect / identify required columns
############################################################

GEA_COV <- find_col(
  gea,
  c(
    "Covariate",
    "Canonical_covariate",
    "covariate"
  )
)


GEA_CHR <- find_col(
  gea,
  c(
    "Scaffold_accession",
    "Scaffold",
    "CHROM",
    "chr"
  )
)


GEA_POS <- find_col(
  gea,
  c(
    "Position",
    "POS",
    "Peak_position",
    "SNP_position",
    "position"
  )
)


gea[
  ,
  Canonical_covariate :=
    as.character(
      get(GEA_COV)
    )
]


gea[
  Canonical_covariate ==
    "TemperaturePC1",
  Canonical_covariate :=
    "TempPC1"
]


gea[
  Canonical_covariate ==
    "PrecipitationPC1",
  Canonical_covariate :=
    "PrecipPC1"
]


gea[
  ,
  `:=`(
    Scaffold_norm =
      normalize_scaffold(
        get(GEA_CHR)
      ),

    Position_final =
      as.integer(
        get(GEA_POS)
      )
  )
]


gea[
  ,
  SNP_ID :=
    make_snp_id(
      Scaffold_norm,
      Position_final
    )
]


gea <- gea[
  Canonical_covariate %in%
    COVARIATE_ORDER
]


############################################################
## Restrict GEA exact SNP sets to SAME revised universe
############################################################

gea[
  ,
  In_final_universe :=
    SNP_ID %in%
      UNIVERSE
]


gea_final <- gea[
  In_final_universe == TRUE
]


############################################################
## 11. BUILD SIX GEA SETS
############################################################

GEA_SETS <- setNames(

  lapply(

    COVARIATE_ORDER,

    function(cov) {

      unique(
        gea_final[
          Canonical_covariate ==
            cov,
          SNP_ID
        ]
      )
    }
  ),

  COVARIATE_ORDER
)


GEA_COUNTS <- data.table(

  Covariate =
    COVARIATE_ORDER,

  N_BF20_original =
    vapply(
      COVARIATE_ORDER,
      function(cov) {
        uniqueN(
          gea[
            Canonical_covariate ==
              cov,
            SNP_ID
          ]
        )
      },
      integer(1)
    ),

  N_BF20_in_final_universe =
    vapply(
      COVARIATE_ORDER,
      function(cov) {
        length(
          GEA_SETS[[cov]]
        )
      },
      integer(1)
    )
)


GEA_POOLED <- unique(
  unlist(
    GEA_SETS,
    use.names = FALSE
  )
)


cat(
  "\n============================================\n"
)

cat(
  "GEA SETS IN REVISED SNP UNIVERSE\n"
)

cat(
  "============================================\n"
)


print(
  GEA_COUNTS
)


cat(
  "\nPooled unique GEA SNPs = ",
  length(GEA_POOLED),
  "\n",
  sep = ""
)


############################################################
## 12. XtX vs C2 EXACT SNP OVERLAP
############################################################

xtx_c2_test <- hyper_test(
  XTX_SET,
  C2_SET,
  UNIVERSE,
  "XtX_POD99.9_vs_C2_POD99.9"
)


xtx_c2_ids <- intersect(
  XTX_SET,
  C2_SET
)


xtx_c2_overlap <- merge(

  xtx_unique[
    SNP_ID %in%
      xtx_c2_ids
  ],

  c2_unique[
    SNP_ID %in%
      xtx_c2_ids
  ],

  by = c(
    "SNP_ID",
    "Scaffold_norm",
    "POS"
  ),

  suffixes = c(
    "_XtX",
    "_C2"
  )

)


############################################################
## 13. XtX/C2 × SIX GEA COVARIATES
############################################################

specific_list <- list()

counter <- 1L


for (
  cov in COVARIATE_ORDER
) {

  specific_list[[counter]] <-
    hyper_test(
      XTX_SET,
      GEA_SETS[[cov]],
      UNIVERSE,
      paste0(
        "XtX_POD99.9_vs_",
        cov
      )
    )


  counter <- counter + 1L


  specific_list[[counter]] <-
    hyper_test(
      C2_SET,
      GEA_SETS[[cov]],
      UNIVERSE,
      paste0(
        "C2_POD99.9_vs_",
        cov
      )
    )


  counter <- counter + 1L
}


specific_tests <- rbindlist(
  specific_list
)


specific_tests[
  ,
  BH_FDR :=
    p.adjust(
      P_hypergeometric,
      method = "BH"
    )
]


setorder(
  specific_tests,
  P_hypergeometric
)


############################################################
## 14. POOLED GEA TESTS
############################################################

pooled_tests <- rbindlist(

  list(

    hyper_test(
      XTX_SET,
      GEA_POOLED,
      UNIVERSE,
      "XtX_POD99.9_vs_GEA_pooled"
    ),

    hyper_test(
      C2_SET,
      GEA_POOLED,
      UNIVERSE,
      "C2_POD99.9_vs_GEA_pooled"
    )

  )
)


pooled_tests[
  ,
  BH_FDR :=
    p.adjust(
      P_hypergeometric,
      method = "BH"
    )
]


############################################################
## 15. EXACT XtX/C2 × GEA OVERLAPS
############################################################

cross_list <- list()

counter <- 1L


for (
  cov in COVARIATE_ORDER
) {

  ##########################################################
  ## XtX
  ##########################################################

  ids_xtx <- intersect(
    XTX_SET,
    GEA_SETS[[cov]]
  )


  if (
    length(ids_xtx) > 0L
  ) {

    tmp <- merge(

      xtx_unique[
        SNP_ID %in%
          ids_xtx
      ],

      unique(
        gea_final[
          Canonical_covariate ==
            cov &
          SNP_ID %in%
            ids_xtx,
          .(
            SNP_ID,
            Covariate =
              Canonical_covariate
          )
        ]
      ),

      by =
        "SNP_ID"
    )


    tmp[
      ,
      Selection_scan :=
        "XtX"
    ]


    cross_list[[counter]] <-
      tmp


    counter <- counter + 1L
  }


  ##########################################################
  ## C2
  ##########################################################

  ids_c2 <- intersect(
    C2_SET,
    GEA_SETS[[cov]]
  )


  if (
    length(ids_c2) > 0L
  ) {

    tmp <- merge(

      c2_unique[
        SNP_ID %in%
          ids_c2
      ],

      unique(
        gea_final[
          Canonical_covariate ==
            cov &
          SNP_ID %in%
            ids_c2,
          .(
            SNP_ID,
            Covariate =
              Canonical_covariate
          )
        ]
      ),

      by =
        "SNP_ID"
    )


    tmp[
      ,
      Selection_scan :=
        "C2"
    ]


    cross_list[[counter]] <-
      tmp


    counter <- counter + 1L
  }
}


if (
  length(cross_list) > 0L
) {

  exact_gea_overlap <- rbindlist(
    cross_list,
    fill = TRUE
  )

} else {

  exact_gea_overlap <- data.table()
}


############################################################
## 16. SUMMARY TABLE
############################################################

SUMMARY <- data.table(

  Metric = c(

    "Final_unique_SNP_universe",

    "XtX_final_localScore_regions",

    "C2_final_localScore_regions",

    "XtX_genomewide_POD99.9_outliers",

    "C2_genomewide_POD99.9_outliers",

    "XtX_POD99.9_SNPs_inside_final_regions",

    "C2_POD99.9_SNPs_inside_final_regions",

    "XtX_regions_with_at_least_one_POD99.9_SNP",

    "C2_regions_with_at_least_one_POD99.9_SNP",

    "XtX_C2_exact_shared_POD99.9_SNPs",

    "Pooled_GEA_BF20_unique_SNPs_in_final_universe"

  ),

  Value = c(

    length(UNIVERSE),

    nrow(xtx_regions),

    nrow(c2_regions),

    nrow(xtx_pod),

    nrow(c2_pod),

    length(XTX_SET),

    length(C2_SET),

    sum(
      xtx_region_support$
        Contains_POD99.9
    ),

    sum(
      c2_region_support$
        Contains_POD99.9
    ),

    length(
      xtx_c2_ids
    ),

    length(
      GEA_POOLED
    )

  )

)


############################################################
## 17. WRITE OUTPUTS
############################################################

fwrite(
  SUMMARY,
  file.path(
    OUTDIR,
    "Table_FINAL_gt1Mb_POD_exactSNP_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  xtx_unique,
  file.path(
    OUTDIR,
    "Table_FINAL_XtX_POD99.9_SNPs_inside_12_regions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  c2_unique,
  file.path(
    OUTDIR,
    "Table_FINAL_C2_POD99.9_SNPs_inside_17_regions.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  xtx_region_support,
  file.path(
    OUTDIR,
    "Table_FINAL_XtX_POD99.9_support_by_region.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  c2_region_support,
  file.path(
    OUTDIR,
    "Table_FINAL_C2_POD99.9_support_by_region.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  xtx_c2_test,
  file.path(
    OUTDIR,
    "Table_FINAL_hypergeometric_XtX_vs_C2.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  specific_tests,
  file.path(
    OUTDIR,
    "Table_FINAL_hypergeometric_XtX_C2_vs_GEA_covariates.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  pooled_tests,
  file.path(
    OUTDIR,
    "Table_FINAL_hypergeometric_XtX_C2_vs_GEA_pooled.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  xtx_c2_overlap,
  file.path(
    OUTDIR,
    "Table_FINAL_exact_SNP_overlap_XtX_vs_C2.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  exact_gea_overlap,
  file.path(
    OUTDIR,
    "Table_FINAL_exact_SNP_overlap_XtX_C2_vs_GEA.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  GEA_COUNTS,
  file.path(
    OUTDIR,
    "Table_FINAL_GEA_BF20_counts_in_gt1Mb_universe.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


############################################################
## 18. FINAL REPORT
############################################################

cat(
  "\n============================================\n"
)

cat(
  "FINAL >1 Mb POD / EXACT SNP SUMMARY\n"
)

cat(
  "============================================\n"
)


print(
  SUMMARY
)


cat(
  "\n============================================\n"
)

cat(
  "XtX vs C2 HYPERGEOMETRIC TEST\n"
)

cat(
  "============================================\n"
)


print(
  xtx_c2_test
)


cat(
  "\nExact XtX-C2 shared SNPs:\n"
)


print(
  xtx_c2_overlap
)


cat(
  "\n============================================\n"
)

cat(
  "XtX/C2 vs GEA — 12 COVARIATE TESTS\n"
)

cat(
  "============================================\n"
)


print(
  specific_tests
)


cat(
  "\n============================================\n"
)

cat(
  "XtX/C2 vs POOLED GEA\n"
)

cat(
  "============================================\n"
)


print(
  pooled_tests
)


cat(
  "\n============================================\n"
)

cat(
  "EXACT XtX/C2 × GEA SNP OVERLAPS\n"
)

cat(
  "============================================\n"
)


if (
  nrow(exact_gea_overlap) > 0L
) {

  print(
    exact_gea_overlap
  )

} else {

  cat(
    "No exact overlaps.\n"
  )
}


cat(
  "\n============================================\n"
)

cat(
  "OUTPUT DIRECTORY\n"
)

cat(
  "============================================\n"
)


cat(
  OUTDIR,
  "\n"
)


cat(
  "\nDONE\n"
)
