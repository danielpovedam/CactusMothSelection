 #!/usr/bin/env Rscript

############################################################
## BAYPASS LOCAL-SCORE ANALYSIS
## Cactoblastis cactorum
##
## Analyses:
##   XtX        -> representative seed 12345
##   C2 Host    -> representative seed 12345
##   GEA BF     -> mean BF(dB) across 3 seeds
##
## Datasets:
##   random_snp
##   all_snps
##
## Local-score parameters:
##   min.maf  = 0.05
##   min.nsnp = 20
##   xi       = 2
##   alpha    = 0.01
############################################################


suppressPackageStartupMessages({
  library(data.table)
})


############################################################
## 0. SETTINGS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

UTILS <- paste0(
  "/home/dpoveda/micetral/selection/baypass/",
  "baypass_public/utils/baypass_utils.R"
)

INBASE <- file.path(
  BASE,
  "selection_scan",
  "prepared_statistics"
)

OUTBASE <- file.path(
  BASE,
  "selection_scan",
  "local_score"
)


DATASETS <- c(
  "random_snp",
  "all_snps"
)


REP_SEED <- "12345"


############################################################
## LOCAL-SCORE PARAMETERS
############################################################

MIN_MAF <- 0.05

# BayPass internally uses:
# chr.names > min.nsnp
#
# Therefore MIN_NSNP = 20 means >20 SNPs/scaffold.
MIN_NSNP <- 20


XI_XTX <- 2

XI_C2 <- 2

XI_BF <- 2


ALPHA_XTX <- 0.01

ALPHA_C2 <- 0.01

ALPHA_BF <- 0.01


############################################################
## 1. LOAD BAYPASS UTILITIES
############################################################

if (!file.exists(UTILS)) {

  stop(
    "BayPass utilities file not found: ",
    UTILS
  )
}


source(
  UTILS
)


if (!exists("compute.local.scores")) {

  stop(
    "compute.local.scores() not found after sourcing BayPass utilities."
  )
}


dir.create(
  OUTBASE,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 2. HELPERS
############################################################

write_if_present <- function(
  x,
  file
) {

  if (
    !is.null(x) &&
    nrow(x) > 0
  ) {

    fwrite(
      as.data.table(x),
      file,
      sep = "\t"
    )

  } else {

    fwrite(
      data.table(),
      file,
      sep = "\t"
    )
  }
}


prepare_positions <- function(dt) {

  dt[, CHROM := as.character(CHROM)]

  dt[, POS := as.integer(POS)]

  setorder(
    dt,
    CHROM,
    POS
  )

  return(dt)
}


############################################################
## 3. LOCAL SCORE USING -log10(p)
##
## Used for XtX and C2
############################################################

run_local_score_pvalue <- function(
  dt,
  dataset,
  analysis,
  xi,
  alpha,
  out_scores,
  out_windows
) {

  dt <- copy(dt)

  dt <- prepare_positions(
    dt
  )


  required <- c(
    "CHROM",
    "POS",
    "M_P",
    "log10_p"
  )


  if (!all(required %in% names(dt))) {

    stop(
      analysis,
      " input missing required columns."
    )
  }


  dt[, M_P := as.numeric(M_P)]

  dt[, log10_p := as.numeric(log10_p)]


  dt <- dt[
    !is.na(CHROM) &
    !is.na(POS) &
    !is.na(M_P) &
    !is.na(log10_p)
  ]


  cat(
    "\n----------------------------------------\n"
  )

  cat(
    "Dataset: ",
    dataset,
    "\n",
    sep = ""
  )

  cat(
    "Analysis: ",
    analysis,
    "\n",
    sep = ""
  )

  cat(
    "Input SNPs: ",
    nrow(dt),
    "\n",
    sep = ""
  )


  res <- compute.local.scores(

    snp.position = data.frame(
      chr = dt$CHROM,
      pos = dt$POS
    ),

    snp.pi = dt$M_P,

    snp.pvalue = dt$log10_p,

    xi = xi,

    pval.local.score.thres = alpha,

    min.maf = MIN_MAF,

    min.nsnp = MIN_NSNP,

    plot.pvalhist = FALSE,

    manplot = FALSE
  )


  ##########################################################
  ## Local-score SNP output
  ##########################################################

  scores <- as.data.table(
    res$res.local.scores
  )


  if (nrow(scores) > 0) {

    scores[
      ,
      Dataset := dataset
    ]

    scores[
      ,
      Analysis := analysis
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


  write_if_present(
    scores,
    out_scores
  )


  ##########################################################
  ## Significant windows
  ##########################################################

  windows <- as.data.table(
    res$significant.windows
  )


  if (nrow(windows) > 0) {

    windows[
      ,
      Dataset := dataset
    ]

    windows[
      ,
      Analysis := analysis
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


  write_if_present(
    windows,
    out_windows
  )


  cat(
    "Local-score SNPs retained: ",
    nrow(scores),
    "\n",
    sep = ""
  )

  cat(
    "Significant windows: ",
    nrow(windows),
    "\n",
    sep = ""
  )


  invisible(
    list(
      scores = scores,
      windows = windows
    )
  )
}


############################################################
## 4. LOCAL SCORE USING BF(dB)
##
## Used for GEA
############################################################

run_local_score_bf <- function(
  dt,
  dataset,
  covariate,
  xi,
  alpha,
  out_scores,
  out_windows
) {

  dt <- copy(dt)

  dt <- prepare_positions(
    dt
  )


  required <- c(
    "CHROM",
    "POS",
    "M_P",
    "BF_mean"
  )


  if (!all(required %in% names(dt))) {

    stop(
      "BF input missing required columns for ",
      covariate
    )
  }


  dt[, M_P := as.numeric(M_P)]

  dt[, BF_mean := as.numeric(BF_mean)]


  dt <- dt[
    !is.na(CHROM) &
    !is.na(POS) &
    !is.na(M_P) &
    !is.na(BF_mean)
  ]


  cat(
    "\n----------------------------------------\n"
  )

  cat(
    "Dataset: ",
    dataset,
    "\n",
    sep = ""
  )

  cat(
    "GEA covariate: ",
    covariate,
    "\n",
    sep = ""
  )

  cat(
    "Input SNPs: ",
    nrow(dt),
    "\n",
    sep = ""
  )


  res <- compute.local.scores(

    snp.position = data.frame(
      chr = dt$CHROM,
      pos = dt$POS
    ),

    snp.pi = dt$M_P,

    snp.bf = dt$BF_mean,

    xi = xi,

    pval.local.score.thres = alpha,

    min.maf = MIN_MAF,

    min.nsnp = MIN_NSNP,

    bf.scale = 1,

    plot.pvalhist = FALSE,

    manplot = FALSE
  )


  ##########################################################
  ## Local-score SNP output
  ##########################################################

  scores <- as.data.table(
    res$res.local.scores
  )


  if (nrow(scores) > 0) {

    scores[
      ,
      Dataset := dataset
    ]

    scores[
      ,
      Analysis := "GEA"
    ]

    scores[
      ,
      Covariate := covariate
    ]


    setcolorder(
      scores,
      c(
        "Dataset",
        "Analysis",
        "Covariate",
        setdiff(
          names(scores),
          c(
            "Dataset",
            "Analysis",
            "Covariate"
          )
        )
      )
    )
  }


  write_if_present(
    scores,
    out_scores
  )


  ##########################################################
  ## Significant windows
  ##########################################################

  windows <- as.data.table(
    res$significant.windows
  )


  if (nrow(windows) > 0) {

    windows[
      ,
      Dataset := dataset
    ]

    windows[
      ,
      Analysis := "GEA"
    ]

    windows[
      ,
      Covariate := covariate
    ]


    setcolorder(
      windows,
      c(
        "Dataset",
        "Analysis",
        "Covariate",
        setdiff(
          names(windows),
          c(
            "Dataset",
            "Analysis",
            "Covariate"
          )
        )
      )
    )
  }


  write_if_present(
    windows,
    out_windows
  )


  cat(
    "Local-score SNPs retained: ",
    nrow(scores),
    "\n",
    sep = ""
  )

  cat(
    "Significant windows: ",
    nrow(windows),
    "\n",
    sep = ""
  )


  invisible(
    list(
      scores = scores,
      windows = windows
    )
  )
}


############################################################
## 5. GLOBAL SUMMARY CONTAINER
############################################################

SUMMARY_LIST <- list()

summary_counter <- 1


############################################################
## 6. MAIN LOOP
############################################################

for (dataset in DATASETS) {

  cat(
    "\n========================================\n"
  )

  cat(
    "PROCESSING DATASET: ",
    dataset,
    "\n",
    sep = ""
  )

  cat(
    "========================================\n"
  )


  in_ds <- file.path(
    INBASE,
    dataset
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
  ## 6.1 XtX
  ##########################################################

  xtx_file <- file.path(
    in_ds,
    paste0(
      "XtX_seed",
      REP_SEED,
      ".tsv"
    )
  )


  if (!file.exists(xtx_file)) {

    stop(
      "Missing XtX input: ",
      xtx_file
    )
  }


  xtx <- fread(
    xtx_file,
    showProgress = FALSE
  )


  xtx_dir <- file.path(
    out_ds,
    "XtX"
  )


  dir.create(
    xtx_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  xtx_result <- run_local_score_pvalue(

    dt = xtx,

    dataset = dataset,

    analysis = "XtX",

    xi = XI_XTX,

    alpha = ALPHA_XTX,

    out_scores = file.path(
      xtx_dir,
      "XtX_local_scores.tsv"
    ),

    out_windows = file.path(
      xtx_dir,
      "XtX_candidate_windows.tsv"
    )
  )


  SUMMARY_LIST[[summary_counter]] <- data.table(

    Dataset = dataset,

    Analysis = "XtX",

    Covariate = NA_character_,

    N_input_SNPs = nrow(xtx),

    N_local_score_SNPs =
      nrow(xtx_result$scores),

    N_candidate_windows =
      nrow(xtx_result$windows)
  )


  summary_counter <- summary_counter + 1


  rm(
    xtx,
    xtx_result
  )


  gc(
    FALSE
  )


  ##########################################################
  ## 6.2 C2 HOST
  ##########################################################

  c2_file <- file.path(
    in_ds,
    paste0(
      "C2_Host_seed",
      REP_SEED,
      ".tsv"
    )
  )


  if (!file.exists(c2_file)) {

    stop(
      "Missing C2 input: ",
      c2_file
    )
  }


  c2 <- fread(
    c2_file,
    showProgress = FALSE
  )


  c2_dir <- file.path(
    out_ds,
    "C2_Host"
  )


  dir.create(
    c2_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  c2_result <- run_local_score_pvalue(

    dt = c2,

    dataset = dataset,

    analysis = "C2_Host",

    xi = XI_C2,

    alpha = ALPHA_C2,

    out_scores = file.path(
      c2_dir,
      "C2_Host_local_scores.tsv"
    ),

    out_windows = file.path(
      c2_dir,
      "C2_Host_candidate_windows.tsv"
    )
  )


  SUMMARY_LIST[[summary_counter]] <- data.table(

    Dataset = dataset,

    Analysis = "C2_Host",

    Covariate = "Host_C2",

    N_input_SNPs = nrow(c2),

    N_local_score_SNPs =
      nrow(c2_result$scores),

    N_candidate_windows =
      nrow(c2_result$windows)
  )


  summary_counter <- summary_counter + 1


  rm(
    c2,
    c2_result
  )


  gc(
    FALSE
  )


  ##########################################################
  ## 6.3 GEA
  ##########################################################

  bf_file <- file.path(
    in_ds,
    "GEA_BF_mean3seeds.tsv"
  )


  if (!file.exists(bf_file)) {

    stop(
      "Missing GEA BF input: ",
      bf_file
    )
  }


  bf <- fread(
    bf_file,
    showProgress = FALSE
  )


  gea_dir <- file.path(
    out_ds,
    "GEA"
  )


  dir.create(
    gea_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  covariates <- unique(
    bf$Covariate
  )


  ##########################################################
  ## Preserve intended covariate order
  ##########################################################

  covariate_order <- c(
    "TempPC1",
    "PrecipPC1",
    "AridityIndex",
    "Elevation",
    "Latitude",
    "Longitude",
    "BIO1",
    "BIO12"
  )


  covariates <- covariate_order[
    covariate_order %in% covariates
  ]


  for (covariate in covariates) {

    bf_cov <- bf[
      Covariate == covariate
    ]


    cov_dir <- file.path(
      gea_dir,
      covariate
    )


    dir.create(
      cov_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )


    bf_result <- run_local_score_bf(

      dt = bf_cov,

      dataset = dataset,

      covariate = covariate,

      xi = XI_BF,

      alpha = ALPHA_BF,

      out_scores = file.path(
        cov_dir,
        paste0(
          covariate,
          "_local_scores.tsv"
        )
      ),

      out_windows = file.path(
        cov_dir,
        paste0(
          covariate,
          "_candidate_windows.tsv"
        )
      )
    )


    SUMMARY_LIST[[summary_counter]] <- data.table(

      Dataset = dataset,

      Analysis = "GEA",

      Covariate = covariate,

      N_input_SNPs =
        nrow(bf_cov),

      N_local_score_SNPs =
        nrow(bf_result$scores),

      N_candidate_windows =
        nrow(bf_result$windows)
    )


    summary_counter <- summary_counter + 1


    rm(
      bf_cov,
      bf_result
    )


    gc(
      FALSE
    )
  }


  rm(
    bf
  )


  gc(
    FALSE
  )
}


############################################################
## 7. WRITE RUN SUMMARY
############################################################

RUN_SUMMARY <- rbindlist(
  SUMMARY_LIST,
  use.names = TRUE,
  fill = TRUE
)


fwrite(
  RUN_SUMMARY,
  file.path(
    OUTBASE,
    "local_score_run_summary.tsv"
  ),
  sep = "\t"
)


############################################################
## 8. SETTINGS TABLE
############################################################

SETTINGS <- data.table(

  Parameter = c(
    "Representative_seed_XtX",
    "Representative_seed_C2",
    "BF_seed_combination",
    "min_maf",
    "min_nsnp",
    "xi_XtX",
    "xi_C2",
    "xi_BF",
    "alpha_XtX",
    "alpha_C2",
    "alpha_BF"
  ),

  Value = c(
    REP_SEED,
    REP_SEED,
    "mean_12345_67890_13579",
    MIN_MAF,
    MIN_NSNP,
    XI_XTX,
    XI_C2,
    XI_BF,
    ALPHA_XTX,
    ALPHA_C2,
    ALPHA_BF
  )
)


fwrite(
  SETTINGS,
  file.path(
    OUTBASE,
    "local_score_settings.tsv"
  ),
  sep = "\t"
)


############################################################
## 9. PRINT SUMMARY
############################################################

cat(
  "\n========================================\n"
)

cat(
  "LOCAL-SCORE RUN SUMMARY\n"
)

cat(
  "========================================\n"
)


print(
  RUN_SUMMARY
)


cat(
  "\n========================================\n"
)

cat(
  "SETTINGS\n"
)

cat(
  "========================================\n"
)


print(
  SETTINGS
)


cat(
  "\nDONE\n"
)
