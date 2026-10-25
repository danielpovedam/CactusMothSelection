#!/usr/bin/env Rscript

############################################################
## 49A_simulate_neutral_POD.R
##
## Cactoblastis cactorum
##
## Generate neutral pseudo-observed data (POD) using:
##
##   - empirical BayPass Omega
##   - empirical beta distribution of ancestral frequency
##   - empirical SNP-specific population sample sizes
##   - 100,000 simulated loci
##
## This POD contains NO simulated environmental or host
## effects. It therefore represents the neutral expectation
## conditional on the fitted population covariance structure.
############################################################


options(stringsAsFactors = FALSE)


############################################################
## 0. PATHS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

UTILS <- paste0(
  "/home/dpoveda/micetral/selection/baypass/",
  "baypass_public/utils/baypass_utils.R"
)

GENO_FILE <- paste0(
  "/home/dpoveda/micetral/Population_snps/",
  "baypass_preparation/baypass/",
  "all_snps.baypass.geno"
)

OMEGA_FILE <- file.path(
  BASE,
  "runs_10x",
  "all_snps",
  "seed_12345",
  "BP_all_snps_seed12345_mat_omega.out"
)

BETA_FILE <- file.path(
  BASE,
  "runs_10x",
  "all_snps",
  "seed_12345",
  "BP_all_snps_seed12345_summary_beta_params.out"
)

POP_FILE <- paste0(
  "/home/dpoveda/micetral/Population_snps/",
  "baypass_preparation/input/",
  "baypass_population_order.txt"
)

OUTDIR <- file.path(
  BASE,
  "POD_analysis",
  "01_simulation"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 1. SETTINGS — FROZEN BEFORE OVERLAP ANALYSIS
############################################################

N_POD <- 100000L

SIMULATION_SEED <- 20260827L

PI_MAF <- 0.05

SUFFIX <- "POD_neutral_100k"


############################################################
## 2. INPUT QC
############################################################

FILES <- c(
  UTILS,
  GENO_FILE,
  OMEGA_FILE,
  BETA_FILE,
  POP_FILE
)

missing_files <- FILES[
  !file.exists(FILES)
]

if (
  length(missing_files) > 0
) {

  stop(
    "\nMissing input file(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


cat(
  "\n============================================\n"
)

cat(
  "49A — BAYPASS NEUTRAL POD SIMULATION\n"
)

cat(
  "============================================\n\n"
)


############################################################
## 3. LOAD BAYPASS UTILITIES
############################################################

source(
  UTILS
)

if (
  !exists(
    "simulate.baypass"
  )
) {

  stop(
    "simulate.baypass() was not loaded."
  )
}


############################################################
## 4. POPULATION ORDER
############################################################

POP <- scan(
  POP_FILE,
  what = character(),
  quiet = TRUE
)


if (
  length(POP) != 28L
) {

  stop(
    "Expected 28 populations; found ",
    length(POP)
  )
}


cat(
  "[49A] Populations: ",
  length(POP),
  "\n",
  sep = ""
)


############################################################
## 5. READ EMPIRICAL OMEGA
############################################################

OMEGA <- as.matrix(
  read.table(
    OMEGA_FILE,
    header = FALSE
  )
)


storage.mode(
  OMEGA
) <- "double"


if (
  !all(
    dim(OMEGA) == c(28L, 28L)
  )
) {

  stop(
    "Omega is not 28 x 28."
  )
}


if (
  max(
    abs(
      OMEGA -
        t(OMEGA)
    )
  ) > 1e-5
) {

  stop(
    "Omega matrix is not symmetric."
  )
}


eigenvalues <- eigen(
  OMEGA,
  symmetric = TRUE,
  only.values = TRUE
)$values


if (
  min(eigenvalues) <= 0
) {

  stop(
    "Omega is not positive definite. Minimum eigenvalue = ",
    min(eigenvalues)
  )
}


dimnames(
  OMEGA
) <- list(
  POP,
  POP
)


cat(
  "[49A] Omega dimensions: ",
  nrow(OMEGA),
  " x ",
  ncol(OMEGA),
  "\n",
  sep = ""
)


cat(
  "[49A] Minimum Omega eigenvalue: ",
  signif(
    min(eigenvalues),
    6
  ),
  "\n",
  sep = ""
)


############################################################
## 6. READ EMPIRICAL BETA PARAMETERS
############################################################

BETA <- read.table(
  BETA_FILE,
  header = TRUE,
  stringsAsFactors = FALSE
)


if (
  !all(
    c(
      "PARAM",
      "Mean"
    ) %in%
      names(BETA)
  )
) {

  stop(
    "Unexpected beta parameter file structure."
  )
}


A_BETA <- BETA$Mean[
  BETA$PARAM ==
    "a_beta_pi"
]


B_BETA <- BETA$Mean[
  BETA$PARAM ==
    "b_beta_pi"
]


if (
  length(A_BETA) != 1L ||
  length(B_BETA) != 1L
) {

  stop(
    "Could not extract a_beta_pi / b_beta_pi."
  )
}


if (
  A_BETA <= 0 ||
  B_BETA <= 0
) {

  stop(
    "Invalid beta distribution parameters."
  )
}


BETA_PI <- c(
  A_BETA,
  B_BETA
)


cat(
  "[49A] a_beta_pi = ",
  A_BETA,
  "\n",
  sep = ""
)


cat(
  "[49A] b_beta_pi = ",
  B_BETA,
  "\n",
  sep = ""
)


############################################################
## 7. READ EMPIRICAL GENOTYPE COUNT MATRIX
############################################################

cat(
  "\n[49A] Reading empirical BayPass count matrix...\n"
)


G <- as.matrix(
  read.table(
    GENO_FILE,
    header = FALSE
  )
)


storage.mode(
  G
) <- "integer"


if (
  nrow(G) != 132283L
) {

  stop(
    "Expected 132,283 empirical SNPs; found ",
    nrow(G)
  )
}


if (
  ncol(G) != 56L
) {

  stop(
    "Expected 56 BayPass count columns; found ",
    ncol(G)
  )
}


if (
  any(
    G < 0L
  )
) {

  stop(
    "Negative allele count detected."
  )
}


############################################################
## 8. BUILD EXACT EMPIRICAL SAMPLE-SIZE MATRIX
##
## For each population:
##
## N_ij = allele1_count + allele2_count
##
## This preserves the empirical distribution of genotype
## missingness among SNPs.
############################################################

NN <- matrix(
  NA_integer_,
  nrow = nrow(G),
  ncol = 28L
)


for (
  j in seq_len(28L)
) {

  c1 <- 2L * j - 1L
  c2 <- 2L * j


  NN[
    ,
    j
  ] <- G[, c1] + G[, c2]
}


colnames(
  NN
) <- POP


if (
  any(
    NN < 0L
  )
) {

  stop(
    "Negative sample size detected."
  )
}


############################################################
## Sample-size QC
############################################################

NN_SUMMARY <- data.frame(

  Population =
    POP,

  Min_haploid_count =
    apply(
      NN,
      2,
      min
    ),

  Median_haploid_count =
    apply(
      NN,
      2,
      median
    ),

  Mean_haploid_count =
    apply(
      NN,
      2,
      mean
    ),

  Max_haploid_count =
    apply(
      NN,
      2,
      max
    )
)


write.table(
  NN_SUMMARY,
  file = file.path(
    OUTDIR,
    "QC_empirical_sample_size_distribution.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "\n[49A] Empirical haploid-count summary:\n"
)


print(
  NN_SUMMARY
)


############################################################
## 9. GENERATE POD
############################################################

set.seed(
  SIMULATION_SEED
)


oldwd <- getwd()


setwd(
  OUTDIR
)


cat(
  "\n[49A] Simulating ",
  format(
    N_POD,
    big.mark = ","
  ),
  " neutral SNPs...\n",
  sep = ""
)


SIM <- simulate.baypass(

  omega.mat =
    OMEGA,

  nsnp =
    N_POD,

  beta.coef =
    NA,

  beta.pi =
    BETA_PI,

  sample.size =
    NN,

  pi.maf =
    PI_MAF,

  suffix =
    SUFFIX,

  print.sim.params.values =
    TRUE,

  remove.fixed.loci =
    FALSE,

  output.bayenv.format =
    FALSE,

  coverage =
    NA
)


setwd(
  oldwd
)


############################################################
## 10. OUTPUT QC
############################################################

POD_FILE <- file.path(
  OUTDIR,
  paste0(
    "G.",
    SUFFIX
  )
)


if (
  !file.exists(
    POD_FILE
  )
) {

  stop(
    "Expected POD genotype file was not created:\n",
    POD_FILE
  )
}


POD_G <- as.matrix(
  read.table(
    POD_FILE,
    header = FALSE
  )
)


if (
  nrow(POD_G) !=
    N_POD
) {

  stop(
    "Expected ",
    N_POD,
    " POD rows; found ",
    nrow(POD_G)
  )
}


if (
  ncol(POD_G) != 56L
) {

  stop(
    "Expected 56 POD columns; found ",
    ncol(POD_G)
  )
}


############################################################
## Global observed MAF of simulated data
############################################################

ALT_TOTAL <- rowSums(
  POD_G[
    ,
    seq(
      1,
      56,
      by = 2
    ),
    drop = FALSE
  ]
)


REF_TOTAL <- rowSums(
  POD_G[
    ,
    seq(
      2,
      56,
      by = 2
    ),
    drop = FALSE
  ]
)


TOTAL <- ALT_TOTAL +
  REF_TOTAL


FREQ <- ALT_TOTAL /
  TOTAL


OBS_MAF <- pmin(
  FREQ,
  1 - FREQ
)


POD_QC <- data.frame(

  Metric = c(
    "N_POD_SNPs",
    "N_populations",
    "N_count_columns",
    "a_beta_pi",
    "b_beta_pi",
    "pi_maf_simulation",
    "simulation_seed",
    "Observed_MAF_min",
    "Observed_MAF_median",
    "Observed_MAF_mean",
    "Observed_MAF_max",
    "N_observed_fixed",
    "N_observed_MAF_lt_0.05"
  ),

  Value = c(
    N_POD,
    28,
    56,
    A_BETA,
    B_BETA,
    PI_MAF,
    SIMULATION_SEED,
    min(OBS_MAF),
    median(OBS_MAF),
    mean(OBS_MAF),
    max(OBS_MAF),
    sum(OBS_MAF == 0),
    sum(OBS_MAF < 0.05)
  )
)


write.table(
  POD_QC,
  file = file.path(
    OUTDIR,
    "QC_POD_simulation_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
## Save simulation object for reproducibility
############################################################

saveRDS(
  SIM,
  file = file.path(
    OUTDIR,
    "POD_simulation_object.rds"
  ),
  compress = "xz"
)


############################################################
## Parameters
############################################################

writeLines(
  c(

    paste0(
      "Simulation date: ",
      Sys.Date()
    ),

    paste0(
      "Simulation seed: ",
      SIMULATION_SEED
    ),

    paste0(
      "POD SNPs: ",
      N_POD
    ),

    paste0(
      "Omega file: ",
      OMEGA_FILE
    ),

    paste0(
      "Beta parameter file: ",
      BETA_FILE
    ),

    paste0(
      "a_beta_pi: ",
      A_BETA
    ),

    paste0(
      "b_beta_pi: ",
      B_BETA
    ),

    paste0(
      "pi.maf: ",
      PI_MAF
    ),

    "sample.size: full empirical 132283 x 28 haploid-count matrix",

    "beta.coef: NA; neutral simulation",

    "remove.fixed.loci: FALSE",

    paste0(
      "Output POD file: ",
      POD_FILE
    )
  ),

  con = file.path(
    OUTDIR,
    "POD_simulation_parameters.txt"
  )
)


capture.output(
  sessionInfo(),
  file = file.path(
    OUTDIR,
    "sessionInfo_49A.txt"
  )
)


cat(
  "\n============================================\n"
)

cat(
  "49A COMPLETE\n"
)

cat(
  "============================================\n\n"
)


cat(
  "POD file:\n",
  POD_FILE,
  "\n\n",
  sep = ""
)


cat(
  "Number of POD SNPs: ",
  nrow(POD_G),
  "\n",
  sep = ""
)


cat(
  "Observed simulated SNPs with MAF < 0.05: ",
  sum(
    OBS_MAF < 0.05
  ),
  "\n",
  sep = ""
)


cat(
  "Observed fixed simulated SNPs: ",
  sum(
    OBS_MAF == 0
  ),
  "\n",
  sep = ""
)


cat(
  "\n[49A] DONE\n"
)
