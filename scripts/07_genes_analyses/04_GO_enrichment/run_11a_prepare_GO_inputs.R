#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

############################################################
## 40A_FINAL_gt1Mb_prepare_GO_inputs.R
##
## FINAL GO INPUT PREPARATION
## Cactoblastis cactorum
##
## Candidate sets:
##   XtX       = peak-associated genes from 12 final >1 Mb regions
##   C2        = peak-associated genes from 17 final >1 Mb regions
##   GEA       = unchanged BF>20 / 1-kb candidate loci
##
## Background:
##   genes overlapping or <=5 kb from any SNP screened in
##   the complete all-SNP BayPass dataset.
##
## POD99.9 is NOT used to redefine GO candidate genes.
############################################################


# ============================================================
# PATHS
# ============================================================

BASE <- "/home/dpoveda/micetral/Population_snps"
BAYPASS <- file.path(BASE, "baypass")
SEL <- file.path(BAYPASS, "selection_scan")

META_FILE <- file.path(
  BASE,
  "baypass_preparation",
  "metadata",
  "all_snps.baypass.snp_metadata.tsv"
)

GFF_FILE <- file.path(
  BASE,
  "Ccac_functional-annotation",
  "02_agat",
  "Ccac",
  "out_abril4",
  "Ccac.annotations.corrected.gff3"
)

GO_FILE <- file.path(
  BASE,
  "Ccac_functional-annotation",
  "02_agat",
  "Ccac",
  "out_abril4",
  "GO.txt"
)

# FINAL >1 Mb candidate tables
XTX_FILE <- file.path(
  SEL,
  "FINAL_gt1Mb",
  "Table_XtX_candidate_regions_allSNP_gt1Mb_FINAL.tsv"
)

C2_FILE <- file.path(
  SEL,
  "FINAL_gt1Mb",
  "Table_C2_Host_candidate_regions_allSNP_gt1Mb_FINAL.tsv"
)

# GEA unchanged
GEA_FILE <- file.path(
  SEL,
  "FINAL_GEA_BF20_1kb",
  "Table_GEA_BF20_candidate_loci_1kb_PRIMARY_FINAL.tsv"
)

# New directory: preserve old GO results
OUTDIR <- file.path(
  SEL,
  "GO_ENRICHMENT_FINAL_gt1Mb"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# SETTINGS
# ============================================================

MAX_DISTANCE <- 5000L

EXPECTED_SNPS <- 132283L
EXPECTED_XTX_REGIONS <- 12L
EXPECTED_C2_REGIONS <- 17L
EXPECTED_GEA_LOCI <- 167L

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

check_file <- function(path) {

  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0
  ) {
    stop("Missing or empty file:\n", path)
  }
}


clean_na <- function(x) {

  x <- trimws(as.character(x))

  x[
    x %in% c(
      "",
      "NA",
      "N/A",
      "<NA>"
    )
  ] <- NA_character_

  x
}


clean_gene_id <- function(x) {

  x <- clean_na(x)

  sub(
    "\\.t[0-9]+$",
    "",
    x
  )
}


valid_gene_id <- function(x) {

  !is.na(x) &
    grepl(
      "^Ccac_g[0-9]+$",
      x
    )
}


extract_gene_id_from_attributes <- function(x) {

  x <- as.character(x)

  m <- regexpr(
    "ID=Ccac_g[0-9]+",
    x,
    perl = TRUE
  )

  out <- rep(
    NA_character_,
    length(x)
  )

  hit <- m > 0

  out[hit] <- regmatches(
    x,
    m
  )[hit]

  out <- sub(
    "^ID=",
    "",
    out
  )

  clean_gene_id(out)
}


write_gene_list <- function(x, path) {

  x <- sort(
    unique(
      x[
        valid_gene_id(x)
      ]
    )
  )

  fwrite(
    data.table(
      Gene_ID = x
    ),
    path,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )
}


# ============================================================
# INPUT CHECKS
# ============================================================

for (
  f in c(
    META_FILE,
    GFF_FILE,
    GO_FILE,
    XTX_FILE,
    C2_FILE,
    GEA_FILE
  )
) {
  check_file(f)
}


cat(
  "\n============================================================\n",
  "40A FINAL gt1Mb GO INPUT PREPARATION\n",
  "============================================================\n",
  sep = ""
)


# ============================================================
# 1. ALL SCREENED BAYPASS SNPs
# ============================================================

snps <- fread(
  META_FILE,
  showProgress = FALSE
)

required_meta <- c(
  "BayPass_row",
  "CHROM",
  "SNP",
  "POS",
  "ALT",
  "REF"
)

miss <- setdiff(
  required_meta,
  names(snps)
)

if (length(miss) > 0L) {

  stop(
    "Metadata missing columns: ",
    paste(
      miss,
      collapse = ", "
    )
  )
}


snps[
  ,
  `:=`(
    BayPass_row =
      as.integer(BayPass_row),

    CHROM =
      as.character(CHROM),

    SNP =
      as.character(SNP),

    POS =
      as.integer(POS)
  )
]


if (nrow(snps) != EXPECTED_SNPS) {

  stop(
    "Expected ",
    EXPECTED_SNPS,
    " SNP records; found ",
    nrow(snps)
  )
}


if (
  uniqueN(
    snps$BayPass_row
  ) != EXPECTED_SNPS
) {

  stop(
    "BayPass_row is not unique."
  )
}


# ============================================================
# 2. GENE COORDINATES FROM GFF3
# ============================================================

awk_cmd <- sprintf(
  "awk -F '\\t' 'BEGIN{OFS=\"\\t\"} $0 !~ /^#/ && $3==\"gene\" {print $1,$4,$5,$9}' %s",
  shQuote(GFF_FILE)
)


genes <- fread(
  cmd = awk_cmd,
  header = FALSE,
  sep = "\t",
  quote = "",
  showProgress = FALSE
)


if (
  ncol(genes) != 4L
) {

  stop(
    "Could not extract four columns from GFF3 gene features."
  )
}


setnames(
  genes,
  c(
    "CHROM",
    "Gene_start",
    "Gene_end",
    "Attributes"
  )
)


genes[
  ,
  `:=`(
    CHROM =
      as.character(CHROM),

    Gene_start =
      as.integer(Gene_start),

    Gene_end =
      as.integer(Gene_end),

    Gene_ID =
      extract_gene_id_from_attributes(
        Attributes
      )
  )
]


genes <- unique(
  genes[
    valid_gene_id(Gene_ID),
    .(
      CHROM,
      Gene_start,
      Gene_end,
      Gene_ID
    )
  ]
)


if (
  nrow(genes) == 0L
) {

  stop(
    "No valid Ccac_g gene features extracted."
  )
}


if (
  any(
    genes$Gene_start >
      genes$Gene_end
  )
) {

  stop(
    "Invalid gene coordinates detected."
  )
}


# ============================================================
# 3. COMMON BAYPASS-ACCESSIBLE GENE BACKGROUND
# ============================================================

snp_intervals <- snps[
  ,
  .(
    CHROM,
    start =
      pmax(
        1L,
        POS - MAX_DISTANCE
      ),
    end =
      POS + MAX_DISTANCE,
    BayPass_row,
    SNP,
    POS
  )
]


gene_intervals <- genes[
  ,
  .(
    CHROM,
    start =
      Gene_start,
    end =
      Gene_end,
    Gene_ID
  )
]


setkey(
  gene_intervals,
  CHROM,
  start,
  end
)


accessible_hits <- foverlaps(
  x = snp_intervals,
  y = gene_intervals,
  by.x = c(
    "CHROM",
    "start",
    "end"
  ),
  by.y = c(
    "CHROM",
    "start",
    "end"
  ),
  type = "any",
  nomatch = 0L
)


if (
  nrow(accessible_hits) == 0L
) {

  stop(
    "No SNP-to-gene hits within 5 kb."
  )
}


background_genes <- sort(
  unique(
    accessible_hits$Gene_ID
  )
)


background_genes <- background_genes[
  valid_gene_id(
    background_genes
  )
]


snp_gene_map <- unique(
  accessible_hits[
    ,
    .(
      BayPass_row,
      SNP,
      POS,
      Gene_ID
    )
  ]
)


setorder(
  snp_gene_map,
  BayPass_row,
  Gene_ID
)


fwrite(
  snp_gene_map,
  file.path(
    OUTDIR,
    "BayPass_allSNP_SNP_to_gene_within5kb.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write_gene_list(
  background_genes,
  file.path(
    OUTDIR,
    "GO_background_gene_list.txt"
  )
)


background_accessibility <- snp_gene_map[
  ,
  .(
    N_accessible_BayPass_SNPs =
      uniqueN(BayPass_row)
  ),
  by =
    Gene_ID
]


setorder(
  background_accessibility,
  -N_accessible_BayPass_SNPs,
  Gene_ID
)


fwrite(
  background_accessibility,
  file.path(
    OUTDIR,
    "GO_background_gene_accessibility.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 4. GENE-LEVEL GO DICTIONARY
# ============================================================

go_raw <- fread(
  GO_FILE,
  header = FALSE,
  sep = "\t",
  fill = TRUE,
  quote = "",
  showProgress = FALSE
)


if (
  ncol(go_raw) < 2L
) {

  stop(
    "GO.txt has fewer than two columns."
  )
}


go_raw <- go_raw[
  ,
  1:2
]


setnames(
  go_raw,
  c(
    "Transcript_ID",
    "GO_string"
  )
)


go_raw[
  ,
  `:=`(
    Transcript_ID =
      clean_na(Transcript_ID),

    Gene_ID =
      clean_gene_id(
        Transcript_ID
      ),

    GO_string =
      clean_na(
        GO_string
      )
  )
]


go_raw <- go_raw[
  valid_gene_id(Gene_ID) &
    !is.na(GO_string)
]


go_long <- go_raw[
  ,
  .(
    GO_ID =
      trimws(
        unlist(
          strsplit(
            GO_string,
            ",",
            fixed = TRUE
          )
        )
      )
  ),
  by =
    Gene_ID
]


go_long[
  ,
  GO_ID :=
    gsub(
      "\\s+",
      "",
      GO_ID
    )
]


go_long <- unique(
  go_long[
    grepl(
      "^GO:[0-9]{7}$",
      GO_ID
    )
  ],
  by = c(
    "Gene_ID",
    "GO_ID"
  )
)


setorder(
  go_long,
  Gene_ID,
  GO_ID
)


if (
  nrow(go_long) == 0L
) {

  stop(
    "No valid GO assignments recovered."
  )
}


fwrite(
  go_long,
  file.path(
    OUTDIR,
    "Ccactorum_gene2GO_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


genes_with_go <- unique(
  go_long$Gene_ID
)


# ============================================================
# 5. FINAL XtX PEAK-ASSOCIATED GENES
# ============================================================

xtx <- fread(
  XTX_FILE,
  showProgress = FALSE
)


if (
  nrow(xtx) !=
    EXPECTED_XTX_REGIONS
) {

  stop(
    "Expected ",
    EXPECTED_XTX_REGIONS,
    " XtX regions; found ",
    nrow(xtx)
  )
}


required_xtx <- c(
  "Region_ID",
  "Gene_ID",
  "Distance_to_gene_bp"
)


miss <- setdiff(
  required_xtx,
  names(xtx)
)


if (
  length(miss) > 0L
) {

  stop(
    "XtX table missing: ",
    paste(
      miss,
      collapse = ", "
    )
  )
}


xtx[
  ,
  `:=`(
    Gene_ID =
      clean_gene_id(
        Gene_ID
      ),

    Distance_to_gene_bp =
      suppressWarnings(
        as.numeric(
          Distance_to_gene_bp
        )
      )
  )
]


xtx_genes <- sort(
  unique(
    xtx[
      valid_gene_id(Gene_ID) &
        !is.na(
          Distance_to_gene_bp
        ) &
        abs(
          Distance_to_gene_bp
        ) <= MAX_DISTANCE,
      Gene_ID
    ]
  )
)


# ============================================================
# 6. FINAL C2 PEAK-ASSOCIATED GENES
# ============================================================

c2 <- fread(
  C2_FILE,
  showProgress = FALSE
)


if (
  nrow(c2) !=
    EXPECTED_C2_REGIONS
) {

  stop(
    "Expected ",
    EXPECTED_C2_REGIONS,
    " C2 regions; found ",
    nrow(c2)
  )
}


required_c2 <- c(
  "Region_ID",
  "Gene_ID",
  "Distance_to_gene_bp"
)


miss <- setdiff(
  required_c2,
  names(c2)
)


if (
  length(miss) > 0L
) {

  stop(
    "C2 table missing: ",
    paste(
      miss,
      collapse = ", "
    )
  )
}


c2[
  ,
  `:=`(
    Gene_ID =
      clean_gene_id(
        Gene_ID
      ),

    Distance_to_gene_bp =
      suppressWarnings(
        as.numeric(
          Distance_to_gene_bp
        )
      )
  )
]


c2_genes <- sort(
  unique(
    c2[
      valid_gene_id(Gene_ID) &
        !is.na(
          Distance_to_gene_bp
        ) &
        abs(
          Distance_to_gene_bp
        ) <= MAX_DISTANCE,
      Gene_ID
    ]
  )
)


# ============================================================
# 7. GEA PEAK-ASSOCIATED GENES
# ============================================================

gea <- fread(
  GEA_FILE,
  showProgress = FALSE
)


if (
  nrow(gea) !=
    EXPECTED_GEA_LOCI
) {

  stop(
    "Expected ",
    EXPECTED_GEA_LOCI,
    " GEA loci; found ",
    nrow(gea)
  )
}


required_gea <- c(
  "Covariate",
  "GEA_locus_ID",
  "Gene_ID",
  "Distance_to_gene_bp"
)


miss <- setdiff(
  required_gea,
  names(gea)
)


if (
  length(miss) > 0L
) {

  stop(
    "GEA table missing: ",
    paste(
      miss,
      collapse = ", "
    )
  )
}


gea[
  ,
  `:=`(
    Covariate =
      as.character(
        Covariate
      ),

    Gene_ID =
      clean_gene_id(
        Gene_ID
      ),

    Distance_to_gene_bp =
      suppressWarnings(
        as.numeric(
          Distance_to_gene_bp
        )
      )
  )
]


gea_rows <- gea[
  Covariate %in%
    PRIMARY_COVARIATES &
    valid_gene_id(
      Gene_ID
    ) &
    !is.na(
      Distance_to_gene_bp
    ) &
    abs(
      Distance_to_gene_bp
    ) <= MAX_DISTANCE
]


gea_pooled_genes <- sort(
  unique(
    gea_rows$Gene_ID
  )
)


# ============================================================
# 8. LONG-FORM CANDIDATE TABLE
# ============================================================

candidate_long <- rbindlist(
  c(
    list(
      data.table(
        Analysis = "XtX",
        Gene_ID = xtx_genes
      ),
      data.table(
        Analysis = "C2_Host",
        Gene_ID = c2_genes
      ),
      data.table(
        Analysis = "GEA_pooled",
        Gene_ID = gea_pooled_genes
      )
    ),
    lapply(
      PRIMARY_COVARIATES,
      function(cv) {

        data.table(
          Analysis =
            paste0(
              "GEA_",
              cv
            ),

          Gene_ID =
            sort(
              unique(
                gea_rows[
                  Covariate ==
                    cv,
                  Gene_ID
                ]
              )
            )
        )
      }
    )
  ),
  use.names = TRUE,
  fill = TRUE
)


candidate_long <- unique(
  candidate_long,
  by = c(
    "Analysis",
    "Gene_ID"
  )
)


candidate_long[
  ,
  Has_any_GO :=
    Gene_ID %in%
      genes_with_go
]


candidate_long[
  ,
  In_background :=
    Gene_ID %in%
      background_genes
]


missing_from_background <- candidate_long[
  In_background == FALSE
]


if (
  nrow(
    missing_from_background
  ) > 0L
) {

  fwrite(
    missing_from_background,
    file.path(
      OUTDIR,
      "ERROR_candidate_genes_missing_from_background.tsv"
    ),
    sep = "\t",
    quote = FALSE
  )

  stop(
    "At least one candidate gene is absent from the accessible background."
  )
}


setorder(
  candidate_long,
  Analysis,
  Gene_ID
)


fwrite(
  candidate_long,
  file.path(
    OUTDIR,
    "GO_candidate_gene_sets_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


for (
  a in unique(
    candidate_long$Analysis
  )
) {

  write_gene_list(
    candidate_long[
      Analysis == a,
      Gene_ID
    ],
    file.path(
      OUTDIR,
      paste0(
        "candidate_genes_",
        a,
        ".txt"
      )
    )
  )
}


# ============================================================
# 9. CANDIDATE SYMBOL DICTIONARY
# ============================================================

dict_parts <- list()


if (
  "Gene_symbol" %in%
    names(xtx)
) {

  dict_parts[[
    length(dict_parts) + 1L
  ]] <- xtx[
    ,
    .(
      Gene_ID =
        clean_gene_id(
          Gene_ID
        ),

      Gene_symbol =
        clean_na(
          Gene_symbol
        )
    )
  ]
}


if (
  "Gene_symbol" %in%
    names(c2)
) {

  dict_parts[[
    length(dict_parts) + 1L
  ]] <- c2[
    ,
    .(
      Gene_ID =
        clean_gene_id(
          Gene_ID
        ),

      Gene_symbol =
        clean_na(
          Gene_symbol
        )
    )
  ]
}


if (
  "Gene_symbol" %in%
    names(gea)
) {

  dict_parts[[
    length(dict_parts) + 1L
  ]] <- gea[
    ,
    .(
      Gene_ID =
        clean_gene_id(
          Gene_ID
        ),

      Gene_symbol =
        clean_na(
          Gene_symbol
        )
    )
  ]
}


if (
  length(dict_parts) > 0L
) {

  gene_dict <- rbindlist(
    dict_parts,
    use.names = TRUE,
    fill = TRUE
  )


  gene_dict <- gene_dict[
    valid_gene_id(
      Gene_ID
    )
  ]


  gene_dict[
    ,
    Biological_symbol :=
      !is.na(
        Gene_symbol
      ) &
      !grepl(
        "^Ccac_g[0-9]+$",
        Gene_symbol
      )
  ]


  setorder(
    gene_dict,
    Gene_ID,
    -Biological_symbol
  )


  gene_dict <- gene_dict[
    ,
    .SD[1L],
    by =
      Gene_ID
  ]


  gene_dict[
    ,
    Display_gene :=
      fifelse(
        Biological_symbol,
        Gene_symbol,
        Gene_ID
      )
  ]


  fwrite(
    gene_dict,
    file.path(
      OUTDIR,
      "GO_candidate_gene_dictionary_FINAL.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
  )
}


# ============================================================
# 10. QC
# ============================================================

background_qc <- data.table(

  N_BayPass_SNPs =
    nrow(snps),

  N_BayPass_SNPs_with_gene_within5kb =
    uniqueN(
      snp_gene_map$BayPass_row
    ),

  Percent_BayPass_SNPs_with_gene_within5kb =
    round(
      100 *
        uniqueN(
          snp_gene_map$BayPass_row
        ) /
        nrow(snps),
      2
    ),

  N_accessible_background_genes =
    length(
      background_genes
    ),

  N_background_genes_with_any_GO =
    sum(
      background_genes %in%
        genes_with_go
    ),

  Percent_background_genes_with_any_GO =
    round(
      100 *
        sum(
          background_genes %in%
            genes_with_go
        ) /
        length(
          background_genes
        ),
      2
    )
)


candidate_qc <- candidate_long[
  ,
  .(
    N_candidate_genes =
      .N,

    N_candidate_genes_with_any_GO =
      sum(
        Has_any_GO
      ),

    Percent_candidate_genes_with_any_GO =
      round(
        100 *
          mean(
            Has_any_GO
          ),
        2
      )
  ),
  by =
    Analysis
]


preferred_order <- c(
  "XtX",
  "C2_Host",
  "GEA_pooled",
  paste0(
    "GEA_",
    PRIMARY_COVARIATES
  )
)


candidate_qc[
  ,
  ord :=
    match(
      Analysis,
      preferred_order
    )
]


setorder(
  candidate_qc,
  ord
)


candidate_qc[
  ,
  ord :=
    NULL
]


fwrite(
  background_qc,
  file.path(
    OUTDIR,
    "QC_GO_background.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


fwrite(
  candidate_qc,
  file.path(
    OUTDIR,
    "QC_GO_candidate_gene_sets.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# FINAL REPORT
# ============================================================

all_selection_union <- unique(
  c(
    xtx_genes,
    c2_genes,
    gea_pooled_genes
  )
)


cat(
  "\n[40A] OK\n"
)

cat(
  "  BayPass SNPs screened: ",
  nrow(snps),
  "\n",
  sep = ""
)

cat(
  "  Accessible background genes: ",
  length(background_genes),
  "\n",
  sep = ""
)

cat(
  "  Background genes with GO: ",
  sum(
    background_genes %in%
      genes_with_go
  ),
  "\n",
  sep = ""
)

cat(
  "  XtX candidate genes: ",
  length(xtx_genes),
  "\n",
  sep = ""
)

cat(
  "  C2 candidate genes: ",
  length(c2_genes),
  "\n",
  sep = ""
)

cat(
  "  GEA pooled candidate genes: ",
  length(gea_pooled_genes),
  "\n",
  sep = ""
)

cat(
  "  All-selection union: ",
  length(all_selection_union),
  "\n",
  sep = ""
)

cat(
  "\nCandidate GO coverage:\n"
)

print(candidate_qc)

cat(
  "\nOutput directory:\n",
  OUTDIR,
  "\n",
  sep = ""
)
