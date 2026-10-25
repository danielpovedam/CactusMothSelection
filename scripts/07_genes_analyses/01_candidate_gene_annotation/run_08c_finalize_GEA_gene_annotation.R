#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

############################################################
## FINAL GEA MANUSCRIPT TABLE
##
## Dataset: all SNPs
## Covariates:
##   TempPC1
##   PrecipPC1
##   AridityIndex
##   Elevation
##   Latitude
##   Longitude
##
## Candidate definition:
##   median BF > 20 dB across 3 BayPass runs
##   BF>20 SNPs <= 1 kb apart grouped into loci
##
## Description hierarchy identical to XtX/C2:
##   1. Swiss-Prot description
##   2. GFF predicted product
##   3. InterPro annotation
##   4. Hypothetical / poorly characterized protein
############################################################


############################################################
## PATHS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

FUNC_BASE <- paste0(
  "/home/dpoveda/micetral/Population_snps/",
  "Ccac_functional-annotation"
)

GEA_FILE <- file.path(
  BASE,
  "selection_scan",
  "FINAL_GEA_BF20_1kb",
  "Table_GEA_BF20_candidate_loci_1kb_PRIMARY_FINAL.tsv"
)

GFF_FILE <- file.path(
  FUNC_BASE,
  "02_agat",
  "Ccac",
  "out_abril4",
  "Ccac.annotations.corrected.gff3"
)

INTERPRO_FILE <- file.path(
  FUNC_BASE,
  "02_agat",
  "Ccac",
  "out_abril4",
  "InterPro.txt"
)

DIAMOND_FILE <- file.path(
  FUNC_BASE,
  "01_mapping2db",
  "diamond_uniprot",
  "Ccac.uniprot.blast.ultra-sensitive"
)

SWISSPROT_FASTA <- file.path(
  FUNC_BASE,
  "01_mapping2db",
  "diamond_uniprot",
  "uniprot_sprot.fasta"
)

OUTDIR <- file.path(
  BASE,
  "selection_scan",
  "FINAL_manuscript_tables"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

OUT_ALL <- file.path(
  OUTDIR,
  "Table_GEA_primary_allSNP_FINAL.tsv"
)

OUT_GENE <- file.path(
  OUTDIR,
  "Table_GEA_primary_gene_associated_allSNP_FINAL.tsv"
)

OUT_ANNOT <- file.path(
  OUTDIR,
  "Table_GEA_primary_gene_annotation_dictionary.tsv"
)


############################################################
## FINAL COVARIATES
############################################################

FINAL_COVARIATES <- c(
  "TempPC1",
  "PrecipPC1",
  "AridityIndex",
  "Elevation",
  "Latitude",
  "Longitude"
)


############################################################
## HELPERS
############################################################

present <- function(x) {

  !is.na(x) &
    trimws(as.character(x)) != ""
}


gene_from_transcript <- function(x) {

  sub(
    "\\.t[0-9]+.*$",
    "",
    as.character(x)
  )
}


extract_attr <- function(x, key) {

  pattern <- paste0(
    "(^|;)",
    key,
    "=([^;]+)"
  )

  m <- regexec(
    pattern,
    as.character(x),
    perl = TRUE
  )

  z <- regmatches(
    as.character(x),
    m
  )

  out <- vapply(
    z,
    function(v) {

      if (length(v) >= 3L) {
        v[3L]
      } else {
        NA_character_
      }

    },
    character(1)
  )

  out[
    out == ""
  ] <- NA_character_

  idx <- !is.na(out)

  if (any(idx)) {

    out[idx] <- vapply(
      out[idx],
      utils::URLdecode,
      character(1)
    )
  }

  out
}


first_valid <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
    trimws(x) != ""
  ]

  if (length(x) == 0L) {

    return(
      NA_character_
    )
  }

  x[1L]
}


collapse_unique <- function(x) {

  x <- unique(
    as.character(x)[
      !is.na(x) &
      trimws(as.character(x)) != ""
    ]
  )

  if (length(x) == 0L) {

    return(
      NA_character_
    )
  }

  paste(
    x,
    collapse = "; "
  )
}


clean_symbol <- function(x, gene_id) {

  x <- as.character(x)
  gene_id <- as.character(gene_id)

  out <- x

  idx <- !is.na(out) &
    grepl(
      ",",
      out,
      fixed = TRUE
    )

  out[idx] <- sub(
    "^.*?,",
    "",
    out[idx]
  )

  out[
    !is.na(out) &
    !is.na(gene_id) &
    out == gene_id
  ] <- NA_character_

  out[
    !is.na(out) &
    trimws(out) == ""
  ] <- NA_character_

  out
}


############################################################
## CHECK FILES
############################################################

required_files <- c(
  GEA_FILE,
  GFF_FILE,
  INTERPRO_FILE,
  DIAMOND_FILE,
  SWISSPROT_FASTA
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files) > 0L) {

  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


############################################################
## LOAD FINAL SIX-COVARIATE GEA TABLE
############################################################

cat(
  "\n=============================================\n",
  "LOADING FINAL SIX-COVARIATE GEA RESULTS\n",
  "=============================================\n",
  sep = ""
)

gea <- fread(
  GEA_FILE,
  showProgress = FALSE
)


############################################################
## HARD FILTER TO FINAL COVARIATES
############################################################

gea <- gea[
  Covariate %in%
    FINAL_COVARIATES
]


unexpected <- setdiff(
  unique(
    gea$Covariate
  ),
  FINAL_COVARIATES
)

if (length(unexpected) > 0L) {

  stop(
    "Unexpected covariates: ",
    paste(
      unexpected,
      collapse = ", "
    )
  )
}


if (any(
  gea$Covariate %in%
    c(
      "BIO1",
      "BIO12"
    )
)) {

  stop(
    "BIO1/BIO12 unexpectedly remain in final GEA table."
  )
}


cat(
  "Final GEA loci: ",
  nrow(gea),
  "\n",
  sep = ""
)

cat(
  "Covariates: ",
  paste(
    unique(
      gea$Covariate
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


if (nrow(gea) != 167L) {

  warning(
    "Expected 167 candidate loci; found ",
    nrow(gea)
  )
}


############################################################
## REQUIRED COLUMNS
############################################################

required_cols <- c(
  "Covariate",
  "GEA_locus_ID",
  "Scaffold",
  "Region_start",
  "Region_end",
  "Region_size_bp",
  "N_BF20_SNPs",
  "Peak_position",
  "Peak_BFdb",
  "Peak_context",
  "Gene_ID",
  "Gene_symbol",
  "Distance_to_gene_bp"
)

missing_cols <- setdiff(
  required_cols,
  names(gea)
)

if (length(missing_cols) > 0L) {

  stop(
    "Missing GEA columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}


############################################################
## CANDIDATE GENES
############################################################

candidate_genes <- unique(
  gea[
    present(
      Gene_ID
    ),
    Gene_ID
  ]
)

cat(
  "Unique candidate gene IDs: ",
  length(
    candidate_genes
  ),
  "\n",
  sep = ""
)


############################################################
## GFF TRANSCRIPT ANNOTATION
############################################################

cat(
  "Reading transcript annotations from GFF...\n"
)

gff_cmd <- sprintf(
  "awk -F '\\t' '$3==\"transcript\"' %s",
  shQuote(
    GFF_FILE
  )
)

gff <- fread(
  cmd = gff_cmd,
  sep = "\t",
  header = FALSE,
  fill = TRUE,
  quote = "",
  showProgress = FALSE
)

if (ncol(gff) < 9L) {

  stop(
    "Failed to read transcript records from GFF."
  )
}

setnames(
  gff,
  1:9,
  c(
    "seqid",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
  )
)


gff[
  ,
  Transcript_ID :=
    extract_attr(
      attributes,
      "ID"
    )
]

gff[
  ,
  Gene_ID :=
    extract_attr(
      attributes,
      "Parent"
    )
]

gff[
  ,
  Product :=
    extract_attr(
      attributes,
      "product"
    )
]


gff <- gff[
  Gene_ID %in%
    candidate_genes
]


gff_gene <- gff[
  ,
  .(
    Transcript_IDs =
      collapse_unique(
        Transcript_ID
      ),

    Product =
      first_valid(
        Product
      )
  ),
  by = Gene_ID
]


############################################################
## INTERPRO
############################################################

cat(
  "Reading InterPro annotations...\n"
)

ipr <- fread(
  INTERPRO_FILE,
  sep = "\t",
  header = FALSE,
  fill = TRUE,
  quote = "",
  showProgress = FALSE
)

if (ncol(ipr) < 2L) {

  stop(
    "InterPro file has fewer than two columns."
  )
}

ipr <- ipr[
  ,
  1:2,
  with = FALSE
]

setnames(
  ipr,
  c(
    "Transcript_ID",
    "InterPro"
  )
)

ipr[
  ,
  Gene_ID :=
    gene_from_transcript(
      Transcript_ID
    )
]

ipr <- ipr[
  Gene_ID %in%
    candidate_genes
]

ipr_gene <- ipr[
  ,
  .(
    InterPro =
      collapse_unique(
        InterPro
      )
  ),
  by = Gene_ID
]


############################################################
## DIAMOND / SWISS-PROT
############################################################

cat(
  "Reading DIAMOND Swiss-Prot hits...\n"
)

diamond <- fread(
  DIAMOND_FILE,
  sep = "\t",
  header = FALSE,
  fill = TRUE,
  quote = "",
  showProgress = FALSE
)

if (ncol(diamond) < 12L) {

  stop(
    "DIAMOND file has fewer than 12 columns."
  )
}

diamond <- diamond[
  ,
  1:12,
  with = FALSE
]

setnames(
  diamond,
  c(
    "Transcript_ID",
    "SwissProt_subject",
    "Percent_identity",
    "Alignment_length",
    "Mismatch",
    "Gapopen",
    "Qstart",
    "Qend",
    "Sstart",
    "Send",
    "Evalue",
    "Bitscore"
  )
)

diamond[
  ,
  Gene_ID :=
    gene_from_transcript(
      Transcript_ID
    )
]

diamond <- diamond[
  Gene_ID %in%
    candidate_genes
]

diamond[
  ,
  `:=`(
    Percent_identity =
      suppressWarnings(
        as.numeric(
          Percent_identity
        )
      ),

    Evalue =
      suppressWarnings(
        as.numeric(
          Evalue
        )
      ),

    Bitscore =
      suppressWarnings(
        as.numeric(
          Bitscore
        )
      )
  )
]


############################################################
## SAME BEST-HIT LOGIC AS PREVIOUS XTX/C2 SCRIPT
############################################################

setorder(
  diamond,
  Gene_ID,
  -Bitscore,
  Evalue,
  -Percent_identity
)

best_diamond <- diamond[
  ,
  .SD[1L],
  by = Gene_ID
]


############################################################
## SWISS-PROT FASTA HEADERS
############################################################

cat(
  "Parsing Swiss-Prot headers...\n"
)

con <- pipe(
  paste(
    "grep '^>'",
    shQuote(
      SWISSPROT_FASTA
    )
  ),
  open = "r"
)

headers <- readLines(
  con,
  warn = FALSE
)

close(
  con
)

sp <- data.table(
  Header = headers
)

sp[
  ,
  SwissProt_subject :=
    sub(
      "^>(\\S+).*$",
      "\\1",
      Header
    )
]

needed_subjects <- unique(
  best_diamond$SwissProt_subject
)

sp <- sp[
  SwissProt_subject %in%
    needed_subjects
]


sp[
  ,
  SwissProt_description :=
    sub(
      "^>\\S+\\s+",
      "",
      Header
    )
]

sp[
  ,
  SwissProt_description :=
    sub(
      "\\s+OS=.*$",
      "",
      SwissProt_description
    )
]


sp[
  ,
  SwissProt_gene_name :=
    fifelse(
      grepl(
        " GN=",
        Header
      ),

      sub(
        ".* GN=([^ ]+).*",
        "\\1",
        Header
      ),

      NA_character_
    )
]


best_diamond <- merge(
  best_diamond,
  sp[
    ,
    .(
      SwissProt_subject,
      SwissProt_gene_name,
      SwissProt_description
    )
  ],
  by = "SwissProt_subject",
  all.x = TRUE,
  sort = FALSE
)


diamond_gene <- best_diamond[
  ,
  .(
    Gene_ID,
    SwissProt_gene_name,
    SwissProt_description,
    Percent_identity,
    Evalue,
    Bitscore
  )
]


############################################################
## BUILD ANNOTATION DICTIONARY
############################################################

ann <- data.table(
  Gene_ID =
    candidate_genes
)

ann <- merge(
  ann,
  gff_gene,
  by = "Gene_ID",
  all.x = TRUE,
  sort = FALSE
)

ann <- merge(
  ann,
  ipr_gene,
  by = "Gene_ID",
  all.x = TRUE,
  sort = FALSE
)

ann <- merge(
  ann,
  diamond_gene,
  by = "Gene_ID",
  all.x = TRUE,
  sort = FALSE
)


############################################################
## DESCRIPTION
##
## IDENTICAL HIERARCHY TO XTX/C2
############################################################

ann[
  ,
  Biological_description :=
    fifelse(

      present(
        SwissProt_description
      ),

      SwissProt_description,

      fifelse(

        present(
          Product
        ),

        Product,

        fifelse(

          present(
            InterPro
          ),

          paste0(
            "Predicted protein; InterPro annotation: ",
            InterPro
          ),

          "Hypothetical / poorly characterized protein"
        )
      )
    )
]


############################################################
## SAVE ANNOTATION DICTIONARY
############################################################

fwrite(
  ann,
  OUT_ANNOT,
  sep = "\t",
  quote = FALSE,
  na = ""
)


############################################################
## MERGE WITH GEA
############################################################

gea[
  ,
  Existing_symbol :=
    clean_symbol(
      Gene_symbol,
      Gene_ID
    )
]

gea <- merge(
  gea,
  ann,
  by = "Gene_ID",
  all.x = TRUE,
  sort = FALSE
)


############################################################
## FINAL GENE SYMBOL
############################################################

gea[
  ,
  Final_symbol :=
    fifelse(

      present(
        Existing_symbol
      ),

      Existing_symbol,

      fifelse(

        present(
          SwissProt_gene_name
        ),

        SwissProt_gene_name,

        NA_character_
      )
    )
]


############################################################
## DISTAL LOCI HAVE NO DESCRIPTION
############################################################

gea[
  !present(
    Gene_ID
  ),
  Biological_description :=
    NA_character_
]


############################################################
## FINAL SIMPLE MANUSCRIPT TABLE
############################################################

final <- gea[
  ,
  .(

    Statistic =
      "BF",

    Covariate,

    Scaffold,

    Window_start =
      Region_start,

    Window_end =
      Region_end,

    Window_size_bp =
      Region_size_bp,

    N_SNPs =
      N_BF20_SNPs,

    BF_peak_position =
      Peak_position,

    BF_peak_value =
      round(
        Peak_BFdb,
        2
      ),

    Peak_context,

    Gene_ID,

    Gene_symbol =
      Final_symbol,

    Distance_to_gene_bp,

    Description =
      Biological_description
  )
]


############################################################
## ORDER
############################################################

final[
  ,
  Cov_order :=
    match(
      Covariate,
      FINAL_COVARIATES
    )
]

setorder(
  final,
  Cov_order,
  Scaffold,
  Window_start
)

final[
  ,
  Cov_order :=
    NULL
]


############################################################
## WRITE ALL 167 LOCI
############################################################

fwrite(
  final,
  OUT_ALL,
  sep = "\t",
  quote = FALSE,
  na = ""
)


############################################################
## WRITE GENE-ASSOCIATED LOCI ONLY
############################################################

gene_final <- final[
  present(
    Gene_ID
  )
]

fwrite(
  gene_final,
  OUT_GENE,
  sep = "\t",
  quote = FALSE,
  na = ""
)


############################################################
## QC
############################################################

cat(
  "\n=============================================\n",
  "FINAL PRIMARY GEA MANUSCRIPT TABLE QC\n",
  "=============================================\n",
  sep = ""
)

cat(
  "Total candidate loci: ",
  nrow(final),
  "\n",
  sep = ""
)

cat(
  "BF>20 SNPs represented: ",
  sum(
    final$N_SNPs
  ),
  "\n",
  sep = ""
)

cat(
  "Gene-associated loci: ",
  nrow(
    gene_final
  ),
  "\n",
  sep = ""
)

cat(
  "Unique genes: ",
  uniqueN(
    gene_final$Gene_ID
  ),
  "\n",
  sep = ""
)

cat(
  "Gene-associated loci with description: ",
  sum(
    present(
      gene_final$Description
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Gene-associated loci missing description: ",
  sum(
    !present(
      gene_final$Description
    )
  ),
  "\n",
  sep = ""
)


cat(
  "\nCounts by covariate:\n"
)

print(
  final[
    ,
    .(
      N_loci = .N,
      N_SNPs = sum(N_SNPs),
      N_gene_associated =
        sum(
          present(
            Gene_ID
          )
        )
    ),
    by = Covariate
  ][
    match(
      FINAL_COVARIATES,
      Covariate
    )
  ]
)


cat(
  "\nDescription examples:\n"
)

print(
  head(
    gene_final[
      ,
      .(
        Covariate,
        Scaffold,
        BF_peak_value,
        Gene_ID,
        Gene_symbol,
        Peak_context,
        Description
      )
    ],
    30
  )
)


cat(
  "\nOutputs:\n",
  OUT_ALL,
  "\n",
  OUT_GENE,
  "\n",
  OUT_ANNOT,
  "\n",
  sep = ""
)

cat(
  "\nDONE\n"
)

