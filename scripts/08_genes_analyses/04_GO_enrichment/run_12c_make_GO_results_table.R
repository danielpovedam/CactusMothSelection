#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

GO_DIR <- file.path(
  BASE,
  "selection_scan",
  "GO_ENRICHMENT_FINAL_gt1Mb"
)

INFILE <- file.path(
  GO_DIR,
  "Table_GO_enrichment_QC_summary.tsv"
)

OUTFILE <- file.path(
  GO_DIR,
  "Table_S12_GO_enrichment_summary_FINAL_gt1Mb.tsv"
)

x <- fread(
  INFILE,
  showProgress = FALSE
)


required <- c(
  "Analysis",
  "Ontology",
  "Candidate_genes_input",
  "Candidate_genes_mapped",
  "Best_weight01_P",
  "Best_weight01_FDR",
  "N_weight01_FDR05"
)


missing <- setdiff(
  required,
  names(x)
)


if (
  length(missing) > 0L
) {

  stop(
    "Missing columns: ",
    paste(
      missing,
      collapse = ", "
    )
  )
}


analysis_labels <- c(

  XtX =
    "XtX",

  C2_Host =
    "C2 host-use contrast",

  GEA_pooled =
    "GEA pooled",

  All_selection =
    "All selection candidates",

  GEA_TempPC1 =
    "TemperaturePC1",

  GEA_PrecipPC1 =
    "PrecipitationPC1",

  GEA_AridityIndex =
    "Aridity index",

  GEA_Elevation =
    "Elevation",

  GEA_Latitude =
    "Latitude",

  GEA_Longitude =
    "Longitude"
)


x[
  ,
  Candidate_set :=
    analysis_labels[
      Analysis
    ]
]


x[
  is.na(
    Candidate_set
  ),
  Candidate_set :=
    Analysis
]


preferred <- c(
  "All_selection",
  "XtX",
  "C2_Host",
  "GEA_pooled",
  "GEA_TempPC1",
  "GEA_PrecipPC1",
  "GEA_AridityIndex",
  "GEA_Elevation",
  "GEA_Latitude",
  "GEA_Longitude"
)


ontology_order <- c(
  "BP",
  "MF",
  "CC"
)


x[
  ,
  `:=`(
    Analysis_order =
      match(
        Analysis,
        preferred
      ),

    Ontology_order =
      match(
        Ontology,
        ontology_order
      )
  )
]


setorder(
  x,
  Analysis_order,
  Ontology_order
)


out <- x[
  ,
  .(
    Candidate_set,
    Ontology,
    Candidate_genes =
      Candidate_genes_input,
    GO_mapped_genes =
      Candidate_genes_mapped,
    Best_P =
      Best_weight01_P,
    Minimum_FDR =
      Best_weight01_FDR,
    GO_terms_FDR_lt_0.05 =
      N_weight01_FDR05
  )
]


fwrite(
  out,
  OUTFILE,
  sep = "\t",
  quote = FALSE,
  na = "-"
)


cat(
  "\nFINAL Table S12\n\n"
)

print(out)


cat(
  "\nWritten to:\n",
  OUTFILE,
  "\n",
  sep = ""
)
