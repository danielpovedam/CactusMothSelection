#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(topGO)
  library(GO.db)
  library(AnnotationDbi)
})

############################################################
## 40B_FINAL_gt1Mb_topGO_Ccactorum.R
##
## Formal GO enrichment using final >1 Mb XtX/C2 sets.
##
## Formal test:
##   Fisher's exact test
##   topGO weight01 algorithm
##
## Multiple testing:
##   BH-FDR within each analysis x ontology
############################################################


# ============================================================
# PATHS
# ============================================================

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

OUTDIR <- file.path(
  BASE,
  "selection_scan",
  "GO_ENRICHMENT_FINAL_gt1Mb"
)

BACKGROUND_FILE <- file.path(
  OUTDIR,
  "GO_background_gene_list.txt"
)

CANDIDATE_FILE <- file.path(
  OUTDIR,
  "GO_candidate_gene_sets_FINAL.tsv"
)

GENE2GO_FILE <- file.path(
  OUTDIR,
  "Ccactorum_gene2GO_FINAL.tsv"
)

DICT_FILE <- file.path(
  OUTDIR,
  "GO_candidate_gene_dictionary_FINAL.tsv"
)


# ============================================================
# SETTINGS
# ============================================================

ONTOLOGIES <- c(
  "BP",
  "MF",
  "CC"
)

PRIMARY_ANALYSES <- c(
  "XtX",
  "C2_Host",
  "GEA_pooled",
  "All_selection"
)

GEA_ANALYSES <- c(
  "GEA_TempPC1",
  "GEA_PrecipPC1",
  "GEA_AridityIndex",
  "GEA_Elevation",
  "GEA_Latitude",
  "GEA_Longitude"
)

EXPECTED_ANALYSES <- c(
  PRIMARY_ANALYSES,
  GEA_ANALYSES
)

NODE_SIZE <- 10L

FDR_CUTOFF <- 0.05

NOMINAL_CUTOFF <- 0.01

TOP_N <- 30L

MIN_MAPPED_CANDIDATES_FLAG <- 5L


# ============================================================
# HELPERS
# ============================================================

check_file <- function(path) {

  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0
  ) {

    stop(
      "Missing or empty required file:\n",
      path
    )
  }
}


safe_min <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (
    length(x) == 0L
  ) {

    return(
      NA_real_
    )
  }

  min(x)
}


for (
  f in c(
    BACKGROUND_FILE,
    CANDIDATE_FILE,
    GENE2GO_FILE,
    DICT_FILE
  )
) {

  check_file(f)
}


cat(
  "\n============================================================\n",
  "40B FINAL gt1Mb topGO ENRICHMENT\n",
  "============================================================\n",
  sep = ""
)


# ============================================================
# BACKGROUND
# ============================================================

background_dt <- fread(
  BACKGROUND_FILE,
  header = FALSE,
  col.names = "Gene_ID",
  showProgress = FALSE
)


background_dt[
  ,
  Gene_ID :=
    trimws(
      as.character(
        Gene_ID
      )
    )
]


background_dt <- unique(
  background_dt[
    !is.na(Gene_ID) &
      Gene_ID != ""
  ]
)


background_genes <- sort(
  unique(
    background_dt$Gene_ID
  )
)


if (
  length(background_genes) == 0L
) {

  stop(
    "Background gene set is empty."
  )
}


# ============================================================
# CANDIDATE SETS
# ============================================================

candidate_dt <- fread(
  CANDIDATE_FILE,
  header = TRUE,
  showProgress = FALSE
)


if (
  !all(
    c(
      "Analysis",
      "Gene_ID"
    ) %in%
      names(candidate_dt)
  )
) {

  stop(
    "Candidate file must contain Analysis and Gene_ID."
  )
}


candidate_dt[
  ,
  `:=`(
    Analysis =
      as.character(
        Analysis
      ),

    Gene_ID =
      trimws(
        as.character(
          Gene_ID
        )
      )
  )
]


candidate_dt <- unique(
  candidate_dt[
    !is.na(Analysis) &
      Analysis != "" &
      !is.na(Gene_ID) &
      Gene_ID != "",
    .(
      Analysis,
      Gene_ID
    )
  ]
)


# ============================================================
# GLOBAL UNION
# ============================================================

all_selection_genes <- sort(
  unique(
    candidate_dt[
      Analysis %in%
        c(
          "XtX",
          "C2_Host",
          "GEA_pooled"
        ),
      Gene_ID
    ]
  )
)


candidate_dt <- rbind(
  candidate_dt,
  data.table(
    Analysis =
      "All_selection",
    Gene_ID =
      all_selection_genes
  ),
  use.names = TRUE,
  fill = TRUE
)


candidate_dt <- unique(
  candidate_dt,
  by = c(
    "Analysis",
    "Gene_ID"
  )
)


missing_analyses <- setdiff(
  EXPECTED_ANALYSES,
  unique(
    candidate_dt$Analysis
  )
)


if (
  length(missing_analyses) > 0L
) {

  stop(
    "Missing candidate sets: ",
    paste(
      missing_analyses,
      collapse = ", "
    )
  )
}


cat(
  "All-selection union: ",
  length(
    all_selection_genes
  ),
  " genes\n",
  sep = ""
)


# ============================================================
# GENE2GO
# ============================================================

go_long <- fread(
  GENE2GO_FILE,
  header = TRUE,
  showProgress = FALSE
)


if (
  !all(
    c(
      "Gene_ID",
      "GO_ID"
    ) %in%
      names(go_long)
  )
) {

  stop(
    "Gene2GO file must contain Gene_ID and GO_ID."
  )
}


go_long[
  ,
  `:=`(
    Gene_ID =
      trimws(
        as.character(
          Gene_ID
        )
      ),

    GO_ID =
      trimws(
        as.character(
          GO_ID
        )
      )
  )
]


go_long <- unique(
  go_long[
    !is.na(Gene_ID) &
      Gene_ID != "" &
      grepl(
        "^GO:[0-9]{7}$",
        GO_ID
      ),
    .(
      Gene_ID,
      GO_ID
    )
  ]
)


# ============================================================
# DICTIONARY
# ============================================================

gene_dictionary <- fread(
  DICT_FILE,
  header = TRUE,
  showProgress = FALSE
)


if (
  !"Gene_ID" %in%
    names(gene_dictionary)
) {

  stop(
    "Gene dictionary lacks Gene_ID."
  )
}


if (
  !"Display_gene" %in%
    names(gene_dictionary)
) {

  gene_dictionary[
    ,
    Display_gene :=
      Gene_ID
  ]
}


gene_dictionary[
  ,
  `:=`(
    Gene_ID =
      trimws(
        as.character(
          Gene_ID
        )
      ),

    Display_gene =
      trimws(
        as.character(
          Display_gene
        )
      )
  )
]


gene_dictionary <- unique(
  gene_dictionary,
  by =
    "Gene_ID"
)


# ============================================================
# VALIDATE GO IDS AGAINST GO.db
# ============================================================

valid_go_ids <- AnnotationDbi::keys(
  GO.db::GO.db,
  keytype = "GOID"
)


invalid_go <- unique(
  go_long[
    !GO_ID %in%
      valid_go_ids,
    .(
      GO_ID
    )
  ]
)


if (
  nrow(invalid_go) > 0L
) {

  fwrite(
    invalid_go,
    file.path(
      OUTDIR,
      "QC_GO_terms_absent_from_GOdb.tsv"
    ),
    sep = "\t",
    quote = FALSE
  )
}


go_long <- go_long[
  GO_ID %in%
    valid_go_ids
]


if (
  nrow(go_long) == 0L
) {

  stop(
    "No valid GO assignments remain."
  )
}


# ============================================================
# GO LOOKUP
# ============================================================

go_lookup <- AnnotationDbi::select(
  GO.db::GO.db,
  keys =
    unique(
      go_long$GO_ID
    ),
  keytype =
    "GOID",
  columns =
    c(
      "TERM",
      "ONTOLOGY"
    )
)


go_lookup <- as.data.table(
  go_lookup
)


go_lookup <- unique(
  go_lookup[
    !is.na(GOID) &
      !is.na(ONTOLOGY),
    .(
      GO_ID =
        GOID,
      Term =
        TERM,
      Ontology =
        ONTOLOGY
    )
  ],
  by =
    "GO_ID"
)


go_long <- merge(
  go_long,
  go_lookup,
  by =
    "GO_ID",
  all.x =
    TRUE,
  sort =
    FALSE
)


go_long <- go_long[
  !is.na(Ontology)
]


gene2GO <- split(
  go_long$GO_ID,
  go_long$Gene_ID
)


gene2GO <- lapply(
  gene2GO,
  unique
)


# ============================================================
# LABEL FUNCTION
# ============================================================

get_gene_label <- function(gene_ids) {

  labels <- gene_dictionary[
    match(
      gene_ids,
      Gene_ID
    ),
    Display_gene
  ]


  bad <- is.na(labels) |
    labels == ""


  labels[bad] <-
    gene_ids[bad]


  labels
}


# ============================================================
# RESULT TABLE
# ============================================================

make_result_table <- function(
  GOdata,
  result_object,
  analysis_name,
  ontology_name,
  algorithm_name,
  candidate_genes
) {

  score_vector <- score(
    result_object
  )


  if (
    length(score_vector) == 0L
  ) {

    return(
      data.table()
    )
  }


  result_dt <- data.table(
    GO_ID =
      names(
        score_vector
      ),
    p_value =
      as.numeric(
        score_vector
      )
  )


  result_dt <- merge(
    result_dt,
    go_lookup[
      ,
      .(
        GO_ID,
        Term,
        Ontology
      )
    ],
    by =
      "GO_ID",
    all.x =
      TRUE,
    sort =
      FALSE
  )


  term_stats <- as.data.table(
    termStat(
      GOdata,
      result_dt$GO_ID
    ),
    keep.rownames =
      "GO_ID"
  )


  result_dt <- merge(
    result_dt,
    term_stats,
    by =
      "GO_ID",
    all.x =
      TRUE,
    sort =
      FALSE
  )


  result_dt[
    ,
    FDR :=
      p.adjust(
        p_value,
        method = "BH"
      )
  ]


  result_dt[
    ,
    Fold_enrichment :=
      fifelse(
        !is.na(Expected) &
          Expected > 0,
        Significant /
          Expected,
        NA_real_
      )
  ]


  term_genes <- genesInTerm(
    GOdata,
    result_dt$GO_ID
  )


  contributing <- rbindlist(
    lapply(
      result_dt$GO_ID,
      function(go_id) {

        genes_here <-
          term_genes[[
            go_id
          ]]


        if (
          is.null(
            genes_here
          )
        ) {

          genes_here <-
            character(0)
        }


        candidate_here <-
          intersect(
            genes_here,
            candidate_genes
          )


        labels_here <-
          get_gene_label(
            candidate_here
          )


        data.table(

          GO_ID =
            go_id,

          Candidate_gene_IDs =
            if (
              length(
                candidate_here
              ) == 0L
            ) {
              NA_character_
            } else {
              paste(
                sort(
                  unique(
                    candidate_here
                  )
                ),
                collapse = ";"
              )
            },

          Candidate_gene_labels =
            if (
              length(
                labels_here
              ) == 0L
            ) {
              NA_character_
            } else {
              paste(
                sort(
                  unique(
                    labels_here
                  )
                ),
                collapse = ";"
              )
            }
        )
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )


  result_dt <- merge(
    result_dt,
    contributing,
    by =
      "GO_ID",
    all.x =
      TRUE,
    sort =
      FALSE
  )


  ontology_genes <- unique(
    go_long[
      Ontology ==
        ontology_name,
      Gene_ID
    ]
  )


  mapped_background <- intersect(
    background_genes,
    ontology_genes
  )


  mapped_candidates <- intersect(
    candidate_genes,
    ontology_genes
  )


  result_dt[
    ,
    `:=`(
      Analysis =
        analysis_name,

      Ontology =
        ontology_name,

      Algorithm =
        algorithm_name,

      Candidate_genes_input =
        length(
          candidate_genes
        ),

      Candidate_genes_mapped_to_ontology =
        length(
          mapped_candidates
        ),

      Background_genes_input =
        length(
          background_genes
        ),

      Background_genes_mapped_to_ontology =
        length(
          mapped_background
        )
    )
  ]


  setorder(
    result_dt,
    p_value,
    GO_ID
  )


  result_dt
}


# ============================================================
# RUN ONE ANALYSIS x ONTOLOGY
# ============================================================

run_one <- function(
  analysis_name,
  ontology_name
) {

  candidate_genes <- sort(
    unique(
      candidate_dt[
        Analysis ==
          analysis_name,
        Gene_ID
      ]
    )
  )


  if (
    length(candidate_genes) ==
      0L
  ) {

    stop(
      "No candidate genes for ",
      analysis_name
    )
  }


  outside <- setdiff(
    candidate_genes,
    background_genes
  )


  if (
    length(outside) >
      0L
  ) {

    stop(
      "Candidate genes outside background for ",
      analysis_name,
      ": ",
      paste(
        outside,
        collapse = ", "
      )
    )
  }


  gene_list <- factor(
    as.integer(
      background_genes %in%
        candidate_genes
    )
  )


  names(gene_list) <-
    background_genes


  GOdata <- new(
    "topGOdata",
    ontology =
      ontology_name,
    allGenes =
      gene_list,
    annot =
      annFUN.gene2GO,
    gene2GO =
      gene2GO,
    nodeSize =
      NODE_SIZE
  )


  ontology_genes <- unique(
    go_long[
      Ontology ==
        ontology_name,
      Gene_ID
    ]
  )


  mapped_candidates <-
    length(
      intersect(
        candidate_genes,
        ontology_genes
      )
    )


  mapped_background <-
    length(
      intersect(
        background_genes,
        ontology_genes
      )
    )


  if (
    mapped_candidates ==
      0L ||
    mapped_background ==
      0L
  ) {

    summary_dt <- data.table(

      Analysis =
        analysis_name,

      Ontology =
        ontology_name,

      Candidate_genes_input =
        length(
          candidate_genes
        ),

      Candidate_genes_mapped =
        mapped_candidates,

      Candidate_GO_coverage_percent =
        round(
          100 *
            mapped_candidates /
            length(
              candidate_genes
            ),
          2
        ),

      Background_genes_input =
        length(
          background_genes
        ),

      Background_genes_mapped =
        mapped_background,

      Background_GO_coverage_percent =
        round(
          100 *
            mapped_background /
            length(
              background_genes
            ),
          2
        ),

      Low_power_flag =
        mapped_candidates <
          MIN_MAPPED_CANDIDATES_FLAG,

      N_terms_tested =
        0L,

      N_weight01_FDR05 =
        0L,

      Best_weight01_P =
        NA_real_,

      Best_weight01_FDR =
        NA_real_
    )


    return(
      list(
        weight01 =
          data.table(),
        summary =
          summary_dt
      )
    )
  }


  weight_result <- runTest(
    GOdata,
    algorithm =
      "weight01",
    statistic =
      "fisher"
  )


  weight_dt <- make_result_table(
    GOdata =
      GOdata,
    result_object =
      weight_result,
    analysis_name =
      analysis_name,
    ontology_name =
      ontology_name,
    algorithm_name =
      "weight01",
    candidate_genes =
      candidate_genes
  )


  summary_dt <- data.table(

    Analysis =
      analysis_name,

    Ontology =
      ontology_name,

    Candidate_genes_input =
      length(
        candidate_genes
      ),

    Candidate_genes_mapped =
      mapped_candidates,

    Candidate_GO_coverage_percent =
      round(
        100 *
          mapped_candidates /
          length(
            candidate_genes
          ),
        2
      ),

    Background_genes_input =
      length(
        background_genes
      ),

    Background_genes_mapped =
      mapped_background,

    Background_GO_coverage_percent =
      round(
        100 *
          mapped_background /
          length(
            background_genes
          ),
        2
      ),

    Low_power_flag =
      mapped_candidates <
        MIN_MAPPED_CANDIDATES_FLAG,

    N_terms_tested =
      nrow(
        weight_dt
      ),

    N_weight01_FDR05 =
      sum(
        weight_dt$FDR <
          FDR_CUTOFF,
        na.rm = TRUE
      ),

    Best_weight01_P =
      safe_min(
        weight_dt$p_value
      ),

    Best_weight01_FDR =
      safe_min(
        weight_dt$FDR
      )
  )


  list(
    weight01 =
      weight_dt,
    summary =
      summary_dt
  )
}


# ============================================================
# RUN ALL TESTS
# ============================================================

all_weight <- list()
all_summary <- list()

iw <- 0L
is <- 0L


for (
  analysis_name in
    EXPECTED_ANALYSES
) {

  analysis_dir <- file.path(
    OUTDIR,
    analysis_name
  )


  dir.create(
    analysis_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  for (
    ontology_name in
      ONTOLOGIES
  ) {

    cat(
      "Running ",
      analysis_name,
      " / ",
      ontology_name,
      "\n",
      sep = ""
    )


    res <- run_one(
      analysis_name,
      ontology_name
    )


    weight_dt <-
      res$weight01

    summary_dt <-
      res$summary


    iw <- iw + 1L
    is <- is + 1L


    all_weight[[iw]] <-
      weight_dt

    all_summary[[is]] <-
      summary_dt


    if (
      nrow(weight_dt) >
        0L
    ) {

      fwrite(
        weight_dt,
        file.path(
          analysis_dir,
          paste0(
            analysis_name,
            "_",
            ontology_name,
            "_weight01_ALL.tsv"
          )
        ),
        sep = "\t",
        quote = FALSE,
        na = ""
      )


      fwrite(
        weight_dt[
          FDR <
            FDR_CUTOFF
        ],
        file.path(
          analysis_dir,
          paste0(
            analysis_name,
            "_",
            ontology_name,
            "_weight01_FDR05.tsv"
          )
        ),
        sep = "\t",
        quote = FALSE,
        na = ""
      )


      fwrite(
        head(
          weight_dt,
          TOP_N
        ),
        file.path(
          analysis_dir,
          paste0(
            analysis_name,
            "_",
            ontology_name,
            "_weight01_TOP30.tsv"
          )
        ),
        sep = "\t",
        quote = FALSE,
        na = ""
      )
    }
  }
}


# ============================================================
# COMBINED OUTPUTS
# ============================================================

formal_results <- rbindlist(
  all_weight,
  use.names = TRUE,
  fill = TRUE
)


summary_results <- rbindlist(
  all_summary,
  use.names = TRUE,
  fill = TRUE
)


fwrite(
  formal_results,
  file.path(
    OUTDIR,
    "Table_GO_FORMAL_weight01_ALL_terms.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


formal_sig <- formal_results[
  FDR <
    FDR_CUTOFF
]


fwrite(
  formal_sig,
  file.path(
    OUTDIR,
    "Table_GO_FORMAL_weight01_FDR05.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  summary_results,
  file.path(
    OUTDIR,
    "Table_GO_enrichment_QC_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)


# ============================================================
# FINAL REPORT
# ============================================================

cat(
  "\n============================================================\n",
  "FINAL GO ENRICHMENT COMPLETE\n",
  "============================================================\n",
  sep = ""
)


cat(
  "Background genes: ",
  length(
    background_genes
  ),
  "\n",
  sep = ""
)


cat(
  "All-selection genes: ",
  length(
    all_selection_genes
  ),
  "\n",
  sep = ""
)


cat(
  "Formal FDR <0.05 terms: ",
  nrow(
    formal_sig
  ),
  "\n\n",
  sep = ""
)


print(
  summary_results[
    ,
    .(
      Analysis,
      Ontology,
      Candidate_genes_input,
      Candidate_genes_mapped,
      Best_weight01_P,
      Best_weight01_FDR,
      N_weight01_FDR05
    )
  ]
)


cat(
  "\nQC table:\n",
  file.path(
    OUTDIR,
    "Table_GO_enrichment_QC_summary.tsv"
  ),
  "\n",
  sep = ""
)
