#!/usr/bin/env Rscript

############################################################
## 46_build_FINAL_gt1Mb_XtX_C2_tables_and_summary.R
##
## Cactoblastis cactorum
##
## PURPOSE
##
## Rebuild final annotated XtX and C2 candidate tables after
## restricting local-score analyses to scaffolds >1 Mb.
##
## IMPORTANT:
## The formal >1 Mb rerun showed:
##
## XtX:
##   12 identical retained regions
##   4 removed short-scaffold regions
##
## C2:
##   17 identical retained regions
##   6 removed short-scaffold regions
##
## Therefore retained rows can safely be filtered directly
## from the existing fully annotated final tables.
##
## This script:
##
## 1. reads true scaffold lengths from FASTA .fai
## 2. filters old annotated tables to >1 Mb scaffolds
## 3. preserves original Region_ID values
## 4. recalculates:
##      - candidate counts
##      - cumulative region span
##      - total regional SNP memberships
##      - POD99.9 support
## 5. recalculates XtX-C2 physical overlap
## 6. recalculates shared peak-associated genes
## 7. exports revised publication-ready tables
##
## OLD FILES ARE NOT OVERWRITTEN.
############################################################


suppressPackageStartupMessages({
  library(data.table)
})


############################################################
## 0. PATHS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps"

BAYPASS <- file.path(
  BASE,
  "baypass"
)


FAI_FILE <- file.path(
  BASE,
  "Ccactorum.fcs_gx_2026_FINAL.fa.fai"
)


OLD_XTX <- file.path(
  BAYPASS,
  "selection_scan",
  "FINAL_manuscript_tables",
  "Table_XtX_candidate_regions_allSNP_FINAL.tsv"
)


OLD_C2 <- file.path(
  BAYPASS,
  "selection_scan",
  "FINAL_manuscript_tables",
  "Table_C2_Host_candidate_regions_allSNP_FINAL.tsv"
)


NEW_WINDOWS_XTX <- file.path(
  BAYPASS,
  "selection_scan",
  "local_score_gt1Mb",
  "all_snps",
  "XtX",
  "XtX_candidate_windows_gt1Mb.tsv"
)


NEW_WINDOWS_C2 <- file.path(
  BAYPASS,
  "selection_scan",
  "local_score_gt1Mb",
  "all_snps",
  "C2_Host",
  "C2_Host_candidate_windows_gt1Mb.tsv"
)


OUTDIR <- file.path(
  BAYPASS,
  "selection_scan",
  "FINAL_gt1Mb"
)


dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 1. SETTINGS
############################################################

MIN_LENGTH <- 1000000L


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


clean_text <- function(x) {

  x <- as.character(x)

  x[
    is.na(x) |
      x == "" |
      x == "NA" |
      x == "." |
      x == "-"
  ] <- NA_character_

  x
}


informative_gene <- function(symbol, gene_id) {

  symbol <- clean_text(symbol)

  gene_id <- clean_text(gene_id)

  good <- (
    !is.na(symbol) &
      !grepl(
        "^Ccac_",
        symbol
      )
  )

  out <- gene_id

  out[good] <- symbol[good]

  out
}


for (
  f in c(
    FAI_FILE,
    OLD_XTX,
    OLD_C2,
    NEW_WINDOWS_XTX,
    NEW_WINDOWS_C2
  )
) {

  check_file(f)

}


############################################################
## 3. READ TRUE SCAFFOLD LENGTHS
############################################################

fai <- fread(
  FAI_FILE,
  header = FALSE,
  showProgress = FALSE
)


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
  Scaffold_length_Mb :=
    Scaffold_length_bp /
    1e6
]


fai[
  ,
  Eligible_gt1Mb :=
    Scaffold_length_bp >
    MIN_LENGTH
]


############################################################
## 4. READ EXISTING ANNOTATED TABLES
############################################################

xtx_old <- fread(
  OLD_XTX,
  showProgress = FALSE
)


c2_old <- fread(
  OLD_C2,
  showProgress = FALSE
)


cat(
  "\nOld XtX regions: ",
  nrow(xtx_old),
  "\n",
  sep = ""
)


cat(
  "Old C2 regions: ",
  nrow(c2_old),
  "\n",
  sep = ""
)


############################################################
## 5. ADD NORMALIZED SCAFFOLD + LENGTH
############################################################

prepare_table <- function(dt) {

  dt <- copy(dt)


  dt[
    ,
    Scaffold_norm :=
      normalize_scaffold(
        Scaffold_accession
      )
  ]


  dt <- merge(

    dt,

    fai[
      ,
      .(
        Scaffold_norm,
        Scaffold_length_bp,
        Scaffold_length_Mb,
        Eligible_gt1Mb
      )
    ],

    by =
      "Scaffold_norm",

    all.x =
      TRUE,

    sort =
      FALSE

  )


  if (
    any(
      is.na(
        dt$Scaffold_length_bp
      )
    )
  ) {

    stop(
      "Some candidate scaffolds failed to map to FAI."
    )

  }


  dt

}


xtx_all <- prepare_table(
  xtx_old
)


c2_all <- prepare_table(
  c2_old
)


############################################################
## 6. FILTER TO >1 Mb
############################################################

xtx <- xtx_all[
  Eligible_gt1Mb == TRUE
]


c2 <- c2_all[
  Eligible_gt1Mb == TRUE
]


if (
  nrow(xtx) != 12L
) {

  stop(
    "Expected 12 retained XtX regions; found ",
    nrow(xtx)
  )

}


if (
  nrow(c2) != 17L
) {

  stop(
    "Expected 17 retained C2 regions; found ",
    nrow(c2)
  )

}


############################################################
## 7. VERIFY AGAINST FORMAL >1 Mb RERUN
############################################################

new_xtx <- fread(
  NEW_WINDOWS_XTX,
  showProgress = FALSE
)


new_c2 <- fread(
  NEW_WINDOWS_C2,
  showProgress = FALSE
)


new_xtx[
  ,
  Scaffold_norm :=
    normalize_scaffold(chr)
]


new_c2[
  ,
  Scaffold_norm :=
    normalize_scaffold(chr)
]


############################################################
## Match retained annotated regions to actual new windows
############################################################

xtx_check <- merge(

  xtx[
    ,
    .(
      Region_ID,
      Scaffold_norm,
      Window_start,
      Window_end
    )
  ],

  new_xtx[
    ,
    .(
      Scaffold_norm,
      Window_start_new = beg,
      Window_end_new = end
    )
  ],

  by = "Scaffold_norm",

  allow.cartesian = TRUE

)


xtx_check[
  ,
  Exact_match :=
    Window_start == Window_start_new &
    Window_end == Window_end_new
]


xtx_match_counts <- xtx_check[
  Exact_match == TRUE,
  .N,
  by = Region_ID
]


if (
  nrow(xtx_match_counts) != 12L ||
  any(
    xtx_match_counts$N != 1L
  )
) {

  stop(
    "XtX retained-table regions do not perfectly match ",
    "the formal >1 Mb rerun."
  )

}


c2_check <- merge(

  c2[
    ,
    .(
      Region_ID,
      Scaffold_norm,
      Window_start,
      Window_end
    )
  ],

  new_c2[
    ,
    .(
      Scaffold_norm,
      Window_start_new = beg,
      Window_end_new = end
    )
  ],

  by = "Scaffold_norm",

  allow.cartesian = TRUE

)


c2_check[
  ,
  Exact_match :=
    Window_start == Window_start_new &
    Window_end == Window_end_new
]


c2_match_counts <- c2_check[
  Exact_match == TRUE,
  .N,
  by = Region_ID
]


if (
  nrow(c2_match_counts) != 17L ||
  any(
    c2_match_counts$N != 1L
  )
) {

  stop(
    "C2 retained-table regions do not perfectly match ",
    "the formal >1 Mb rerun."
  )

}


############################################################
## 8. REMOVE INTERNAL QC COLUMNS FROM PUBLICATION TABLE
############################################################

publication_drop <- c(
  "Scaffold_norm",
  "Scaffold_length_bp",
  "Scaffold_length_Mb",
  "Eligible_gt1Mb"
)


xtx_public <- copy(xtx)

c2_public <- copy(c2)


xtx_public[
  ,
  (publication_drop) := NULL
]


c2_public[
  ,
  (publication_drop) := NULL
]


############################################################
## 9. WRITE UPDATED FINAL TABLES
############################################################

XTX_OUT <- file.path(
  OUTDIR,
  "Table_XtX_candidate_regions_allSNP_gt1Mb_FINAL.tsv"
)


C2_OUT <- file.path(
  OUTDIR,
  "Table_C2_Host_candidate_regions_allSNP_gt1Mb_FINAL.tsv"
)


fwrite(
  xtx_public,
  XTX_OUT,
  sep = "\t",
  quote = FALSE,
  na = ""
)


fwrite(
  c2_public,
  C2_OUT,
  sep = "\t",
  quote = FALSE,
  na = ""
)


############################################################
## 10. BASIC SUMMARY STATISTICS
############################################################

summarize_regions <- function(
  dt,
  analysis
) {

  out <- data.table(

    Analysis =
      analysis,

    N_regions =
      nrow(dt),

    N_scaffolds =
      uniqueN(
        dt$Scaffold_norm
      ),

    Cumulative_region_span_bp =
      sum(
        dt$Window_size_bp,
        na.rm = TRUE
      ),

    Cumulative_region_span_kb =
      sum(
        dt$Window_size_bp,
        na.rm = TRUE
      ) / 1000,

    N_region_SNP_memberships =
      sum(
        dt$N_SNPs,
        na.rm = TRUE
      )

  )


  ##########################################################
  ## POD columns if available
  ##########################################################

  pod_col <- grep(
    "N.*POD99\\.9.*SNP",
    names(dt),
    value = TRUE,
    ignore.case = TRUE
  )


  stat_col <- grep(
    "Statistic.*peak.*POD99\\.9",
    names(dt),
    value = TRUE,
    ignore.case = TRUE
  )


  ls_col <- grep(
    "Local.*Score.*peak.*POD99\\.9",
    names(dt),
    value = TRUE,
    ignore.case = TRUE
  )


  if (
    length(pod_col) == 1L
  ) {

    out[
      ,
      POD99.9_SNP_memberships :=
        sum(
          as.numeric(
            dt[[pod_col]]
          ),
          na.rm = TRUE
        )
    ]

  }


  if (
    length(stat_col) == 1L
  ) {

    x <- tolower(
      as.character(
        dt[[stat_col]]
      )
    )

    out[
      ,
      N_statistic_peaks_gt_POD99.9 :=
        sum(
          x %in% c(
            "yes",
            "y",
            "true",
            "1"
          )
        )
    ]

  }


  if (
    length(ls_col) == 1L
  ) {

    x <- tolower(
      as.character(
        dt[[ls_col]]
      )
    )

    out[
      ,
      N_local_score_peaks_gt_POD99.9 :=
        sum(
          x %in% c(
            "yes",
            "y",
            "true",
            "1"
          )
        )
    ]

  }


  out

}


summary_xtx <- summarize_regions(
  xtx,
  "XtX"
)


summary_c2 <- summarize_regions(
  c2,
  "C2"
)


summary_regions <- rbindlist(
  list(
    summary_xtx,
    summary_c2
  ),
  fill = TRUE
)


############################################################
## 11. GENE DISPLAY NAMES
############################################################

xtx[
  ,
  Candidate_gene :=
    informative_gene(
      Gene_symbol,
      Gene_ID
    )
]


c2[
  ,
  Candidate_gene :=
    informative_gene(
      Gene_symbol,
      Gene_ID
    )
]


############################################################
## 12. XtX-C2 PHYSICAL OVERLAPS
############################################################

overlap_list <- list()

k <- 1L


for (
  i in seq_len(
    nrow(xtx)
  )
) {

  for (
    j in seq_len(
      nrow(c2)
    )
  ) {

    if (
      xtx$Scaffold_norm[i] !=
      c2$Scaffold_norm[j]
    ) {

      next

    }


    overlap_start <- max(
      xtx$Window_start[i],
      c2$Window_start[j]
    )


    overlap_end <- min(
      xtx$Window_end[i],
      c2$Window_end[j]
    )


    if (
      overlap_start <= overlap_end
    ) {

      overlap_list[[k]] <- data.table(

        XtX_Region_ID =
          xtx$Region_ID[i],

        C2_Region_ID =
          c2$Region_ID[j],

        Scaffold =
          xtx$Scaffold_norm[i],

        XtX_start =
          xtx$Window_start[i],

        XtX_end =
          xtx$Window_end[i],

        C2_start =
          c2$Window_start[j],

        C2_end =
          c2$Window_end[j],

        Overlap_start =
          overlap_start,

        Overlap_end =
          overlap_end,

        Overlap_size_bp =
          overlap_end -
          overlap_start +
          1L,

        XtX_gene =
          xtx$Candidate_gene[i],

        C2_gene =
          c2$Candidate_gene[j]

      )


      k <- k + 1L

    }

  }

}


if (
  length(overlap_list) > 0L
) {

  physical_overlap <- rbindlist(
    overlap_list
  )

} else {

  physical_overlap <- data.table()
}


############################################################
## 13. PEAK-ASSOCIATED GENE OVERLAP
############################################################

xtx_gene_set <- unique(
  na.omit(
    xtx$Gene_ID
  )
)


c2_gene_set <- unique(
  na.omit(
    c2$Gene_ID
  )
)


shared_gene_ids <- intersect(
  xtx_gene_set,
  c2_gene_set
)


shared_genes <- unique(

  rbindlist(

    list(

      xtx[
        Gene_ID %in%
          shared_gene_ids,
        .(
          Gene_ID,
          Gene_symbol,
          Candidate_gene,
          XtX_Region_ID =
            Region_ID
        )
      ],

      c2[
        Gene_ID %in%
          shared_gene_ids,
        .(
          Gene_ID,
          Gene_symbol,
          Candidate_gene,
          C2_Region_ID =
            Region_ID
        )
      ]

    ),

    fill = TRUE

  )

)


############################################################
## Better one-row-per-shared-gene table
############################################################

shared_gene_summary <- data.table(
  Gene_ID = shared_gene_ids
)


if (
  nrow(shared_gene_summary) > 0L
) {

  shared_gene_summary <- merge(

    shared_gene_summary,

    xtx[
      Gene_ID %in%
        shared_gene_ids,
      .(
        Gene_ID,
        Gene_symbol_XtX =
          Gene_symbol,
        Gene_name_XtX =
          Candidate_gene,
        XtX_Region_ID =
          Region_ID
      )
    ],

    by = "Gene_ID",

    all.x = TRUE

  )


  shared_gene_summary <- merge(

    shared_gene_summary,

    c2[
      Gene_ID %in%
        shared_gene_ids,
      .(
        Gene_ID,
        Gene_symbol_C2 =
          Gene_symbol,
        Gene_name_C2 =
          Candidate_gene,
        C2_Region_ID =
          Region_ID
      )
    ],

    by = "Gene_ID",

    all.x = TRUE

  )

}


############################################################
## 14. REMOVED REGIONS FOR TRACEABILITY
############################################################

removed_xtx <- xtx_all[
  Eligible_gt1Mb == FALSE
]


removed_c2 <- c2_all[
  Eligible_gt1Mb == FALSE
]


removed_regions <- rbindlist(

  list(

    removed_xtx[
      ,
      .(
        Analysis =
          "XtX",

        Region_ID,

        Scaffold =
          Scaffold_norm,

        Scaffold_length_bp,

        Scaffold_length_Mb,

        Window_start,

        Window_end,

        Gene_ID,

        Gene_symbol,

        Candidate_gene =
          informative_gene(
            Gene_symbol,
            Gene_ID
          )
      )
    ],


    removed_c2[
      ,
      .(
        Analysis =
          "C2",

        Region_ID,

        Scaffold =
          Scaffold_norm,

        Scaffold_length_bp,

        Scaffold_length_Mb,

        Window_start,

        Window_end,

        Gene_ID,

        Gene_symbol,

        Candidate_gene =
          informative_gene(
            Gene_symbol,
            Gene_ID
          )
      )
    ]

  ),

  fill = TRUE

)


############################################################
## 15. WRITE SUMMARY TABLES
############################################################

fwrite(

  summary_regions,

  file.path(
    OUTDIR,
    "Table_FINAL_gt1Mb_candidate_summary.tsv"
  ),

  sep = "\t",
  quote = FALSE,
  na = ""

)


fwrite(

  physical_overlap,

  file.path(
    OUTDIR,
    "Table_FINAL_gt1Mb_XtX_vs_C2_physical_overlaps.tsv"
  ),

  sep = "\t",
  quote = FALSE,
  na = ""

)


fwrite(

  shared_gene_summary,

  file.path(
    OUTDIR,
    "Table_FINAL_gt1Mb_XtX_vs_C2_shared_peak_genes.tsv"
  ),

  sep = "\t",
  quote = FALSE,
  na = ""

)


fwrite(

  removed_regions,

  file.path(
    OUTDIR,
    "Table_REMOVED_short_scaffold_candidate_regions.tsv"
  ),

  sep = "\t",
  quote = FALSE,
  na = ""

)


############################################################
## 16. CONSOLE SUMMARY
############################################################

cat(
  "\n============================================\n"
)

cat(
  "FINAL >1 Mb CANDIDATE SUMMARY\n"
)

cat(
  "============================================\n\n"
)


print(
  summary_regions
)


cat(
  "\n============================================\n"
)

cat(
  "XtX-C2 PHYSICAL OVERLAPS\n"
)

cat(
  "============================================\n\n"
)


if (
  nrow(physical_overlap) > 0L
) {

  print(
    physical_overlap
  )

} else {

  cat(
    "No physical XtX-C2 overlaps.\n"
  )

}


cat(
  "\n============================================\n"
)

cat(
  "SHARED PEAK-ASSOCIATED GENES\n"
)

cat(
  "============================================\n\n"
)


if (
  nrow(shared_gene_summary) > 0L
) {

  print(
    shared_gene_summary
  )

} else {

  cat(
    "No shared peak-associated genes.\n"
  )

}


cat(
  "\n============================================\n"
)

cat(
  "REMOVED SHORT-SCAFFOLD REGIONS\n"
)

cat(
  "============================================\n\n"
)


print(
  removed_regions
)


cat(
  "\n============================================\n"
)

cat(
  "OUTPUTS\n"
)

cat(
  "============================================\n"
)


cat(
  "\n",
  XTX_OUT,
  "\n",
  C2_OUT,
  "\n",
  sep = ""
)


cat(
  "\nDONE\n"
)
