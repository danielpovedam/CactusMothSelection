#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

############################################################
## 28B_FINAL_gt1Mb_VEP_candidate_SNP_classification.R
##
## FINAL >1 Mb VEP consequence summary
##
## NO VEP rerun required.
##
## The revised XtX/C2 candidate sets are subsets of the old
## candidate set. We therefore:
##
## 1. read the new final candidate master;
## 2. read the already finalized VEP classifications;
## 3. retain only positions belonging to the new master;
## 4. replace old XtX/C2/GEA membership flags;
## 5. recompute consequence summaries.
############################################################


# ============================================================
# PATHS
# ============================================================

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"


NEW_DIR <- file.path(
  BASE,
  "selection_scan",
  "VEP_FINAL_gt1Mb_selection_candidates"
)


NEW_MASTER <- file.path(
  NEW_DIR,
  "ALL_selection_candidate_SNPs_master_gt1Mb_FINAL.tsv"
)


OLD_FINAL <- file.path(
  BASE,
  "selection_scan",
  "VEP_ALL_selection_candidates",
  "Table_ALL_candidate_SNPs_VEP_primary_consequence_FINAL.tsv"
)


OUT_FINAL <- file.path(
  NEW_DIR,
  "Table_ALL_candidate_SNPs_VEP_primary_consequence_gt1Mb_FINAL.tsv"
)


OUT_SUMMARY <- file.path(
  NEW_DIR,
  "Table_VEP_primary_consequence_summary_gt1Mb_FINAL.tsv"
)


OUT_ANALYSIS <- file.path(
  NEW_DIR,
  "Table_VEP_primary_consequence_by_analysis_gt1Mb_FINAL.tsv"
)


OUT_COMBINATION <- file.path(
  NEW_DIR,
  "Table_VEP_primary_consequence_by_analysis_combination_gt1Mb_FINAL.tsv"
)


OUT_PLOT_DATA <- file.path(
  NEW_DIR,
  "Figure_VEP_candidate_SNP_consequences_data_gt1Mb_FINAL.tsv"
)


OUT_PDF <- file.path(
  NEW_DIR,
  "Figure_VEP_candidate_SNP_consequences_gt1Mb_FINAL.pdf"
)


OUT_PNG <- file.path(
  NEW_DIR,
  "Figure_VEP_candidate_SNP_consequences_gt1Mb_FINAL.png"
)


# ============================================================
# CHECKS
# ============================================================

for (
  f in c(
    NEW_MASTER,
    OLD_FINAL
  )
) {

  if (
    !file.exists(f) ||
    file.info(f)$size <= 0
  ) {

    stop(
      "Missing or empty file:\n",
      f
    )
  }
}


# ============================================================
# CONSEQUENCE HIERARCHY
# ============================================================

class_levels <- c(
  "Protein-altering variants",
  "UTR variants",
  "Synonymous variants",
  "Non-coding exon/transcript variants",
  "Intronic or splice-region variants",
  "Upstream/downstream variants",
  "Intergenic variants"
)


# ============================================================
# READ NEW MASTER
# ============================================================

master <- fread(
  NEW_MASTER,
  showProgress = FALSE
)


master[
  ,
  POS :=
    as.integer(POS)
]


N_FINAL <- nrow(
  master
)


if (
  uniqueN(
    master,
    by = c(
      "CHROM",
      "POS"
    )
  ) != N_FINAL
) {

  stop(
    "New master contains duplicated CHROM/POS positions."
  )
}


cat(
  "\n=============================================\n"
)

cat(
  "FINAL >1 Mb VEP CLASSIFICATION\n"
)

cat(
  "=============================================\n"
)


cat(
  "New candidate positions: ",
  N_FINAL,
  "\n",
  sep = ""
)


# ============================================================
# READ OLD FINAL VEP CLASSIFICATION
# ============================================================

old <- fread(
  OLD_FINAL,
  showProgress = FALSE
)


old[
  ,
  POS :=
    as.integer(POS)
]


if (
  uniqueN(
    old,
    by = c(
      "CHROM",
      "POS"
    )
  ) != nrow(old)
) {

  stop(
    "Old VEP table contains duplicated positions."
  )
}


############################################################
## Confirm every new position was already VEP classified
############################################################

missing <- fsetdiff(

  master[
    ,
    .(
      CHROM,
      POS
    )
  ],

  old[
    ,
    .(
      CHROM,
      POS
    )
  ]

)


if (
  nrow(missing) > 0L
) {

  print(missing)

  stop(
    "Some FINAL candidate SNPs were not present in the old ",
    "VEP dataset. A VEP rerun would be required."
  )
}


cat(
  "All new candidate positions already have VEP classifications.\n"
)


# ============================================================
# REMOVE OLD MEMBERSHIP COLUMNS
# ============================================================

membership_cols <- c(
  "In_XtX",
  "In_C2",
  "In_GEA",
  "N_analyses",
  "Analysis_combination",
  "GEA_covariates",
  "GEA_max_BFdb"
)


old_annotation_cols <- setdiff(
  names(old),
  membership_cols
)


old_annotation <- old[
  ,
  ..old_annotation_cols
]


# ============================================================
# MERGE NEW MEMBERSHIP WITH EXISTING VEP ANNOTATIONS
# ============================================================

final <- merge(

  master,

  old_annotation,

  by = c(
    "CHROM",
    "POS"
  ),

  all.x = TRUE,

  sort = FALSE

)


if (
  nrow(final) != N_FINAL
) {

  stop(
    "Final merge produced ",
    nrow(final),
    " rows; expected ",
    N_FINAL
  )
}


if (
  anyNA(
    final$Primary_consequence_class
  )
) {

  bad <- final[
    is.na(
      Primary_consequence_class
    ),
    .(
      CHROM,
      POS
    )
  ]

  print(bad)

  stop(
    "Some final candidate positions lack VEP consequence classification."
  )
}


setorder(
  final,
  CHROM,
  POS
)


# ============================================================
# WRITE FINAL POSITION TABLE
# ============================================================

fwrite(
  final,
  OUT_FINAL,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# GLOBAL SUMMARY
# ============================================================

summary_template <- data.table(
  Primary_consequence_class =
    class_levels
)


summary_observed <- final[
  ,
  .(
    N = .N
  ),
  by =
    Primary_consequence_class
]


summary <- merge(

  summary_template,

  summary_observed,

  by =
    "Primary_consequence_class",

  all.x =
    TRUE,

  sort =
    FALSE

)


summary[
  is.na(N),
  N :=
    0L
]


summary[
  ,
  Percent :=
    100 *
    N /
    N_FINAL
]


summary[
  ,
  Order :=
    match(
      Primary_consequence_class,
      class_levels
    )
]


setorder(
  summary,
  Order
)


summary[
  ,
  Order :=
    NULL
]


if (
  sum(
    summary$N
  ) != N_FINAL
) {

  stop(
    "Global consequence counts do not sum to final candidate total."
  )
}


fwrite(
  summary,
  OUT_SUMMARY,
  sep = "\t",
  quote = FALSE
)


# ============================================================
# SUMMARY BY ANALYSIS
# ============================================================

analysis_definitions <- list(
  XtX = "In_XtX",
  C2 = "In_C2",
  GEA = "In_GEA"
)


analysis_tables <- list()


for (
  analysis_name in names(
    analysis_definitions
  )
) {

  flag <- analysis_definitions[[
    analysis_name
  ]]


  dt <- final[
    get(flag) == TRUE
  ]


  total_analysis <- nrow(
    dt
  )


  observed <- dt[
    ,
    .(
      N = .N
    ),
    by =
      Primary_consequence_class
  ]


  tmp <- merge(

    summary_template,

    observed,

    by =
      "Primary_consequence_class",

    all.x =
      TRUE,

    sort =
      FALSE

  )


  tmp[
    is.na(N),
    N :=
      0L
  ]


  tmp[
    ,
    Percent :=
      if (
        total_analysis > 0
      ) {

        100 *
        N /
        total_analysis

      } else {

        0
      }
  ]


  tmp[
    ,
    `:=`(
      Analysis =
        analysis_name,

      Total_candidate_positions =
        total_analysis,

      Order =
        match(
          Primary_consequence_class,
          class_levels
        )
    )
  ]


  setorder(
    tmp,
    Order
  )


  tmp[
    ,
    Order :=
      NULL
  ]


  analysis_tables[[
    analysis_name
  ]] <- tmp
}


analysis_summary <- rbindlist(
  analysis_tables,
  use.names = TRUE,
  fill = TRUE
)


setcolorder(
  analysis_summary,
  c(
    "Analysis",
    "Total_candidate_positions",
    "Primary_consequence_class",
    "N",
    "Percent"
  )
)


fwrite(
  analysis_summary,
  OUT_ANALYSIS,
  sep = "\t",
  quote = FALSE
)


# ============================================================
# SUMMARY BY ANALYSIS COMBINATION
# ============================================================

combination_summary <- final[
  ,
  .(
    N = .N
  ),
  by = .(
    Analysis_combination,
    Primary_consequence_class
  )
]


combination_summary[
  ,
  Total_combination :=
    sum(N),
  by =
    Analysis_combination
]


combination_summary[
  ,
  Percent :=
    100 *
    N /
    Total_combination
]


combination_summary[
  ,
  Class_order :=
    match(
      Primary_consequence_class,
      class_levels
    )
]


setorder(
  combination_summary,
  Analysis_combination,
  Class_order
)


combination_summary[
  ,
  Class_order :=
    NULL
]


fwrite(
  combination_summary,
  OUT_COMBINATION,
  sep = "\t",
  quote = FALSE
)


# ============================================================
# PLOT DATA
# ============================================================

plot_dt <- copy(
  summary
)


plot_dt[
  ,
  Primary_consequence_class :=
    factor(
      Primary_consequence_class,
      levels =
        class_levels
    )
]


fwrite(
  plot_dt,
  OUT_PLOT_DATA,
  sep = "\t",
  quote = FALSE
)


plot_nonzero <- plot_dt[
  N > 0
]


plot_nonzero[
  ,
  Label :=
    paste0(
      N,
      "\n",
      sprintf(
        "%.1f%%",
        Percent
      )
    )
]


# ============================================================
# PIE CHART
# ============================================================

p <- ggplot(

  plot_nonzero,

  aes(
    x = "",
    y = N,
    fill =
      Primary_consequence_class
  )

) +

  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.35
  ) +

  coord_polar(
    theta = "y"
  ) +

  geom_text(

    aes(
      label = Label
    ),

    position =
      position_stack(
        vjust = 0.5
      ),

    size = 3.0

  ) +

  labs(
    fill = NULL
  ) +

  theme_void(
    base_size = 11
  ) +

  theme(

    legend.position =
      "right",

    legend.text =
      element_text(
        size = 9.5
      ),

    plot.margin =
      margin(
        10,
        10,
        10,
        10
      )

  )


ggsave(
  OUT_PDF,
  p,
  width = 7.4,
  height = 5.2,
  units = "in"
)


ggsave(
  OUT_PNG,
  p,
  width = 7.4,
  height = 5.2,
  units = "in",
  dpi = 600
)


# ============================================================
# FINAL REPORT
# ============================================================

cat(
  "\nGLOBAL FINAL CONSEQUENCE SUMMARY\n\n"
)


print(
  summary
)


cat(
  "\nBY ANALYSIS\n\n"
)


print(
  analysis_summary
)


cat(
  "\nAnalysis combinations:\n"
)


print(
  final[
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


cat(
  "\nFinal classified positions: ",
  N_FINAL,
  "\n",
  sep = ""
)


cat(
  "\nOutputs:\n",
  OUT_FINAL,
  "\n",
  OUT_SUMMARY,
  "\n",
  OUT_ANALYSIS,
  "\n",
  sep = ""
)


cat(
  "\nDONE\n"
)
