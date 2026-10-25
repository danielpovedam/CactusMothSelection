#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

# ============================================================
# SETTINGS
# ============================================================

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

# ------------------------------------------------------------
# GEA BF >20 / 1-kb clustering
# ------------------------------------------------------------

CLUSTER_DIR <- file.path(
  BASE,
  "selection_scan",
  "GEA_BF20_distance_clustering_diagnostic"
)

GEA_REGION_FILE <- file.path(
  CLUSTER_DIR,
  "GEA_BF20_candidate_regions_gap_1000bp.tsv"
)

GEA_SNP_FILE <- file.path(
  CLUSTER_DIR,
  "GEA_BF20_SNPs_clustered_gap_1000bp.tsv"
)

# ------------------------------------------------------------
# Existing final XtX / C2 tables
# ------------------------------------------------------------

FINAL_TABLE_DIR <- file.path(
  BASE,
  "selection_scan",
  "FINAL_manuscript_tables"
)

XTX_FILE <- file.path(
  FINAL_TABLE_DIR,
  "Table_XtX_candidate_regions_allSNP_FINAL.tsv"
)

C2_FILE <- file.path(
  FINAL_TABLE_DIR,
  "Table_C2_Host_candidate_regions_allSNP_FINAL.tsv"
)

# ------------------------------------------------------------
# Cactus-moth GFF
# ------------------------------------------------------------

GFF_FILE <- file.path(
  "/home/dpoveda/micetral/Population_snps",
  "Ccac_functional-annotation",
  "02_agat",
  "Ccac",
  "out_abril4",
  "Ccac.annotations.corrected.gff3"
)

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUTDIR <- file.path(
  BASE,
  "selection_scan",
  "FINAL_GEA_BF20_1kb"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Final analysis definitions
# ------------------------------------------------------------

BF_THRESHOLD <- 20
CLUSTER_GAP_BP <- 1000L
GENE_PROXIMITY_BP <- 5000L

PRIMARY_COVARIATES <- c(
  "TempPC1",
  "PrecipPC1",
  "AridityIndex",
  "Elevation",
  "Latitude",
  "Longitude"
)

SENSITIVITY_COVARIATES <- c(
  "BIO1",
  "BIO12"
)

EXPECTED_COVARIATES <- c(
  PRIMARY_COVARIATES,
  SENSITIVITY_COVARIATES
)

# ============================================================
# HELPERS
# ============================================================

clean_scaffold <- function(x) {

  x <- as.character(x)

  x <- sub("^lcl\\|", "", x)
  x <- sub("l_1_1$", "", x)
  x <- sub("l_1$", "", x)

  x
}


extract_attr <- function(x, key) {

  pattern <- paste0(
    "(^|;)",
    key,
    "=([^;]+)"
  )

  m <- regexec(
    pattern,
    x,
    perl = TRUE
  )

  z <- regmatches(
    x,
    m
  )

  out <- vapply(
    z,
    function(xx) {

      if (length(xx) >= 3L) {
        xx[3L]
      } else {
        NA_character_
      }

    },
    character(1)
  )

  out[
    out == ""
  ] <- NA_character_

  out
}


first_nonmissing <- function(...) {

  values <- list(...)

  n <- max(
    vapply(
      values,
      length,
      integer(1)
    )
  )

  out <- rep(
    NA_character_,
    n
  )

  for (v in values) {

    v <- as.character(v)

    if (
      length(v) == 1L &&
      n > 1L
    ) {

      v <- rep(
        v,
        n
      )
    }

    take <-
      is.na(out) &
      !is.na(v) &
      v != ""

    out[take] <- v[take]
  }

  out
}


collapse_unique <- function(x) {

  x <- as.character(x)

  x <- unique(
    x[
      !is.na(x) &
      x != ""
    ]
  )

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(
    sort(x),
    collapse = ";"
  )
}


safe_uniqueN <- function(x) {

  x <- x[
    !is.na(x)
  ]

  if (is.character(x)) {

    x <- x[
      x != ""
    ]
  }

  if (length(x) == 0L) {
    return(0L)
  }

  as.integer(
    uniqueN(x)
  )
}


find_column <- function(
  dt,
  candidates,
  description
) {

  hit <- candidates[
    candidates %in%
    names(dt)
  ]

  if (length(hit) == 0L) {

    stop(
      "Could not identify ",
      description,
      ". Tried: ",
      paste(
        candidates,
        collapse = ", "
      )
    )
  }

  hit[1L]
}


classify_peak_context <- function(
  annotation_available,
  in_cds,
  in_exon,
  in_gene,
  near_gene
) {

  if (!isTRUE(annotation_available)) {
    return("annotation_unavailable")
  }

  if (isTRUE(in_cds)) {
    return("coding")
  }

  if (isTRUE(in_exon)) {
    return("exonic_non_coding")
  }

  if (isTRUE(in_gene)) {
    return("intronic")
  }

  if (isTRUE(near_gene)) {
    return("gene_proximal_5kb")
  }

  "distal_intergenic"
}


standardize_region_table <- function(
  dt,
  test_name
) {

  chr_col <- find_column(
    dt,
    c(
      "Scaffold",
      "chr",
      "CHROM"
    ),
    paste0(
      test_name,
      " scaffold column"
    )
  )

  start_col <- find_column(
    dt,
    c(
      "Window_start",
      "beg",
      "start",
      "Region_start"
    ),
    paste0(
      test_name,
      " start column"
    )
  )

  end_col <- find_column(
    dt,
    c(
      "Window_end",
      "end",
      "Region_end"
    ),
    paste0(
      test_name,
      " end column"
    )
  )

  out <- data.table(

    Region_index =
      seq_len(
        nrow(dt)
      ),

    Scaffold =
      clean_scaffold(
        dt[[chr_col]]
      ),

    Region_start =
      as.integer(
        dt[[start_col]]
      ),

    Region_end =
      as.integer(
        dt[[end_col]]
      )
  )

  out <- out[
    !is.na(Scaffold) &
    !is.na(Region_start) &
    !is.na(Region_end) &
    Region_end >= Region_start
  ]

  out[
    ,
    Test :=
      test_name
  ]

  out
}


empty_overlap_table <- function(prefix) {

  out <- data.table(

    Analysis_class =
      character(),

    Covariate =
      character(),

    GEA_locus_ID =
      character(),

    Scaffold =
      character(),

    GEA_start =
      integer(),

    GEA_end =
      integer(),

    GEA_peak_position =
      integer(),

    GEA_peak_BFdb =
      numeric(),

    GEA_Gene_ID =
      character(),

    GEA_Gene_symbol =
      character(),

    GEA_Peak_context =
      character(),

    Other_region_index =
      integer(),

    Other_region_start =
      integer(),

    Other_region_end =
      integer()
  )

  setnames(
    out,
    c(
      "Other_region_index",
      "Other_region_start",
      "Other_region_end"
    ),
    c(
      paste0(
        prefix,
        "_region_index"
      ),
      paste0(
        prefix,
        "_region_start"
      ),
      paste0(
        prefix,
        "_region_end"
      )
    )
  )

  out
}


find_region_overlaps <- function(
  gea_dt,
  other_dt,
  prefix
) {

  template <- empty_overlap_table(
    prefix
  )

  if (
    nrow(gea_dt) == 0L ||
    nrow(other_dt) == 0L
  ) {

    return(
      template
    )
  }

  shared_scaffolds <- intersect(
    unique(
      gea_dt$Scaffold
    ),
    unique(
      other_dt$Scaffold
    )
  )

  if (length(shared_scaffolds) == 0L) {

    return(
      template
    )
  }

  gea_sub <- gea_dt[
    Scaffold %in%
    shared_scaffolds
  ]

  other_sub <- other_dt[
    Scaffold %in%
    shared_scaffolds
  ]

  if (
    nrow(gea_sub) == 0L ||
    nrow(other_sub) == 0L
  ) {

    return(
      template
    )
  }

  gea_gr <- GRanges(
    seqnames =
      gea_sub$Scaffold,
    ranges =
      IRanges(
        start =
          gea_sub$Region_start,
        end =
          gea_sub$Region_end
      )
  )

  other_gr <- GRanges(
    seqnames =
      other_sub$Scaffold,
    ranges =
      IRanges(
        start =
          other_sub$Region_start,
        end =
          other_sub$Region_end
      )
  )

  hits <- suppressWarnings(
    findOverlaps(
      gea_gr,
      other_gr,
      ignore.strand = TRUE
    )
  )

  if (length(hits) == 0L) {

    return(
      template
    )
  }

  q <- queryHits(
    hits
  )

  s <- subjectHits(
    hits
  )

  out <- data.table(

    Analysis_class =
      gea_sub$Analysis_class[q],

    Covariate =
      gea_sub$Covariate[q],

    GEA_locus_ID =
      gea_sub$GEA_locus_ID[q],

    Scaffold =
      gea_sub$Scaffold[q],

    GEA_start =
      gea_sub$Region_start[q],

    GEA_end =
      gea_sub$Region_end[q],

    GEA_peak_position =
      gea_sub$Peak_position[q],

    GEA_peak_BFdb =
      gea_sub$Peak_BFdb[q],

    GEA_Gene_ID =
      gea_sub$Gene_ID[q],

    GEA_Gene_symbol =
      gea_sub$Gene_symbol[q],

    GEA_Peak_context =
      gea_sub$Peak_context[q],

    Other_region_index =
      other_sub$Region_index[s],

    Other_region_start =
      other_sub$Region_start[s],

    Other_region_end =
      other_sub$Region_end[s]
  )

  setnames(
    out,
    c(
      "Other_region_index",
      "Other_region_start",
      "Other_region_end"
    ),
    c(
      paste0(
        prefix,
        "_region_index"
      ),
      paste0(
        prefix,
        "_region_start"
      ),
      paste0(
        prefix,
        "_region_end"
      )
    )
  )

  unique(
    out
  )
}

# ============================================================
# INPUT CHECK
# ============================================================

required_files <- c(
  GEA_REGION_FILE,
  GEA_SNP_FILE,
  XTX_FILE,
  C2_FILE,
  GFF_FILE
)

cat("\n=============================================\n")
cat("INPUT CHECK\n")
cat("=============================================\n")

for (f in required_files) {

  if (!file.exists(f)) {

    stop(
      "Missing required file: ",
      f
    )
  }

  cat(
    "OK: ",
    f,
    "\n",
    sep = ""
  )
}

# ============================================================
# LOAD GEA LOCI
# ============================================================

cat("\n=============================================\n")
cat("LOAD FINAL BF >20 / 1-kb GEA LOCI\n")
cat("=============================================\n")

gea <- fread(
  GEA_REGION_FILE,
  showProgress = FALSE
)

required_gea <- c(
  "COVARIABLE",
  "Covariate",
  "Scaffold",
  "Cluster_ID",
  "Cluster_number",
  "Region_start",
  "Region_end",
  "Region_size_bp",
  "N_BF20_SNPs",
  "Peak_position",
  "Peak_MRK",
  "Peak_BFdb",
  "Peak_REF",
  "Peak_ALT"
)

missing_gea <- setdiff(
  required_gea,
  names(gea)
)

if (length(missing_gea) > 0L) {

  stop(
    "GEA table missing columns: ",
    paste(
      missing_gea,
      collapse = ", "
    )
  )
}

gea[
  ,
  Scaffold :=
    clean_scaffold(
      Scaffold
    )
]

integer_columns <- c(
  "COVARIABLE",
  "Cluster_number",
  "Region_start",
  "Region_end",
  "Region_size_bp",
  "N_BF20_SNPs",
  "Peak_position",
  "Peak_MRK"
)

for (cc in integer_columns) {

  set(
    gea,
    j = cc,
    value =
      as.integer(
        gea[[cc]]
      )
  )
}

gea[
  ,
  Peak_BFdb :=
    as.numeric(
      Peak_BFdb
    )
]

gea[
  ,
  Analysis_class :=
    fcase(

      Covariate %in%
      PRIMARY_COVARIATES,
      "Primary",

      Covariate %in%
      SENSITIVITY_COVARIATES,
      "Sensitivity",

      default =
        "Unexpected"
    )
]

if (
  any(
    gea$Analysis_class ==
    "Unexpected"
  )
) {

  bad_cov <- unique(
    gea[
      Analysis_class ==
      "Unexpected",
      Covariate
    ]
  )

  stop(
    "Unexpected GEA covariates: ",
    paste(
      bad_cov,
      collapse = ", "
    )
  )
}

gea[
  ,
  GEA_locus_ID :=
    paste0(
      "GEA_",
      Covariate,
      "_",
      Scaffold,
      "_",
      Region_start,
      "_",
      Region_end
    )
]

if (
  anyDuplicated(
    gea$GEA_locus_ID
  ) > 0L
) {

  stop(
    "GEA_locus_ID is not unique."
  )
}

cat(
  "Total GEA loci = ",
  nrow(gea),
  "\n",
  sep = ""
)

cat(
  "Primary GEA loci = ",
  nrow(
    gea[
      Analysis_class ==
      "Primary"
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "Sensitivity loci = ",
  nrow(
    gea[
      Analysis_class ==
      "Sensitivity"
    ]
  ),
  "\n",
  sep = ""
)

# ============================================================
# LOAD BF >20 SNP TABLE
# ============================================================

cat("\n=============================================\n")
cat("LOAD BF >20 SNP TABLE\n")
cat("=============================================\n")

gea_snps <- fread(
  GEA_SNP_FILE,
  showProgress = FALSE
)

required_snps <- c(
  "COVARIABLE",
  "Covariate",
  "Scaffold",
  "pos",
  "MRK",
  "REF",
  "ALT",
  "BFdb",
  "Cluster_ID"
)

missing_snps <- setdiff(
  required_snps,
  names(gea_snps)
)

if (length(missing_snps) > 0L) {

  stop(
    "GEA SNP table missing columns: ",
    paste(
      missing_snps,
      collapse = ", "
    )
  )
}

gea_snps[
  ,
  Scaffold :=
    clean_scaffold(
      Scaffold
    )
]

gea_snps[
  ,
  pos :=
    as.integer(
      pos
    )
]

gea_snps[
  ,
  MRK :=
    as.integer(
      MRK
    )
]

gea_snps[
  ,
  COVARIABLE :=
    as.integer(
      COVARIABLE
    )
]

gea_snps[
  ,
  BFdb :=
    as.numeric(
      BFdb
    )
]

gea_snps[
  ,
  Analysis_class :=
    fcase(

      Covariate %in%
      PRIMARY_COVARIATES,
      "Primary",

      Covariate %in%
      SENSITIVITY_COVARIATES,
      "Sensitivity",

      default =
        "Unexpected"
    )
]

if (
  any(
    gea_snps$Analysis_class ==
    "Unexpected"
  )
) {

  stop(
    "Unexpected covariate in GEA SNP table."
  )
}

if (
  any(
    gea_snps$BFdb <=
    BF_THRESHOLD
  )
) {

  stop(
    "BF >20 SNP table contains BF <= ",
    BF_THRESHOLD
  )
}

cat(
  "BF >20 SNP-covariate associations = ",
  nrow(gea_snps),
  "\n",
  sep = ""
)

# ============================================================
# LOAD GFF
# ============================================================

cat("\n=============================================\n")
cat("LOAD GFF\n")
cat("=============================================\n")

gff <- fread(
  GFF_FILE,
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  quote = "",
  fill = TRUE,
  showProgress = TRUE
)

if (ncol(gff) < 9L) {

  stop(
    "Expected at least 9 GFF columns."
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
  Scaffold :=
    clean_scaffold(
      seqid
    )
]

gff[
  ,
  start :=
    as.integer(
      start
    )
]

gff[
  ,
  end :=
    as.integer(
      end
    )
]

gff <- gff[
  !is.na(Scaffold) &
  !is.na(start) &
  !is.na(end) &
  end >= start
]

gff_scaffolds <- unique(
  gff$Scaffold
)

gea[
  ,
  Annotation_available :=
    Scaffold %in%
    gff_scaffolds
]

cat(
  "GEA loci with annotation available = ",
  sum(
    gea$Annotation_available ==
    TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "GEA loci with annotation unavailable = ",
  sum(
    gea$Annotation_available ==
    FALSE
  ),
  "\n",
  sep = ""
)

missing_annotation_scaffolds <- unique(
  gea[
    Annotation_available ==
    FALSE,
    Scaffold
  ]
)

if (
  length(
    missing_annotation_scaffolds
  ) > 0L
) {

  cat(
    "GEA scaffolds absent from GFF:\n"
  )

  print(
    missing_annotation_scaffolds
  )
}

# ============================================================
# BUILD GENE TABLE
# ============================================================

genes <- gff[
  feature ==
  "gene"
]

genes[
  ,
  attr_ID :=
    extract_attr(
      attributes,
      "ID"
    )
]

genes[
  ,
  attr_Name :=
    extract_attr(
      attributes,
      "Name"
    )
]

genes[
  ,
  attr_gene :=
    extract_attr(
      attributes,
      "gene"
    )
]

genes[
  ,
  attr_gene_id :=
    extract_attr(
      attributes,
      "gene_id"
    )
]

genes[
  ,
  attr_description :=
    extract_attr(
      attributes,
      "description"
    )
]

genes[
  ,
  attr_product :=
    extract_attr(
      attributes,
      "product"
    )
]

genes[
  ,
  Gene_ID :=
    first_nonmissing(
      attr_gene_id,
      attr_ID
    )
]

genes[
  ,
  Gene_symbol :=
    first_nonmissing(
      attr_gene,
      attr_Name
    )
]

genes[
  ,
  Description :=
    first_nonmissing(
      attr_description,
      attr_product
    )
]

genes[
  Gene_symbol ==
  Gene_ID,
  Gene_symbol :=
    NA_character_
]

cat(
  "Genes loaded = ",
  nrow(genes),
  "\n",
  sep = ""
)

# ============================================================
# INITIALIZE ANNOTATION COLUMNS
# ============================================================

gea[
  ,
  Region_overlaps_gene :=
    FALSE
]

gea[
  ,
  Region_overlapping_gene_IDs :=
    NA_character_
]

gea[
  ,
  Region_overlapping_gene_symbols :=
    NA_character_
]

gea[
  ,
  Peak_in_gene :=
    FALSE
]

gea[
  ,
  Peak_in_exon :=
    FALSE
]

gea[
  ,
  Peak_in_CDS :=
    FALSE
]

gea[
  ,
  Gene_ID :=
    NA_character_
]

gea[
  ,
  Gene_symbol :=
    NA_character_
]

gea[
  ,
  Description :=
    NA_character_
]

gea[
  ,
  Distance_to_gene_bp :=
    NA_integer_
]

# ============================================================
# ANNOTATION
# ============================================================

cat("\n=============================================\n")
cat("ANNOTATE GEA LOCI\n")
cat("=============================================\n")

annot_idx <- which(
  gea$Annotation_available ==
  TRUE
)

if (length(annot_idx) > 0L) {

  gea_annot <- gea[
    annot_idx
  ]

  relevant_scaffolds <- unique(
    gea_annot$Scaffold
  )

  genes_sub <- genes[
    Scaffold %in%
    relevant_scaffolds
  ]

  exon_sub <- gff[
    feature ==
    "exon" &
    Scaffold %in%
    relevant_scaffolds
  ]

  cds_sub <- gff[
    feature ==
    "CDS" &
    Scaffold %in%
    relevant_scaffolds
  ]

  region_gr <- GRanges(
    seqnames =
      gea_annot$Scaffold,
    ranges =
      IRanges(
        start =
          gea_annot$Region_start,
        end =
          gea_annot$Region_end
      )
  )

  peak_gr <- GRanges(
    seqnames =
      gea_annot$Scaffold,
    ranges =
      IRanges(
        start =
          gea_annot$Peak_position,
        end =
          gea_annot$Peak_position
      )
  )

  gene_gr <- GRanges(
    seqnames =
      genes_sub$Scaffold,
    ranges =
      IRanges(
        start =
          genes_sub$start,
        end =
          genes_sub$end
      )
  )

  mcols(gene_gr)$Gene_ID <-
    genes_sub$Gene_ID

  mcols(gene_gr)$Gene_symbol <-
    genes_sub$Gene_symbol

  mcols(gene_gr)$Description <-
    genes_sub$Description

  exon_gr <- GRanges(
    seqnames =
      exon_sub$Scaffold,
    ranges =
      IRanges(
        start =
          exon_sub$start,
        end =
          exon_sub$end
      )
  )

  cds_gr <- GRanges(
    seqnames =
      cds_sub$Scaffold,
    ranges =
      IRanges(
        start =
          cds_sub$start,
        end =
          cds_sub$end
      )
  )

  # ----------------------------------------------------------
  # REGION -> GENE OVERLAPS
  # ----------------------------------------------------------

  if (length(gene_gr) > 0L) {

    region_gene_hits <- suppressWarnings(
      findOverlaps(
        region_gr,
        gene_gr,
        ignore.strand = TRUE
      )
    )

    region_bool <- suppressWarnings(
      overlapsAny(
        region_gr,
        gene_gr,
        ignore.strand = TRUE
      )
    )

    gea[
      annot_idx,
      Region_overlaps_gene :=
        region_bool
    ]

    if (length(region_gene_hits) > 0L) {

      region_hit_dt <- data.table(

        local_idx =
          queryHits(
            region_gene_hits
          ),

        Gene_ID =
          mcols(gene_gr)$Gene_ID[
            subjectHits(
              region_gene_hits
            )
          ],

        Gene_symbol =
          mcols(gene_gr)$Gene_symbol[
            subjectHits(
              region_gene_hits
            )
          ]
      )

      region_hit_sum <- region_hit_dt[
        ,
        .(

          Region_overlapping_gene_IDs =
            collapse_unique(
              Gene_ID
            ),

          Region_overlapping_gene_symbols =
            collapse_unique(
              Gene_symbol
            )
        ),
        by =
          local_idx
      ]

      global_idx <- annot_idx[
        region_hit_sum$local_idx
      ]

      gea[
        global_idx,
        Region_overlapping_gene_IDs :=
          region_hit_sum$Region_overlapping_gene_IDs
      ]

      gea[
        global_idx,
        Region_overlapping_gene_symbols :=
          region_hit_sum$Region_overlapping_gene_symbols
      ]
    }
  }

  # ----------------------------------------------------------
  # PEAK INSIDE GENE
  # ----------------------------------------------------------

  if (length(gene_gr) > 0L) {

    peak_gene_bool <- suppressWarnings(
      overlapsAny(
        peak_gr,
        gene_gr,
        ignore.strand = TRUE
      )
    )

    gea[
      annot_idx,
      Peak_in_gene :=
        peak_gene_bool
    ]
  }

  # ----------------------------------------------------------
  # PEAK INSIDE EXON
  # ----------------------------------------------------------

  if (length(exon_gr) > 0L) {

    peak_exon_bool <- suppressWarnings(
      overlapsAny(
        peak_gr,
        exon_gr,
        ignore.strand = TRUE
      )
    )

    gea[
      annot_idx,
      Peak_in_exon :=
        peak_exon_bool
    ]
  }

  # ----------------------------------------------------------
  # PEAK INSIDE CDS
  # ----------------------------------------------------------

  if (length(cds_gr) > 0L) {

    peak_cds_bool <- suppressWarnings(
      overlapsAny(
        peak_gr,
        cds_gr,
        ignore.strand = TRUE
      )
    )

    gea[
      annot_idx,
      Peak_in_CDS :=
        peak_cds_bool
    ]
  }

  # ----------------------------------------------------------
  # NEAREST GENE <= 5 kb
  # ----------------------------------------------------------

  if (
    length(peak_gr) > 0L &&
    length(gene_gr) > 0L
  ) {

    near_hits <- suppressWarnings(
      distanceToNearest(
        peak_gr,
        gene_gr,
        ignore.strand = TRUE
      )
    )

    if (length(near_hits) > 0L) {

      q <- queryHits(
        near_hits
      )

      s <- subjectHits(
        near_hits
      )

      dd <- as.integer(
        mcols(
          near_hits
        )$distance
      )

      keep <-
        !is.na(dd) &
        dd <=
        GENE_PROXIMITY_BP

      q <- q[keep]
      s <- s[keep]
      dd <- dd[keep]

      if (length(q) > 0L) {

        global_idx <- annot_idx[
          q
        ]

        gea[
          global_idx,
          Gene_ID :=
            mcols(gene_gr)$Gene_ID[
              s
            ]
        ]

        gea[
          global_idx,
          Gene_symbol :=
            mcols(gene_gr)$Gene_symbol[
              s
            ]
        ]

        gea[
          global_idx,
          Description :=
            mcols(gene_gr)$Description[
              s
            ]
        ]

        gea[
          global_idx,
          Distance_to_gene_bp :=
            dd
        ]
      }
    }
  }
}

# ============================================================
# CLEAN GENE IDS
# ============================================================

gea[
  ,
  Gene_ID :=
    sub(
      "^gene-",
      "",
      Gene_ID
    )
]

# ============================================================
# FINAL PEAK CONTEXT
# ============================================================

gea[
  ,
  Peak_context :=
    mapply(
      classify_peak_context,

      Annotation_available,

      Peak_in_CDS,

      Peak_in_exon,

      Peak_in_gene,

      !is.na(
        Gene_ID
      ),

      USE.NAMES =
        FALSE
    )
]

# ============================================================
# QC CONTEXT
# ============================================================

if (
  any(
    gea[
      Annotation_available ==
      FALSE,
      Peak_context
    ] !=
    "annotation_unavailable"
  )
) {

  stop(
    "Annotation-unavailable scaffold was assigned another context."
  )
}

# ============================================================
# COVARIATE ORDER
# ============================================================

gea[
  ,
  Covariate_order :=
    match(
      Covariate,
      EXPECTED_COVARIATES
    )
]

setorder(
  gea,
  Covariate_order,
  Scaffold,
  Region_start
)

# ============================================================
# SUMMARY BY COVARIATE
# ============================================================

summary_by_covariate <- gea[
  ,
  .(

    N_BF20_SNPs =
      as.integer(
        sum(
          N_BF20_SNPs
        )
      ),

    N_candidate_loci =
      as.integer(
        .N
      ),

    N_single_SNP_loci =
      as.integer(
        sum(
          N_BF20_SNPs ==
          1
        )
      ),

    N_multi_SNP_loci =
      as.integer(
        sum(
          N_BF20_SNPs >
          1
        )
      ),

    N_scaffolds =
      as.integer(
        uniqueN(
          Scaffold
        )
      ),

    N_annotation_available =
      as.integer(
        sum(
          Annotation_available ==
          TRUE
        )
      ),

    N_annotation_unavailable =
      as.integer(
        sum(
          Annotation_available ==
          FALSE
        )
      ),

    N_peak_gene_associated_5kb =
      as.integer(
        sum(
          !is.na(
            Gene_ID
          )
        )
      ),

    N_unique_peak_genes =
      safe_uniqueN(
        Gene_ID
      ),

    Max_BFdb =
      as.numeric(
        max(
          Peak_BFdb
        )
      ),

    Median_peak_BFdb =
      as.numeric(
        median(
          Peak_BFdb
        )
      )
  ),
  by = .(
    Analysis_class,
    COVARIABLE,
    Covariate
  )
]

summary_by_covariate[
  ,
  Covariate_order :=
    match(
      Covariate,
      EXPECTED_COVARIATES
    )
]

setorder(
  summary_by_covariate,
  Covariate_order
)

# ============================================================
# SUMMARY PRIMARY VS SENSITIVITY
# ============================================================

summary_by_analysis <- gea[
  ,
  .(

    N_candidate_loci =
      as.integer(
        .N
      ),

    N_BF20_SNPs =
      as.integer(
        sum(
          N_BF20_SNPs
        )
      ),

    N_scaffolds =
      as.integer(
        uniqueN(
          Scaffold
        )
      ),

    N_single_SNP_loci =
      as.integer(
        sum(
          N_BF20_SNPs ==
          1
        )
      ),

    N_multi_SNP_loci =
      as.integer(
        sum(
          N_BF20_SNPs >
          1
        )
      ),

    N_annotation_available =
      as.integer(
        sum(
          Annotation_available ==
          TRUE
        )
      ),

    N_annotation_unavailable =
      as.integer(
        sum(
          Annotation_available ==
          FALSE
        )
      ),

    N_gene_associated_peaks_5kb =
      as.integer(
        sum(
          !is.na(
            Gene_ID
          )
        )
      ),

    N_unique_peak_genes =
      safe_uniqueN(
        Gene_ID
      )
  ),
  by =
    Analysis_class
]

# ============================================================
# CONTEXT SUMMARY
# ============================================================

context_summary <- gea[
  ,
  .(
    N_loci =
      as.integer(
        .N
      )
  ),
  by = .(
    Analysis_class,
    Covariate,
    Peak_context
  )
]

# ============================================================
# LOAD XtX / C2
# ============================================================

cat("\n=============================================\n")
cat("LOAD XtX AND C2 CANDIDATE REGIONS\n")
cat("=============================================\n")

xtx_raw <- fread(
  XTX_FILE,
  showProgress = FALSE
)

c2_raw <- fread(
  C2_FILE,
  showProgress = FALSE
)

xtx_regions <- standardize_region_table(
  xtx_raw,
  "XtX"
)

c2_regions <- standardize_region_table(
  c2_raw,
  "C2"
)

cat(
  "XtX regions = ",
  nrow(
    xtx_regions
  ),
  "\n",
  sep = ""
)

cat(
  "C2 regions = ",
  nrow(
    c2_regions
  ),
  "\n",
  sep = ""
)

# ============================================================
# PHYSICAL GEA LOCUS x XtX OVERLAP
# ============================================================

gea_xtx <- find_region_overlaps(
  gea,
  xtx_regions,
  "XtX"
)

gea_xtx_primary <- gea_xtx[
  Analysis_class ==
  "Primary"
]

gea_xtx_sensitivity <- gea_xtx[
  Analysis_class ==
  "Sensitivity"
]

# ============================================================
# PHYSICAL GEA LOCUS x C2 OVERLAP
# ============================================================

gea_c2 <- find_region_overlaps(
  gea,
  c2_regions,
  "C2"
)

gea_c2_primary <- gea_c2[
  Analysis_class ==
  "Primary"
]

gea_c2_sensitivity <- gea_c2[
  Analysis_class ==
  "Sensitivity"
]

# ============================================================
# ADD OVERLAP STATUS TO GEA
# ============================================================

gea[
  ,
  Overlaps_XtX :=
    GEA_locus_ID %in%
    gea_xtx$GEA_locus_ID
]

gea[
  ,
  Overlaps_C2 :=
    GEA_locus_ID %in%
    gea_c2$GEA_locus_ID
]

gea[
  ,
  Overlaps_XtX_and_C2 :=
    Overlaps_XtX &
    Overlaps_C2
]

# ============================================================
# COLLAPSE XtX REGION INFORMATION
# ============================================================

if (nrow(gea_xtx) > 0L) {

  xtx_collapsed <- gea_xtx[
    ,
    .(

      XtX_region_indices =
        collapse_unique(
          XtX_region_index
        ),

      XtX_regions =
        collapse_unique(
          paste0(
            Scaffold,
            ":",
            XtX_region_start,
            "-",
            XtX_region_end
          )
        )
    ),
    by =
      GEA_locus_ID
  ]

  gea <- merge(
    gea,
    xtx_collapsed,
    by =
      "GEA_locus_ID",
    all.x = TRUE,
    sort = FALSE
  )

} else {

  gea[
    ,
    XtX_region_indices :=
      NA_character_
  ]

  gea[
    ,
    XtX_regions :=
      NA_character_
  ]
}

# ============================================================
# COLLAPSE C2 REGION INFORMATION
# ============================================================

if (nrow(gea_c2) > 0L) {

  c2_collapsed <- gea_c2[
    ,
    .(

      C2_region_indices =
        collapse_unique(
          C2_region_index
        ),

      C2_regions =
        collapse_unique(
          paste0(
            Scaffold,
            ":",
            C2_region_start,
            "-",
            C2_region_end
          )
        )
    ),
    by =
      GEA_locus_ID
  ]

  gea <- merge(
    gea,
    c2_collapsed,
    by =
      "GEA_locus_ID",
    all.x = TRUE,
    sort = FALSE
  )

} else {

  gea[
    ,
    C2_region_indices :=
      NA_character_
  ]

  gea[
    ,
    C2_regions :=
      NA_character_
  ]
}

# ============================================================
# TRIPLE OVERLAP
# ============================================================

triple_loci <- gea[
  Overlaps_XtX_and_C2 ==
  TRUE
]

triple_primary <- triple_loci[
  Analysis_class ==
  "Primary"
]

triple_sensitivity <- triple_loci[
  Analysis_class ==
  "Sensitivity"
]

# ============================================================
# OVERLAP SUMMARY BY COVARIATE
# ============================================================

overlap_summary_by_covariate <- gea[
  ,
  .(

    N_GEA_loci =
      as.integer(
        .N
      ),

    N_GEA_loci_overlapping_XtX =
      as.integer(
        sum(
          Overlaps_XtX ==
          TRUE
        )
      ),

    Percent_GEA_loci_overlapping_XtX =
      as.numeric(
        100 *
        mean(
          Overlaps_XtX ==
          TRUE
        )
      ),

    N_GEA_loci_overlapping_C2 =
      as.integer(
        sum(
          Overlaps_C2 ==
          TRUE
        )
      ),

    Percent_GEA_loci_overlapping_C2 =
      as.numeric(
        100 *
        mean(
          Overlaps_C2 ==
          TRUE
        )
      ),

    N_GEA_loci_overlapping_XtX_and_C2 =
      as.integer(
        sum(
          Overlaps_XtX_and_C2 ==
          TRUE
        )
      )
  ),
  by = .(
    Analysis_class,
    Covariate
  )
]

# ============================================================
# GLOBAL OVERLAP SUMMARY
# ============================================================

primary_gea <- gea[
  Analysis_class ==
  "Primary"
]

sensitivity_gea <- gea[
  Analysis_class ==
  "Sensitivity"
]

overlap_global_summary <- data.table(

  Metric = c(

    "Primary_GEA_loci",

    "Primary_GEA_loci_overlapping_XtX",

    "Unique_XtX_regions_overlapping_primary_GEA",

    "Primary_GEA_loci_overlapping_C2",

    "Unique_C2_regions_overlapping_primary_GEA",

    "Primary_GEA_loci_overlapping_XtX_and_C2",

    "Sensitivity_GEA_loci",

    "Sensitivity_GEA_loci_overlapping_XtX",

    "Unique_XtX_regions_overlapping_sensitivity_GEA",

    "Sensitivity_GEA_loci_overlapping_C2",

    "Unique_C2_regions_overlapping_sensitivity_GEA",

    "Sensitivity_GEA_loci_overlapping_XtX_and_C2"
  ),

  Value = as.integer(
    c(

      nrow(
        primary_gea
      ),

      safe_uniqueN(
        gea_xtx_primary$GEA_locus_ID
      ),

      safe_uniqueN(
        gea_xtx_primary$XtX_region_index
      ),

      safe_uniqueN(
        gea_c2_primary$GEA_locus_ID
      ),

      safe_uniqueN(
        gea_c2_primary$C2_region_index
      ),

      nrow(
        triple_primary
      ),

      nrow(
        sensitivity_gea
      ),

      safe_uniqueN(
        gea_xtx_sensitivity$GEA_locus_ID
      ),

      safe_uniqueN(
        gea_xtx_sensitivity$XtX_region_index
      ),

      safe_uniqueN(
        gea_c2_sensitivity$GEA_locus_ID
      ),

      safe_uniqueN(
        gea_c2_sensitivity$C2_region_index
      ),

      nrow(
        triple_sensitivity
      )
    )
  )
)

# ============================================================
# RESTORE FINAL ORDER AFTER MERGES
# ============================================================

gea[
  ,
  Covariate_order :=
    match(
      Covariate,
      EXPECTED_COVARIATES
    )
]

setorder(
  gea,
  Covariate_order,
  Scaffold,
  Region_start
)

# ============================================================
# FINAL TABLE
# ============================================================

final_columns <- c(

  "Analysis_class",
  "Covariate",
  "COVARIABLE",

  "GEA_locus_ID",

  "Scaffold",

  "Region_start",
  "Region_end",
  "Region_size_bp",

  "N_BF20_SNPs",

  "Peak_position",
  "Peak_MRK",
  "Peak_BFdb",
  "Peak_REF",
  "Peak_ALT",

  "Annotation_available",

  "Peak_context",

  "Gene_ID",
  "Gene_symbol",
  "Description",

  "Distance_to_gene_bp",

  "Region_overlaps_gene",

  "Region_overlapping_gene_IDs",
  "Region_overlapping_gene_symbols",

  "Overlaps_XtX",
  "XtX_region_indices",
  "XtX_regions",

  "Overlaps_C2",
  "C2_region_indices",
  "C2_regions",

  "Overlaps_XtX_and_C2"
)

final_columns <- final_columns[
  final_columns %in%
  names(gea)
]

gea_final <- gea[
  ,
  ..final_columns
]

# ============================================================
# FINAL QC BEFORE WRITE
# ============================================================

if (
  nrow(
    gea_final
  ) !=
  210L
) {

  warning(
    "Expected 210 total GEA loci; observed ",
    nrow(
      gea_final
    )
  )
}

if (
  nrow(
    gea_final[
      Analysis_class ==
      "Primary"
    ]
  ) !=
  167L
) {

  warning(
    "Expected 167 primary GEA loci; observed ",
    nrow(
      gea_final[
        Analysis_class ==
        "Primary"
      ]
    )
  )
}

if (
  nrow(
    gea_final[
      Analysis_class ==
      "Sensitivity"
    ]
  ) !=
  43L
) {

  warning(
    "Expected 43 sensitivity GEA loci; observed ",
    nrow(
      gea_final[
        Analysis_class ==
        "Sensitivity"
      ]
    )
  )
}

if (
  sum(
    gea$N_BF20_SNPs
  ) !=
  nrow(
    gea_snps
  )
) {

  stop(
    "Number of BF>20 SNP-covariate associations does not match locus table."
  )
}

# ============================================================
# WRITE FINAL GEA TABLES
# ============================================================

fwrite(
  gea_final,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_candidate_loci_1kb_ALL_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  gea_final[
    Analysis_class ==
    "Primary"
  ],
  file.path(
    OUTDIR,
    "Table_GEA_BF20_candidate_loci_1kb_PRIMARY_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  gea_final[
    Analysis_class ==
    "Sensitivity"
  ],
  file.path(
    OUTDIR,
    "Table_GEA_BF20_candidate_loci_1kb_SENSITIVITY_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

# ============================================================
# WRITE SUMMARY TABLES
# ============================================================

summary_write <- copy(
  summary_by_covariate
)

summary_write[
  ,
  Covariate_order :=
    NULL
]

fwrite(
  summary_write,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_summary_by_covariate_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  summary_by_analysis,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_summary_primary_vs_sensitivity_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  context_summary,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_peak_context_summary_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

# ============================================================
# WRITE XtX OVERLAPS
# ============================================================

fwrite(
  gea_xtx,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_XtX_overlap_ALL_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  gea_xtx_primary,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_XtX_overlap_PRIMARY_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  gea_xtx_sensitivity,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_XtX_overlap_SENSITIVITY_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

# ============================================================
# WRITE C2 OVERLAPS
# ============================================================

fwrite(
  gea_c2,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_C2_overlap_ALL_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  gea_c2_primary,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_C2_overlap_PRIMARY_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  gea_c2_sensitivity,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_C2_overlap_SENSITIVITY_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

# ============================================================
# WRITE TRIPLE OVERLAPS
# ============================================================

fwrite(
  triple_loci,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_XtX_C2_triple_overlap_ALL_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  triple_primary,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_XtX_C2_triple_overlap_PRIMARY_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  triple_sensitivity,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_XtX_C2_triple_overlap_SENSITIVITY_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

# ============================================================
# WRITE OVERLAP SUMMARIES
# ============================================================

fwrite(
  overlap_summary_by_covariate,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_overlap_summary_by_covariate_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  overlap_global_summary,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_overlap_global_summary_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

# ============================================================
# ADD FINAL LOCUS INFO TO SNP TABLE
# ============================================================

gea_snps_final <- merge(
  gea_snps,
  gea[
    ,
    .(

      Covariate,

      Scaffold,

      Cluster_ID,

      GEA_locus_ID,

      Peak_position,

      Peak_BFdb,

      Annotation_available,

      Peak_context,

      Gene_ID,

      Gene_symbol,

      Overlaps_XtX,

      Overlaps_C2,

      Overlaps_XtX_and_C2
    )
  ],
  by = c(
    "Covariate",
    "Scaffold",
    "Cluster_ID"
  ),
  all.x = TRUE,
  sort = FALSE
)

fwrite(
  gea_snps_final,
  file.path(
    OUTDIR,
    "Table_GEA_BF20_all_candidate_SNPs_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = ""
)

# ============================================================
# FINAL QC OUTPUT
# ============================================================

cat("\n=============================================\n")
cat("FINAL GEA BF >20 / 1-kb RESULTS\n")
cat("=============================================\n")

cat(
  "BF threshold: >",
  BF_THRESHOLD,
  " dB\n"
)

cat(
  "Clustering gap:",
  CLUSTER_GAP_BP,
  "bp\n"
)

cat(
  "Gene proximity:",
  GENE_PROXIMITY_BP,
  "bp\n"
)

cat(
  "\nTotal GEA loci:",
  nrow(
    gea_final
  ),
  "\n"
)

cat(
  "Primary GEA loci:",
  nrow(
    gea_final[
      Analysis_class ==
      "Primary"
    ]
  ),
  "\n"
)

cat(
  "Sensitivity loci:",
  nrow(
    gea_final[
      Analysis_class ==
      "Sensitivity"
    ]
  ),
  "\n"
)

cat(
  "Total BF>20 SNP-covariate associations:",
  nrow(
    gea_snps_final
  ),
  "\n"
)

# ============================================================
# ANNOTATION AVAILABILITY
# ============================================================

cat("\n---------------------------------------------\n")
cat("ANNOTATION AVAILABILITY\n")
cat("---------------------------------------------\n")

annotation_qc <- gea[
  ,
  .(
    N_loci =
      .N
  ),
  by = .(
    Analysis_class,
    Annotation_available
  )
]

print(
  annotation_qc
)

# ============================================================
# PEAK CONTEXT
# ============================================================

cat("\n---------------------------------------------\n")
cat("PEAK CONTEXT\n")
cat("---------------------------------------------\n")

print(
  context_summary[
    order(
      Analysis_class,
      Covariate,
      Peak_context
    )
  ]
)

# ============================================================
# SUMMARY BY COVARIATE
# ============================================================

cat("\n---------------------------------------------\n")
cat("SUMMARY BY COVARIATE\n")
cat("---------------------------------------------\n")

print(
  summary_write
)

# ============================================================
# PRIMARY VS SENSITIVITY
# ============================================================

cat("\n---------------------------------------------\n")
cat("PRIMARY VS SENSITIVITY\n")
cat("---------------------------------------------\n")

print(
  summary_by_analysis
)

# ============================================================
# GLOBAL OVERLAPS
# ============================================================

cat("\n---------------------------------------------\n")
cat("GEA x XtX x C2 GLOBAL OVERLAPS\n")
cat("---------------------------------------------\n")

print(
  overlap_global_summary
)

# ============================================================
# BY-COVARIATE OVERLAPS
# ============================================================

cat("\n---------------------------------------------\n")
cat("OVERLAPS BY COVARIATE\n")
cat("---------------------------------------------\n")

print(
  overlap_summary_by_covariate[
    order(
      Analysis_class,
      Covariate
    )
  ]
)

# ============================================================
# PRIMARY XtX OVERLAPS
# ============================================================

cat("\n---------------------------------------------\n")
cat("PRIMARY GEA x XtX OVERLAPS\n")
cat("---------------------------------------------\n")

if (nrow(gea_xtx_primary) == 0L) {

  cat("NONE\n")

} else {

  print(
    gea_xtx_primary[
      order(
        Scaffold,
        GEA_start
      )
    ]
  )
}

# ============================================================
# PRIMARY C2 OVERLAPS
# ============================================================

cat("\n---------------------------------------------\n")
cat("PRIMARY GEA x C2 OVERLAPS\n")
cat("---------------------------------------------\n")

if (nrow(gea_c2_primary) == 0L) {

  cat("NONE\n")

} else {

  print(
    gea_c2_primary[
      order(
        Scaffold,
        GEA_start
      )
    ]
  )
}

# ============================================================
# TRIPLE OVERLAPS
# ============================================================

cat("\n---------------------------------------------\n")
cat("PRIMARY GEA x XtX x C2 TRIPLE OVERLAPS\n")
cat("---------------------------------------------\n")

if (nrow(triple_primary) == 0L) {

  cat("NONE\n")

} else {

  print(
    triple_primary[
      ,
      .(

        Covariate,

        Scaffold,

        GEA_locus_ID,

        Region_start,

        Region_end,

        Peak_position,

        Peak_BFdb,

        Gene_ID,

        Gene_symbol,

        Peak_context,

        XtX_regions,

        C2_regions
      )
    ]
  )
}

# ============================================================
# TOP 25 PRIMARY GEA LOCI
# ============================================================

cat("\n---------------------------------------------\n")
cat("TOP 25 PRIMARY GEA LOCI BY BF\n")
cat("---------------------------------------------\n")

top_primary <- gea_final[
  Analysis_class ==
  "Primary"
][
  order(
    -Peak_BFdb
  )
]

top_primary <- head(
  top_primary,
  25
)

print(
  top_primary[
    ,
    .(

      Covariate,

      Scaffold,

      Region_start,

      Region_end,

      N_BF20_SNPs,

      Peak_position,

      Peak_BFdb,

      Annotation_available,

      Peak_context,

      Gene_ID,

      Gene_symbol,

      Distance_to_gene_bp,

      Overlaps_XtX,

      Overlaps_C2,

      Overlaps_XtX_and_C2
    )
  ]
)

cat("\n=============================================\n")
cat("OUTPUT DIRECTORY\n")
cat("=============================================\n")

cat(
  OUTDIR,
  "\n"
)

cat("\nDONE\n")

