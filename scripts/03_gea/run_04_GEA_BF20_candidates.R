#!/usr/bin/env Rscript

############################################################
## run_04_GEA_BF20_candidates.R
##
## Cactoblastis cactorum
## Final genotype-environment association analysis
##
## Dataset:
##   all-SNP dataset only
##
## Input statistic:
##   median BF(dB) across three converged BayPass runs
##
## Candidate definition:
##   BF(dB) > 20
##
## Candidate loci:
##   significant SNPs for the same covariate and scaffold
##   are grouped when separated by <= 1 kb.
##
## GEA is NOT analysed using local score.
##
## Final covariates:
##   TempPC1
##   PrecipPC1
##   AridityIndex
##   Elevation
##   Latitude
##   Longitude
############################################################


suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(IRanges)
})


############################################################
## 0. SETTINGS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

BF_FILE <- file.path(
    BASE,
    "selection_scan",
    "prepared_statistics",
    "all_snps",
    "GEA_BF_median3seeds.tsv"
)

GFF_FILE <- file.path(
    "/home/dpoveda/micetral/Population_snps",
    "Ccac_functional-annotation",
    "02_agat",
    "Ccac",
    "out_abril4",
    "Ccac.annotations.corrected.gff3"
)

## Final XtX/C2 candidate tables generated upstream
XTX_C2_DIR <- file.path(
    BASE,
    "selection_scan",
    "candidate_annotation",
    "XtX_C2"
)

XTX_C2_FILE <- file.path(
    XTX_C2_DIR,
    "XtX_C2_candidate_regions_annotation.tsv"
)

OUTDIR <- file.path(
    BASE,
    "selection_scan",
    "GEA_BF20_final"
)

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


############################################################
## 1. FINAL ANALYSIS DEFINITIONS
############################################################

BF_THRESHOLD <- 20

CLUSTER_GAP_BP <- 1000L

GENE_PROXIMITY_BP <- 5000L


COVARIATES <- c(
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

clean_scaffold <- function(x) {

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

    out[out == ""] <- NA_character_

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

    x <- unique(
        as.character(
            x[
                !is.na(x) &
                x != ""
            ]
        )
    )

    if (length(x) == 0L) {
        return(NA_character_)
    }

    paste(
        sort(x),
        collapse = ";"
    )
}


classify_peak_context <- function(
    in_cds,
    in_exon,
    in_gene,
    near_gene
) {

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


############################################################
## 3. READ FINAL MEDIAN BF TABLE
############################################################

if (!file.exists(BF_FILE)) {

    stop(
        "Final GEA BF table not found:\n",
        BF_FILE
    )
}


bf <- fread(
    BF_FILE,
    showProgress = FALSE
)


required <- c(
    "BayPass_row",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "M_P",
    "COVARIABLE",
    "Covariate",
    "BF_median"
)


missing_cols <- setdiff(
    required,
    names(bf)
)


if (length(missing_cols) > 0L) {

    stop(
        "Missing BF columns: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


############################################################
## 4. RETAIN FINAL SIX COVARIATES
############################################################

bf <- bf[
    Covariate %in%
        COVARIATES
]


if (
    !setequal(
        unique(bf$Covariate),
        COVARIATES
    )
) {

    stop(
        "Final covariate set does not match expected manuscript covariates."
    )
}


bf[
    ,
    Scaffold :=
        clean_scaffold(
            CHROM
        )
]

bf[
    ,
    POS :=
        as.integer(POS)
]

bf[
    ,
    BayPass_row :=
        as.integer(BayPass_row)
]

bf[
    ,
    BF_median :=
        as.numeric(BF_median)
]


############################################################
## 5. BF > 20 CANDIDATE SNPs
############################################################

bf20 <- bf[
    BF_median >
        BF_THRESHOLD
]


bf20[
    ,
    Covariate_order :=
        match(
            Covariate,
            COVARIATES
        )
]


setorder(
    bf20,
    Covariate_order,
    Scaffold,
    POS
)


cat(
    "BF >20 SNP-covariate associations: ",
    nrow(bf20),
    "\n",
    sep = ""
)


############################################################
## 6. CLUSTER BF >20 SNPs WITHIN 1 kb
##
## Clustering is performed independently for each
## covariate and scaffold.
############################################################

bf20[
    ,
    Gap_from_previous :=
        c(
            NA_integer_,
            diff(POS)
        ),
    by = .(
        Covariate,
        Scaffold
    )
]


bf20[
    ,
    New_cluster :=
        fifelse(
            is.na(
                Gap_from_previous
            ) |
            Gap_from_previous >
                CLUSTER_GAP_BP,
            1L,
            0L
        )
]


bf20[
    ,
    Cluster_number :=
        cumsum(
            New_cluster
        ),
    by = .(
        Covariate,
        Scaffold
    )
]


bf20[
    ,
    Cluster_ID :=
        paste(
            Covariate,
            Scaffold,
            Cluster_number,
            sep = "_"
        )
]


############################################################
## 7. CREATE CANDIDATE LOCI
############################################################

gea <- bf20[
    ,
    {

        peak_index <- which.max(
            BF_median
        )

        .(
            Region_start =
                min(POS),

            Region_end =
                max(POS),

            Region_size_bp =
                max(POS) -
                min(POS) +
                1L,

            N_BF20_SNPs =
                .N,

            Peak_position =
                POS[
                    peak_index
                ],

            Peak_MRK =
                BayPass_row[
                    peak_index
                ],

            Peak_BFdb =
                BF_median[
                    peak_index
                ],

            Peak_REF =
                REF[
                    peak_index
                ],

            Peak_ALT =
                ALT[
                    peak_index
                ]
        )
    },
    by = .(
        COVARIABLE,
        Covariate,
        Scaffold,
        Cluster_ID,
        Cluster_number
    )
]


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
    )
) {

    stop(
        "GEA_locus_ID values are not unique."
    )
}


############################################################
## 8. READ REFERENCE ANNOTATION
############################################################

if (!file.exists(GFF_FILE)) {

    stop(
        "GFF3 file not found:\n",
        GFF_FILE
    )
}


gff <- fread(
    GFF_FILE,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    quote = "",
    fill = TRUE,
    showProgress = FALSE
)


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
        as.integer(start)
]

gff[
    ,
    end :=
        as.integer(end)
]


gff <- gff[
    !is.na(Scaffold) &
    !is.na(start) &
    !is.na(end) &
    end >= start
]


############################################################
## 9. BUILD GENE TABLE
############################################################

genes <- gff[
    feature == "gene"
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


############################################################
## 10. GENOMIC RANGES
############################################################

gene_gr <- GRanges(
    seqnames =
        genes$Scaffold,
    ranges =
        IRanges(
            start =
                genes$start,
            end =
                genes$end
        )
)

mcols(gene_gr)$Gene_ID <-
    genes$Gene_ID

mcols(gene_gr)$Gene_symbol <-
    genes$Gene_symbol

mcols(gene_gr)$Description <-
    genes$Description


exons <- gff[
    feature == "exon"
]

exon_gr <- GRanges(
    seqnames =
        exons$Scaffold,
    ranges =
        IRanges(
            start =
                exons$start,
            end =
                exons$end
        )
)


cds <- gff[
    feature == "CDS"
]

cds_gr <- GRanges(
    seqnames =
        cds$Scaffold,
    ranges =
        IRanges(
            start =
                cds$start,
            end =
                cds$end
        )
)


region_gr <- GRanges(
    seqnames =
        gea$Scaffold,
    ranges =
        IRanges(
            start =
                gea$Region_start,
            end =
                gea$Region_end
        )
)


peak_gr <- GRanges(
    seqnames =
        gea$Scaffold,
    ranges =
        IRanges(
            start =
                gea$Peak_position,
            end =
                gea$Peak_position
        )
)


############################################################
## 11. REGION-GENE OVERLAPS
############################################################

gea[
    ,
    Region_overlaps_gene :=
        overlapsAny(
            region_gr,
            gene_gr,
            ignore.strand = TRUE
        )
]


region_hits <- findOverlaps(
    region_gr,
    gene_gr,
    ignore.strand = TRUE
)


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


if (length(region_hits) > 0L) {

    hit_dt <- data.table(

        locus_index =
            queryHits(
                region_hits
            ),

        Gene_ID =
            mcols(gene_gr)$Gene_ID[
                subjectHits(
                    region_hits
                )
            ],

        Gene_symbol =
            mcols(gene_gr)$Gene_symbol[
                subjectHits(
                    region_hits
                )
            ]
    )


    hit_summary <- hit_dt[
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
        by = locus_index
    ]


    gea[
        hit_summary$locus_index,
        Region_overlapping_gene_IDs :=
            hit_summary$
                Region_overlapping_gene_IDs
    ]


    gea[
        hit_summary$locus_index,
        Region_overlapping_gene_symbols :=
            hit_summary$
                Region_overlapping_gene_symbols
    ]
}


############################################################
## 12. PEAK GENOMIC CONTEXT
############################################################

gea[
    ,
    Peak_in_gene :=
        overlapsAny(
            peak_gr,
            gene_gr,
            ignore.strand = TRUE
        )
]

gea[
    ,
    Peak_in_exon :=
        overlapsAny(
            peak_gr,
            exon_gr,
            ignore.strand = TRUE
        )
]

gea[
    ,
    Peak_in_CDS :=
        overlapsAny(
            peak_gr,
            cds_gr,
            ignore.strand = TRUE
        )
]


############################################################
## 13. NEAREST GENE TO PEAK
############################################################

gea[, Gene_ID := NA_character_]
gea[, Gene_symbol := NA_character_]
gea[, Description := NA_character_]
gea[, Distance_to_gene_bp := NA_integer_]


near_hits <- distanceToNearest(
    peak_gr,
    gene_gr,
    ignore.strand = TRUE
)


if (length(near_hits) > 0L) {

    q <- queryHits(
        near_hits
    )

    s <- subjectHits(
        near_hits
    )

    d <- as.integer(
        mcols(
            near_hits
        )$distance
    )


    keep <-
        !is.na(d) &
        d <=
            GENE_PROXIMITY_BP


    q <- q[keep]
    s <- s[keep]
    d <- d[keep]


    gea[
        q,
        Gene_ID :=
            mcols(gene_gr)$Gene_ID[
                s
            ]
    ]

    gea[
        q,
        Gene_symbol :=
            mcols(gene_gr)$Gene_symbol[
                s
            ]
    ]

    gea[
        q,
        Description :=
            mcols(gene_gr)$Description[
                s
            ]
    ]

    gea[
        q,
        Distance_to_gene_bp :=
            d
    ]
}


gea[
    ,
    Gene_ID :=
        sub(
            "^gene-",
            "",
            Gene_ID
        )
]


############################################################
## 14. PEAK CONTEXT
############################################################

gea[
    ,
    Peak_context :=
        mapply(
            classify_peak_context,

            Peak_in_CDS,

            Peak_in_exon,

            Peak_in_gene,

            !is.na(
                Gene_ID
            ),

            USE.NAMES = FALSE
        )
]


############################################################
## 15. ORDER FINAL LOCI
############################################################

gea[
    ,
    Covariate_order :=
        match(
            Covariate,
            COVARIATES
        )
]

setorder(
    gea,
    Covariate_order,
    Scaffold,
    Region_start
)


############################################################
## 16. LOAD FINAL XTX / C2 REGIONS
############################################################

if (!file.exists(XTX_C2_FILE)) {

    stop(
        "XtX/C2 candidate annotation table not found:\n",
        XTX_C2_FILE
    )
}


selection_regions <- fread(
    XTX_C2_FILE,
    showProgress = FALSE
)


xtx <- selection_regions[
    Analysis == "XtX"
]

c2 <- selection_regions[
    Analysis == "C2_Host"
]


############################################################
## 17. PHYSICAL OVERLAP WITH XTX AND C2
############################################################

gea[
    ,
    Overlaps_XtX := FALSE
]

gea[
    ,
    Overlaps_C2 := FALSE
]


gea_gr <- GRanges(
    seqnames =
        gea$Scaffold,
    ranges =
        IRanges(
            start =
                gea$Region_start,
            end =
                gea$Region_end
        )
)


if (nrow(xtx) > 0L) {

    xtx_gr <- GRanges(
        seqnames =
            clean_scaffold(
                xtx$CHROM
            ),
        ranges =
            IRanges(
                start =
                    xtx$Start,
                end =
                    xtx$End
            )
    )

    xtx_hits <- findOverlaps(
        gea_gr,
        xtx_gr,
        ignore.strand = TRUE
    )

    if (length(xtx_hits) > 0L) {

        gea[
            unique(
                queryHits(
                    xtx_hits
                )
            ),
            Overlaps_XtX :=
                TRUE
        ]
    }
}


if (nrow(c2) > 0L) {

    c2_gr <- GRanges(
        seqnames =
            clean_scaffold(
                c2$CHROM
            ),
        ranges =
            IRanges(
                start =
                    c2$Start,
                end =
                    c2$End
            )
    )

    c2_hits <- findOverlaps(
        gea_gr,
        c2_gr,
        ignore.strand = TRUE
    )

    if (length(c2_hits) > 0L) {

        gea[
            unique(
                queryHits(
                    c2_hits
                )
            ),
            Overlaps_C2 :=
                TRUE
        ]
    }
}


gea[
    ,
    Overlaps_XtX_and_C2 :=
        Overlaps_XtX &
        Overlaps_C2
]


############################################################
## 18. SUMMARY BY COVARIATE
############################################################

summary_by_covariate <- gea[
    ,
    .(
        N_BF20_SNPs =
            sum(
                N_BF20_SNPs
            ),

        N_candidate_loci =
            .N,

        N_scaffolds =
            uniqueN(
                Scaffold
            ),

        N_gene_associated_peaks_5kb =
            sum(
                !is.na(
                    Gene_ID
                )
            ),

        N_XtX_overlaps =
            sum(
                Overlaps_XtX
            ),

        N_C2_overlaps =
            sum(
                Overlaps_C2
            ),

        Max_BFdb =
            max(
                Peak_BFdb
            )
    ),
    by = .(
        COVARIABLE,
        Covariate
    )
]


summary_by_covariate[
    ,
    Covariate_order :=
        match(
            Covariate,
            COVARIATES
        )
]

setorder(
    summary_by_covariate,
    Covariate_order
)


############################################################
## 19. FINAL OUTPUT TABLE
############################################################

final_columns <- c(

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

    "Peak_context",

    "Gene_ID",
    "Gene_symbol",
    "Description",

    "Distance_to_gene_bp",

    "Region_overlaps_gene",

    "Region_overlapping_gene_IDs",
    "Region_overlapping_gene_symbols",

    "Overlaps_XtX",
    "Overlaps_C2",
    "Overlaps_XtX_and_C2"
)


gea_final <- gea[
    ,
    ..final_columns
]


############################################################
## 20. WRITE OUTPUTS
############################################################

fwrite(
    bf20,
    file.path(
        OUTDIR,
        "GEA_BF20_candidate_SNPs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
)


fwrite(
    gea_final,
    file.path(
        OUTDIR,
        "GEA_BF20_candidate_loci.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
)


fwrite(
    summary_by_covariate[
        ,
        Covariate_order := NULL
    ],
    file.path(
        OUTDIR,
        "GEA_BF20_summary_by_covariate.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
)


fwrite(
    gea_final[
        Overlaps_XtX == TRUE
    ],
    file.path(
        OUTDIR,
        "GEA_BF20_overlap_XtX.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
)


fwrite(
    gea_final[
        Overlaps_C2 == TRUE
    ],
    file.path(
        OUTDIR,
        "GEA_BF20_overlap_C2.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
)


############################################################
## 21. FINAL QC
############################################################

cat(
    "\n====================================================\n"
)

cat(
    "FINAL GEA BF ANALYSIS COMPLETE\n"
)

cat(
    "====================================================\n"
)

cat(
    "\nBF threshold: >",
    BF_THRESHOLD,
    " dB\n",
    sep = ""
)

cat(
    "Clustering distance: <= ",
    CLUSTER_GAP_BP,
    " bp\n",
    sep = ""
)

cat(
    "BF >20 SNP-covariate associations: ",
    nrow(bf20),
    "\n",
    sep = ""
)

cat(
    "Candidate loci: ",
    nrow(gea_final),
    "\n",
    sep = ""
)

cat(
    "\nSummary by covariate:\n"
)

print(
    summary_by_covariate
)

cat(
    "\nOutput directory:\n",
    OUTDIR,
    "\n"
)

cat(
    "\nDONE\n"
)
