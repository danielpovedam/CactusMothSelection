#!/usr/bin/env Rscript

############################################################
## run_03_annotate_XtX_C2_candidates.R
##
## Cactoblastis cactorum
## Functional annotation of final XtX and C2 candidate regions
##
## Dataset:
##   all-SNP dataset only
##
## Analyses:
##   XtX
##   C2 host-use contrast
##
## GEA candidates are annotated separately and are not
## processed in this script.
############################################################


suppressPackageStartupMessages({
    library(data.table)
    library(stringr)
    library(GenomicRanges)
})


############################################################
## 0. SETTINGS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

INPUT_FILE <- file.path(
    BASE,
    "selection_scan/candidate_regions/candidate_regions_detailed.tsv"
)

FUN_BASE <- paste0(
    "/home/dpoveda/micetral/Population_snps/",
    "Ccac_functional-annotation"
)

GFF_FILE <- file.path(
    FUN_BASE,
    "02_agat/Ccac/out_abril4/Ccac.annotations.corrected.gff3"
)

AGAT_DIR <- file.path(
    FUN_BASE,
    "02_agat/Ccac/out_abril4"
)

FANTASIA_GO <- file.path(
    FUN_BASE,
    "01_mapping2db/FANTASIA/Ccac_fantasia/Ccac.proteins_topgo.txt"
)

OUT_DIR <- file.path(
    BASE,
    "selection_scan/candidate_annotation/XtX_C2"
)

dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

MAX_DIST <- 5000L


############################################################
## 1. HELPERS
############################################################

extract_attr <- function(x, key) {

    m <- str_match(
        x,
        paste0("(^|;)", key, "=([^;]+)")
    )

    ans <- m[, 3]

    ans[
        is.na(ans) |
        ans == ""
    ] <- NA_character_

    ans
}


collapse_unique <- function(x) {

    x <- unique(
        x[
            !is.na(x) &
            x != "" &
            x != "-"
        ]
    )

    if (length(x) == 0) {
        return(NA_character_)
    }

    paste(
        sort(x),
        collapse = ";"
    )
}


gene_from_transcript <- function(x) {

    sub(
        "\\.t[0-9]+$",
        "",
        x
    )
}


read_annotation_map <- function(
    file,
    column_name
) {

    if (!file.exists(file)) {

        return(
            data.table(
                gene_id = character(),
                value = character()
            )
        )
    }

    x <- fread(
        file,
        sep = "\t",
        header = FALSE,
        fill = TRUE,
        quote = "",
        showProgress = FALSE
    )

    if (ncol(x) < 2) {

        return(
            data.table(
                gene_id = character(),
                value = character()
            )
        )
    }

    x <- x[, 1:2]

    setnames(
        x,
        c(
            "transcript_id",
            "value"
        )
    )

    x[
        ,
        transcript_id :=
            as.character(transcript_id)
    ]

    x[
        ,
        value :=
            as.character(value)
    ]

    x[
        ,
        gene_id :=
            gene_from_transcript(
                transcript_id
            )
    ]

    ans <- x[
        ,
        .(
            value =
                collapse_unique(value)
        ),
        by = gene_id
    ]

    setnames(
        ans,
        "value",
        column_name
    )

    ans
}


############################################################
## 2. INPUT CHECKS
############################################################

if (!file.exists(INPUT_FILE)) {

    stop(
        "Candidate-region file not found:\n",
        INPUT_FILE
    )
}

if (!file.exists(GFF_FILE)) {

    stop(
        "Functional GFF3 not found:\n",
        GFF_FILE
    )
}


############################################################
## 3. READ FUNCTIONALLY ANNOTATED GFF3
############################################################

cat(
    "Loading functionally annotated GFF3...\n"
)

gff <- fread(
    GFF_FILE,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    quote = "",
    fill = TRUE,
    showProgress = FALSE,
    col.names = c(
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

gff[, start := as.integer(start)]
gff[, end   := as.integer(end)]

gff <- gff[
    !is.na(seqid) &
    !is.na(start) &
    !is.na(end) &
    end >= start
]


############################################################
## 4. GENES AND TRANSCRIPTS
############################################################

genes <- gff[
    feature == "gene"
]

genes[
    ,
    gene_id :=
        extract_attr(
            attributes,
            "ID"
        )
]

genes[
    ,
    gene_name_gff :=
        extract_attr(
            attributes,
            "Name"
        )
]

genes[
    is.na(gene_name_gff),
    gene_name_gff := gene_id
]


tx <- gff[
    feature == "transcript"
]

tx[
    ,
    transcript_id :=
        extract_attr(
            attributes,
            "ID"
        )
]

tx[
    ,
    gene_id :=
        extract_attr(
            attributes,
            "Parent"
        )
]

tx[
    ,
    product :=
        extract_attr(
            attributes,
            "product"
        )
]

tx[
    ,
    dbxref :=
        extract_attr(
            attributes,
            "Dbxref"
        )
]

tx[
    ,
    ontology_term :=
        extract_attr(
            attributes,
            "Ontology_term"
        )
]


############################################################
## 5. COLLAPSE TRANSCRIPT ANNOTATION TO GENE
############################################################

tx_gene <- tx[
    ,
    .(
        transcript_ids =
            collapse_unique(
                transcript_id
            ),

        products =
            collapse_unique(
                product
            ),

        dbxref =
            collapse_unique(
                dbxref
            ),

        ontology_terms =
            collapse_unique(
                ontology_term
            )
    ),
    by = gene_id
]


tx_product <- tx[
    !is.na(product) &
    product != "" &
    product != "hypothetical protein",
    .(
        functional_product =
            collapse_unique(
                product
            )
    ),
    by = gene_id
]

tx_gene <- merge(
    tx_gene,
    tx_product,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)

tx_gene[
    is.na(functional_product) &
    !is.na(products) &
    grepl(
        "hypothetical protein",
        products,
        fixed = TRUE
    ),
    functional_product :=
        "hypothetical protein"
]


############################################################
## 6. FUNCTIONAL ANNOTATION MAPS
############################################################

GO_MAP <- read_annotation_map(
    file.path(
        AGAT_DIR,
        "GO.txt"
    ),
    "GO"
)

INTERPRO_MAP <- read_annotation_map(
    file.path(
        AGAT_DIR,
        "InterPro.txt"
    ),
    "InterPro"
)

PFAM_MAP <- read_annotation_map(
    file.path(
        AGAT_DIR,
        "Pfam.txt"
    ),
    "Pfam"
)

PANTHER_MAP <- read_annotation_map(
    file.path(
        AGAT_DIR,
        "PANTHER.txt"
    ),
    "PANTHER"
)

SIGNALP_MAP <- read_annotation_map(
    file.path(
        AGAT_DIR,
        "SignalP_EUK.txt"
    ),
    "SignalP"
)

TMHMM_MAP <- read_annotation_map(
    file.path(
        AGAT_DIR,
        "TMHMM.txt"
    ),
    "TMHMM"
)

GENE3D_MAP <- read_annotation_map(
    file.path(
        AGAT_DIR,
        "Gene3D.txt"
    ),
    "Gene3D"
)

FANTASIA_MAP <- read_annotation_map(
    FANTASIA_GO,
    "FANTASIA_GO"
)


############################################################
## 7. MASTER GENE ANNOTATION TABLE
############################################################

gene_info <- genes[
    ,
    .(
        gene_id,
        gene_name_gff,
        seqid,
        gene_start = start,
        gene_end = end,
        gene_strand = strand
    )
]

maps <- list(
    tx_gene,
    GO_MAP,
    INTERPRO_MAP,
    PFAM_MAP,
    PANTHER_MAP,
    SIGNALP_MAP,
    TMHMM_MAP,
    GENE3D_MAP,
    FANTASIA_MAP
)

for (m in maps) {

    gene_info <- merge(
        gene_info,
        m,
        by = "gene_id",
        all.x = TRUE,
        sort = FALSE
    )
}


############################################################
## 8. GENOMIC RANGES
############################################################

gene_gr <- GRanges(
    seqnames = genes$seqid,
    ranges = IRanges(
        start = genes$start,
        end = genes$end
    ),
    strand = ifelse(
        genes$strand %in% c("+", "-"),
        genes$strand,
        "*"
    )
)

mcols(gene_gr)$gene_id <-
    genes$gene_id


exons <- gff[
    feature == "exon"
]

exon_gr <- GRanges(
    seqnames = exons$seqid,
    ranges = IRanges(
        start = exons$start,
        end = exons$end
    )
)


cds <- gff[
    feature == "CDS"
]

cds_gr <- GRanges(
    seqnames = cds$seqid,
    ranges = IRanges(
        start = cds$start,
        end = cds$end
    )
)


introns <- gff[
    feature == "intron"
]

intron_gr <- GRanges(
    seqnames = introns$seqid,
    ranges = IRanges(
        start = introns$start,
        end = introns$end
    )
)


############################################################
## 9. READ FINAL CANDIDATE REGIONS
############################################################

regions <- fread(
    INPUT_FILE,
    showProgress = FALSE
)


############################################################
## IMPORTANT:
## retain only the final all-SNP XtX and C2 analyses
############################################################

regions <- regions[
    Dataset == "all_snps" &
    Analysis %in% c(
        "XtX",
        "C2_Host"
    )
]

if (nrow(regions) == 0) {

    stop(
        paste0(
            "No all-SNP XtX/C2 candidate regions ",
            "were found in the input table."
        )
    )
}


required_cols <- c(
    "Dataset",
    "Analysis",
    "Region_ID",
    "CHROM",
    "Start",
    "End",
    "Window_size_bp",
    "N_SNPs",
    "Statistic_peak_position"
)

missing_cols <- setdiff(
    required_cols,
    names(regions)
)

if (length(missing_cols) > 0) {

    stop(
        "Missing required candidate-region columns: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


regions[, Start := as.integer(Start)]
regions[, End := as.integer(End)]

regions[
    ,
    Window_size_bp :=
        as.numeric(Window_size_bp)
]

regions[
    ,
    N_SNPs :=
        as.integer(N_SNPs)
]

regions[
    ,
    Statistic_peak_position :=
        as.integer(
            Statistic_peak_position
        )
]


cat(
    "Final XtX/C2 candidate regions: ",
    nrow(regions),
    "\n",
    sep = ""
)


############################################################
## 10. VERIFY SCAFFOLD NAMES
############################################################

missing_chr <- setdiff(
    unique(regions$CHROM),
    unique(gff$seqid)
)

if (length(missing_chr) > 0) {

    stop(
        "Candidate scaffolds missing from GFF3: ",
        paste(
            missing_chr,
            collapse = ", "
        )
    )
}


############################################################
## 11. CANDIDATE-WINDOW RANGES
############################################################

window_gr <- GRanges(
    seqnames = regions$CHROM,
    ranges = IRanges(
        start = regions$Start,
        end = regions$End
    )
)


regions[
    ,
    Window_in_gene :=
        overlapsAny(
            window_gr,
            gene_gr,
            ignore.strand = TRUE
        )
]

regions[
    ,
    Window_in_exon :=
        overlapsAny(
            window_gr,
            exon_gr,
            ignore.strand = TRUE
        )
]

regions[
    ,
    Window_in_CDS :=
        overlapsAny(
            window_gr,
            cds_gr,
            ignore.strand = TRUE
        )
]

regions[
    ,
    Window_in_intron :=
        overlapsAny(
            window_gr,
            intron_gr,
            ignore.strand = TRUE
        )
]


############################################################
## 12. GENES OVERLAPPING EACH CANDIDATE REGION
############################################################

hits <- findOverlaps(
    window_gr,
    gene_gr,
    ignore.strand = TRUE
)


if (length(hits) > 0) {

    region_gene <- data.table(
        region_index =
            queryHits(hits),

        gene_index =
            subjectHits(hits)
    )


    region_gene[
        ,
        Region_ID :=
            regions$Region_ID[
                region_index
            ]
    ]

    region_gene[
        ,
        Analysis :=
            regions$Analysis[
                region_index
            ]
    ]

    region_gene[
        ,
        CHROM :=
            regions$CHROM[
                region_index
            ]
    ]

    region_gene[
        ,
        Region_start :=
            regions$Start[
                region_index
            ]
    ]

    region_gene[
        ,
        Region_end :=
            regions$End[
                region_index
            ]
    ]

    region_gene[
        ,
        Window_size_bp :=
            regions$Window_size_bp[
                region_index
            ]
    ]

    region_gene[
        ,
        N_SNPs :=
            regions$N_SNPs[
                region_index
            ]
    ]

    region_gene[
        ,
        Statistic_peak_position :=
            regions$Statistic_peak_position[
                region_index
            ]
    ]


    if (
        "Statistic_peak_value" %in%
        names(regions)
    ) {

        region_gene[
            ,
            Statistic_peak_value :=
                regions$Statistic_peak_value[
                    region_index
                ]
        ]
    }


    if (
        "Peak_SNP" %in%
        names(regions)
    ) {

        region_gene[
            ,
            Peak_SNP :=
                regions$Peak_SNP[
                    region_index
                ]
        ]
    }


    region_gene[
        ,
        gene_id :=
            mcols(gene_gr)$gene_id[
                gene_index
            ]
    ]


    region_gene <- merge(
        region_gene,
        gene_info,
        by = "gene_id",
        all.x = TRUE,
        sort = FALSE
    )


    region_gene[
        ,
        overlap_start :=
            pmax(
                Region_start,
                gene_start
            )
    ]

    region_gene[
        ,
        overlap_end :=
            pmin(
                Region_end,
                gene_end
            )
    ]

    region_gene[
        ,
        overlap_bp :=
            overlap_end -
            overlap_start +
            1L
    ]

    region_gene[
        ,
        Fraction_gene_overlapped :=
            overlap_bp /
            (
                gene_end -
                gene_start +
                1
            )
    ]

    region_gene[
        ,
        Fraction_window_overlapping_gene :=
            overlap_bp /
            Window_size_bp
    ]


    region_gene[
        ,
        c(
            "region_index",
            "gene_index"
        ) := NULL
    ]

} else {

    region_gene <- data.table()
}


############################################################
## 13. COLLAPSE GENES PER REGION
############################################################

if (nrow(region_gene) > 0) {

    collapsed <- region_gene[
        ,
        .(
            N_overlapping_genes =
                uniqueN(gene_id),

            Overlapping_gene_IDs =
                collapse_unique(
                    gene_id
                ),

            Functional_products =
                collapse_unique(
                    functional_product
                )
        ),
        by = Region_ID
    ]

    regions <- merge(
        regions,
        collapsed,
        by = "Region_ID",
        all.x = TRUE,
        sort = FALSE
    )

} else {

    regions[
        ,
        N_overlapping_genes := 0L
    ]

    regions[
        ,
        Overlapping_gene_IDs :=
            NA_character_
    ]

    regions[
        ,
        Functional_products :=
            NA_character_
    ]
}


regions[
    is.na(N_overlapping_genes),
    N_overlapping_genes := 0L
]


############################################################
## 14. NEAREST GENE TO EACH REGION
############################################################

nearest_win <- distanceToNearest(
    window_gr,
    gene_gr,
    ignore.strand = TRUE
)

regions[
    ,
    Nearest_gene_window :=
        NA_character_
]

regions[
    ,
    Distance_to_nearest_gene_window :=
        NA_integer_
]


if (length(nearest_win) > 0) {

    regions$Nearest_gene_window[
        queryHits(nearest_win)
    ] <- mcols(gene_gr)$gene_id[
        subjectHits(nearest_win)
    ]

    regions$Distance_to_nearest_gene_window[
        queryHits(nearest_win)
    ] <- as.integer(
        mcols(nearest_win)$distance
    )
}


regions[
    ,
    Window_has_gene_within_5kb :=
        (
            N_overlapping_genes > 0
        ) |
        (
            !is.na(
                Distance_to_nearest_gene_window
            ) &
            Distance_to_nearest_gene_window <=
                MAX_DIST
        )
]


############################################################
## 15. ANNOTATE CANDIDATE PEAK POSITION
############################################################

regions[, Peak_context := NA_character_]
regions[, Peak_nearest_gene := NA_character_]
regions[, Peak_distance_to_gene := NA_integer_]

regions[, Peak_in_gene := FALSE]
regions[, Peak_in_exon := FALSE]
regions[, Peak_in_CDS := FALSE]
regions[, Peak_in_intron := FALSE]


peak_ok <- !is.na(
    regions$Statistic_peak_position
)


if (any(peak_ok)) {

    peak_gr <- GRanges(

        seqnames =
            regions$CHROM[
                peak_ok
            ],

        ranges = IRanges(

            start =
                regions$Statistic_peak_position[
                    peak_ok
                ],

            end =
                regions$Statistic_peak_position[
                    peak_ok
                ]
        )
    )


    in_gene <- overlapsAny(
        peak_gr,
        gene_gr,
        ignore.strand = TRUE
    )

    in_exon <- overlapsAny(
        peak_gr,
        exon_gr,
        ignore.strand = TRUE
    )

    in_cds <- overlapsAny(
        peak_gr,
        cds_gr,
        ignore.strand = TRUE
    )

    in_intron <- overlapsAny(
        peak_gr,
        intron_gr,
        ignore.strand = TRUE
    )


    nearest_peak <- distanceToNearest(
        peak_gr,
        gene_gr,
        ignore.strand = TRUE
    )


    idx <- which(
        peak_ok
    )


    regions$Peak_in_gene[idx] <-
        in_gene

    regions$Peak_in_exon[idx] <-
        in_exon

    regions$Peak_in_CDS[idx] <-
        in_cds

    regions$Peak_in_intron[idx] <-
        in_intron


    d <- rep(
        NA_integer_,
        length(idx)
    )

    nearest_gene <- rep(
        NA_character_,
        length(idx)
    )


    if (length(nearest_peak) > 0) {

        d[
            queryHits(nearest_peak)
        ] <- as.integer(
            mcols(nearest_peak)$distance
        )

        nearest_gene[
            queryHits(nearest_peak)
        ] <- mcols(gene_gr)$gene_id[
            subjectHits(nearest_peak)
        ]
    }


    ########################################################
    ## Peak-context hierarchy
    ########################################################

    context <- rep(
        "distal_intergenic",
        length(idx)
    )

    context[
        !is.na(d) &
        d <= MAX_DIST
    ] <- "gene_proximal_5kb"

    context[
        in_gene
    ] <- "genic_other"

    context[
        in_intron
    ] <- "intronic"

    context[
        in_exon
    ] <- "exonic_non_coding"

    context[
        in_cds
    ] <- "coding"


    regions$Peak_context[idx] <-
        context

    regions$Peak_nearest_gene[idx] <-
        nearest_gene

    regions$Peak_distance_to_gene[idx] <-
        d
}


############################################################
## 16. CANDIDATE-GENE TABLE
############################################################

if (nrow(region_gene) > 0) {

    keep_gene_cols <- c(

        "Analysis",
        "Region_ID",
        "CHROM",
        "Region_start",
        "Region_end",
        "Window_size_bp",
        "N_SNPs",
        "Statistic_peak_position",

        "gene_id",
        "gene_name_gff",
        "gene_start",
        "gene_end",
        "gene_strand",

        "transcript_ids",
        "functional_product",
        "products",
        "dbxref",
        "ontology_terms",

        "GO",
        "FANTASIA_GO",
        "InterPro",
        "Pfam",
        "PANTHER",
        "Gene3D",
        "SignalP",
        "TMHMM",

        "overlap_start",
        "overlap_end",
        "overlap_bp",
        "Fraction_gene_overlapped",
        "Fraction_window_overlapping_gene"
    )

    keep_gene_cols <- keep_gene_cols[
        keep_gene_cols %in%
        names(region_gene)
    ]

    candidate_genes <- unique(
        region_gene[
            ,
            ..keep_gene_cols
        ]
    )

} else {

    candidate_genes <-
        data.table()
}


############################################################
## 17. FUNCTIONAL-ANNOTATION FLAGS
############################################################

if (nrow(candidate_genes) > 0) {

    candidate_genes[
        ,
        Has_named_product :=
            !is.na(functional_product) &
            functional_product != "" &
            functional_product !=
                "hypothetical protein"
    ]

    candidate_genes[
        ,
        Has_any_functional_annotation :=
            Has_named_product |
            !is.na(GO) |
            !is.na(FANTASIA_GO) |
            !is.na(InterPro) |
            !is.na(Pfam) |
            !is.na(PANTHER) |
            !is.na(Gene3D)
    ]
}


############################################################
## 18. PEAK-ANNOTATION TABLE
############################################################

peak_cols <- c(
    "Analysis",
    "Region_ID",
    "CHROM",
    "Start",
    "End",
    "Window_size_bp",
    "N_SNPs",
    "Statistic_peak_position",
    "Statistic_peak_value",
    "Peak_SNP",
    "Peak_context",
    "Peak_nearest_gene",
    "Peak_distance_to_gene",
    "Peak_in_gene",
    "Peak_in_exon",
    "Peak_in_CDS",
    "Peak_in_intron",
    "N_overlapping_genes",
    "Overlapping_gene_IDs",
    "Functional_products"
)

peak_cols <- peak_cols[
    peak_cols %in%
    names(regions)
]

peak_table <- regions[
    ,
    ..peak_cols
]


############################################################
## 19. SUMMARY BY ANALYSIS
############################################################

summary_regions <- regions[
    ,
    .(
        N_regions = .N,

        N_scaffolds =
            uniqueN(CHROM),

        Total_span_bp =
            sum(
                Window_size_bp,
                na.rm = TRUE
            ),

        Total_span_kb =
            sum(
                Window_size_bp,
                na.rm = TRUE
            ) / 1000,

        Min_region_bp =
            min(
                Window_size_bp,
                na.rm = TRUE
            ),

        Max_region_bp =
            max(
                Window_size_bp,
                na.rm = TRUE
            ),

        N_regions_overlapping_gene =
            sum(
                N_overlapping_genes > 0,
                na.rm = TRUE
            ),

        N_regions_gene_within_5kb =
            sum(
                Window_has_gene_within_5kb,
                na.rm = TRUE
            ),

        N_peaks_coding =
            sum(
                Peak_context == "coding",
                na.rm = TRUE
            ),

        N_peaks_intronic =
            sum(
                Peak_context == "intronic",
                na.rm = TRUE
            ),

        N_peaks_exonic_non_coding =
            sum(
                Peak_context ==
                    "exonic_non_coding",
                na.rm = TRUE
            ),

        N_peaks_gene_proximal_5kb =
            sum(
                Peak_context ==
                    "gene_proximal_5kb",
                na.rm = TRUE
            ),

        N_peaks_distal_intergenic =
            sum(
                Peak_context ==
                    "distal_intergenic",
                na.rm = TRUE
            )
    ),
    by = Analysis
]


############################################################
## 20. WRITE FINAL OUTPUTS
############################################################

setorder(
    regions,
    Analysis,
    CHROM,
    Start
)


fwrite(
    regions,
    file.path(
        OUT_DIR,
        "XtX_C2_candidate_regions_annotation.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


fwrite(
    region_gene,
    file.path(
        OUT_DIR,
        "XtX_C2_region_gene_overlaps.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


fwrite(
    candidate_genes,
    file.path(
        OUT_DIR,
        "XtX_C2_candidate_genes.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


fwrite(
    peak_table,
    file.path(
        OUT_DIR,
        "XtX_C2_peak_annotation.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


fwrite(
    summary_regions,
    file.path(
        OUT_DIR,
        "XtX_C2_annotation_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


############################################################
## 21. PRINT SUMMARY
############################################################

cat(
    "\n====================================================\n"
)

cat(
    "XTX / C2 CANDIDATE ANNOTATION COMPLETE\n"
)

cat(
    "====================================================\n\n"
)

print(
    summary_regions
)

cat(
    "\nOutput directory:\n",
    OUT_DIR,
    "\n"
)

cat(
    "\nNOTE:\n"
)

cat(
    paste0(
        "Only final all-SNP XtX and C2 candidate regions ",
        "were annotated in this script. GEA candidate ",
        "annotation is handled separately.\n"
    )
)

cat(
    "\nDONE\n"
)
