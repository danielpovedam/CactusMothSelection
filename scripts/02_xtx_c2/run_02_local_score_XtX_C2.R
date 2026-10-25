#!/usr/bin/env Rscript

############################################################
## run_02_local_score_XtX_C2.R
##
## Cactoblastis cactorum
## Final local-score analysis
##
## Dataset:
##   all-SNP dataset only
##
## Analyses:
##   XtX     -> representative converged seed 12345
##   C2 Host -> representative converged seed 12345
##
## Scaffold eligibility:
##   scaffold length > 1 Mb
##   AND BayPass min.nsnp = 20
##
## Local-score parameters:
##   min.maf = 0.05
##   xi      = 2
##   alpha   = 0.01
##
## GEA is NOT analysed using local score.
############################################################


suppressPackageStartupMessages({
    library(data.table)
})


############################################################
## 0. SETTINGS
############################################################

BASE <- "/home/dpoveda/micetral/Population_snps/baypass"

POPBASE <- "/home/dpoveda/micetral/Population_snps"

UTILS <- paste0(
    "/home/dpoveda/micetral/selection/baypass/",
    "baypass_public/utils/baypass_utils.R"
)

FAI_FILE <- file.path(
    POPBASE,
    "Ccactorum.fcs_gx_2026_FINAL.fa.fai"
)

INBASE <- file.path(
    BASE,
    "selection_scan",
    "prepared_statistics",
    "all_snps"
)

OUTBASE <- file.path(
    BASE,
    "selection_scan",
    "local_score_gt1Mb",
    "all_snps"
)

REP_SEED <- "12345"


############################################################
## 1. LOCAL-SCORE PARAMETERS
############################################################

MIN_MAF <- 0.05

MIN_NSNP <- 20

MIN_SCAFFOLD_LENGTH_BP <- 1000000L

XI_XTX <- 2
XI_C2  <- 2

ALPHA_XTX <- 0.01
ALPHA_C2  <- 0.01


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
            "Missing or empty input file:\n",
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


write_if_present <- function(
    x,
    file
) {

    if (
        !is.null(x) &&
        nrow(x) > 0
    ) {

        fwrite(
            as.data.table(x),
            file,
            sep = "\t",
            quote = FALSE,
            na = ""
        )

    } else {

        fwrite(
            data.table(),
            file,
            sep = "\t",
            quote = FALSE,
            na = ""
        )
    }
}


prepare_positions <- function(dt) {

    dt <- copy(dt)

    dt[
        ,
        CHROM :=
            as.character(CHROM)
    ]

    dt[
        ,
        POS :=
            as.integer(POS)
    ]

    dt[
        ,
        Scaffold_norm :=
            normalize_scaffold(
                CHROM
            )
    ]

    setorder(
        dt,
        CHROM,
        POS
    )

    dt
}


############################################################
## 3. INPUT FILES
############################################################

XTX_FILE <- file.path(
    INBASE,
    paste0(
        "XtX_seed",
        REP_SEED,
        ".tsv"
    )
)

C2_FILE <- file.path(
    INBASE,
    paste0(
        "C2_Host_seed",
        REP_SEED,
        ".tsv"
    )
)


for (
    f in c(
        UTILS,
        FAI_FILE,
        XTX_FILE,
        C2_FILE
    )
) {

    check_file(f)
}


############################################################
## 4. LOAD BAYPASS UTILITIES
############################################################

source(
    UTILS
)

if (
    !exists(
        "compute.local.scores"
    )
) {

    stop(
        paste0(
            "compute.local.scores() was not found ",
            "after sourcing BayPass utilities."
        )
    )
}


############################################################
## 5. OUTPUT DIRECTORIES
############################################################

dir.create(
    OUTBASE,
    recursive = TRUE,
    showWarnings = FALSE
)

XTX_OUTDIR <- file.path(
    OUTBASE,
    "XtX"
)

C2_OUTDIR <- file.path(
    OUTBASE,
    "C2_Host"
)

dir.create(
    XTX_OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    C2_OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)


############################################################
## 6. READ TRUE SCAFFOLD LENGTHS
############################################################

fai <- fread(
    FAI_FILE,
    header = FALSE,
    showProgress = FALSE
)

if (
    ncol(fai) < 2L
) {

    stop(
        "FAI file has fewer than two columns."
    )
}


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
            as.character(
                FASTA_scaffold
            ),

        Scaffold_length_bp =
            as.integer(
                Scaffold_length_bp
            )
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
        MIN_SCAFFOLD_LENGTH_BP
]


############################################################
## Check normalized scaffold IDs
############################################################

dup_fai <- fai[
    ,
    .N,
    by = Scaffold_norm
][
    N > 1
]


if (
    nrow(dup_fai) > 0L
) {

    stop(
        "Duplicated normalized scaffold names detected in FAI."
    )
}


eligible_scaffolds <- fai[
    Eligible_gt1Mb == TRUE,
    Scaffold_norm
]


cat(
    "\nTotal FASTA scaffolds: ",
    nrow(fai),
    "\n",
    sep = ""
)

cat(
    "Scaffolds >1 Mb: ",
    length(
        eligible_scaffolds
    ),
    "\n",
    sep = ""
)


############################################################
## 7. LOCAL-SCORE FUNCTION
############################################################

run_local_score_gt1Mb <- function(
    dt,
    analysis,
    xi,
    alpha,
    out_scores,
    out_windows,
    out_input,
    out_scaffold_qc
) {

    dt <- prepare_positions(
        dt
    )


    required <- c(
        "CHROM",
        "POS",
        "M_P",
        "log10_p"
    )


    missing_cols <- setdiff(
        required,
        names(dt)
    )


    if (
        length(missing_cols) > 0L
    ) {

        stop(
            analysis,
            " input missing required columns: ",
            paste(
                missing_cols,
                collapse = ", "
            )
        )
    }


    dt[
        ,
        M_P :=
            as.numeric(M_P)
    ]

    dt[
        ,
        log10_p :=
            as.numeric(log10_p)
    ]


    dt <- dt[
        !is.na(CHROM) &
        !is.na(POS) &
        !is.na(M_P) &
        !is.na(log10_p)
    ]


    ########################################################
    ## Attach true scaffold lengths
    ########################################################

    dt <- merge(

        dt,

        fai[
            ,
            .(
                Scaffold_norm,
                FASTA_scaffold,
                Scaffold_length_bp,
                Scaffold_length_Mb,
                Eligible_gt1Mb
            )
        ],

        by = "Scaffold_norm",

        all.x = TRUE,

        sort = FALSE
    )


    ########################################################
    ## Validate scaffold matching
    ########################################################

    unmapped <- dt[
        is.na(
            Scaffold_length_bp
        )
    ]


    if (
        nrow(unmapped) > 0L
    ) {

        stop(
            analysis,
            ": ",
            nrow(unmapped),
            " SNPs could not be mapped to the FASTA index."
        )
    }


    original_n <-
        nrow(dt)

    original_scaffolds <-
        uniqueN(
            dt$Scaffold_norm
        )


    ########################################################
    ## Keep ONLY scaffolds >1 Mb
    ########################################################

    dt_gt1Mb <- dt[
        Eligible_gt1Mb == TRUE
    ]


    gt1Mb_n <-
        nrow(
            dt_gt1Mb
        )

    gt1Mb_scaffolds <-
        uniqueN(
            dt_gt1Mb$Scaffold_norm
        )


    ########################################################
    ## Scaffold QC
    ##
    ## BayPass compute.local.scores() additionally applies
    ## chr.names > min.nsnp.
    ##
    ## With min.nsnp = 20, scaffolds must contain >20 SNPs.
    ########################################################

    scaffold_qc <- dt_gt1Mb[
        ,
        .(
            N_input_SNPs =
                .N,

            Min_position =
                min(
                    POS,
                    na.rm = TRUE
                ),

            Max_position =
                max(
                    POS,
                    na.rm = TRUE
                ),

            Scaffold_length_bp =
                unique(
                    Scaffold_length_bp
                ),

            Scaffold_length_Mb =
                unique(
                    Scaffold_length_Mb
                )
        ),
        by = .(
            CHROM,
            Scaffold_norm
        )
    ]


    scaffold_qc[
        ,
        Eligible_gt20SNP :=
            N_input_SNPs >
            MIN_NSNP
    ]


    setorder(
        scaffold_qc,
        -Scaffold_length_bp
    )


    fwrite(
        scaffold_qc,
        out_scaffold_qc,
        sep = "\t",
        quote = FALSE,
        na = ""
    )


    ########################################################
    ## Save exact input used
    ########################################################

    fwrite(
        dt_gt1Mb,
        out_input,
        sep = "\t",
        quote = FALSE,
        na = ""
    )


    ########################################################
    ## Report filtering
    ########################################################

    cat(
        "\n--------------------------------------------\n"
    )

    cat(
        "Analysis: ",
        analysis,
        "\n",
        sep = ""
    )

    cat(
        "Original SNPs: ",
        original_n,
        "\n",
        sep = ""
    )

    cat(
        "Original scaffolds: ",
        original_scaffolds,
        "\n",
        sep = ""
    )

    cat(
        "SNPs on >1 Mb scaffolds: ",
        gt1Mb_n,
        "\n",
        sep = ""
    )

    cat(
        "Scaffolds >1 Mb represented: ",
        gt1Mb_scaffolds,
        "\n",
        sep = ""
    )

    cat(
        "Scaffolds >1 Mb AND >20 SNPs: ",
        sum(
            scaffold_qc$
                Eligible_gt20SNP
        ),
        "\n",
        sep = ""
    )


    ########################################################
    ## Run BayPass local score
    ########################################################

    res <- compute.local.scores(

        snp.position = data.frame(
            chr =
                dt_gt1Mb$CHROM,

            pos =
                dt_gt1Mb$POS
        ),

        snp.pi =
            dt_gt1Mb$M_P,

        snp.pvalue =
            dt_gt1Mb$log10_p,

        xi =
            xi,

        pval.local.score.thres =
            alpha,

        min.maf =
            MIN_MAF,

        min.nsnp =
            MIN_NSNP,

        plot.pvalhist =
            FALSE,

        manplot =
            FALSE
    )


    ########################################################
    ## SNP-level local scores
    ########################################################

    scores <- as.data.table(
        res$res.local.scores
    )


    if (
        nrow(scores) > 0L
    ) {

        scores[
            ,
            Dataset :=
                "all_snps"
        ]

        scores[
            ,
            Analysis :=
                analysis
        ]

        setcolorder(
            scores,
            c(
                "Dataset",
                "Analysis",
                setdiff(
                    names(scores),
                    c(
                        "Dataset",
                        "Analysis"
                    )
                )
            )
        )
    }


    write_if_present(
        scores,
        out_scores
    )


    ########################################################
    ## Significant windows
    ########################################################

    windows <- as.data.table(
        res$significant.windows
    )


    if (
        nrow(windows) > 0L
    ) {

        windows[
            ,
            Dataset :=
                "all_snps"
        ]

        windows[
            ,
            Analysis :=
                analysis
        ]

        setcolorder(
            windows,
            c(
                "Dataset",
                "Analysis",
                setdiff(
                    names(windows),
                    c(
                        "Dataset",
                        "Analysis"
                    )
                )
            )
        )
    }


    write_if_present(
        windows,
        out_windows
    )


    cat(
        "Local-score SNPs returned: ",
        nrow(scores),
        "\n",
        sep = ""
    )

    cat(
        "Candidate windows returned: ",
        nrow(windows),
        "\n",
        sep = ""
    )


    invisible(
        list(
            scores = scores,
            windows = windows,
            scaffold_qc = scaffold_qc,
            N_original = original_n,
            N_gt1Mb = gt1Mb_n
        )
    )
}


############################################################
## 8. XtX LOCAL SCORE
############################################################

cat(
    "\n============================================\n"
)

cat(
    "RUNNING XtX LOCAL SCORE: >1 Mb SCAFFOLDS\n"
)

cat(
    "============================================\n"
)


xtx <- fread(
    XTX_FILE,
    showProgress = FALSE
)


xtx_result <- run_local_score_gt1Mb(

    dt =
        xtx,

    analysis =
        "XtX",

    xi =
        XI_XTX,

    alpha =
        ALPHA_XTX,

    out_scores =
        file.path(
            XTX_OUTDIR,
            "XtX_local_scores_gt1Mb.tsv"
        ),

    out_windows =
        file.path(
            XTX_OUTDIR,
            "XtX_candidate_windows_gt1Mb.tsv"
        ),

    out_input =
        file.path(
            XTX_OUTDIR,
            "XtX_input_gt1Mb.tsv"
        ),

    out_scaffold_qc =
        file.path(
            XTX_OUTDIR,
            "XtX_scaffold_QC_gt1Mb.tsv"
        )
)


rm(xtx)

gc(FALSE)


############################################################
## 9. C2 HOST LOCAL SCORE
############################################################

cat(
    "\n============================================\n"
)

cat(
    "RUNNING C2 LOCAL SCORE: >1 Mb SCAFFOLDS\n"
)

cat(
    "============================================\n"
)


c2 <- fread(
    C2_FILE,
    showProgress = FALSE
)


c2_result <- run_local_score_gt1Mb(

    dt =
        c2,

    analysis =
        "C2_Host",

    xi =
        XI_C2,

    alpha =
        ALPHA_C2,

    out_scores =
        file.path(
            C2_OUTDIR,
            "C2_Host_local_scores_gt1Mb.tsv"
        ),

    out_windows =
        file.path(
            C2_OUTDIR,
            "C2_Host_candidate_windows_gt1Mb.tsv"
        ),

    out_input =
        file.path(
            C2_OUTDIR,
            "C2_Host_input_gt1Mb.tsv"
        ),

    out_scaffold_qc =
        file.path(
            C2_OUTDIR,
            "C2_Host_scaffold_QC_gt1Mb.tsv"
        )
)


rm(c2)

gc(FALSE)


############################################################
## 10. RUN SUMMARY
############################################################

RUN_SUMMARY <- rbindlist(

    list(

        data.table(

            Dataset =
                "all_snps",

            Analysis =
                "XtX",

            Representative_seed =
                REP_SEED,

            N_original_SNPs =
                xtx_result$N_original,

            N_SNPs_gt1Mb =
                xtx_result$N_gt1Mb,

            Percent_SNPs_retained =
                100 *
                xtx_result$N_gt1Mb /
                xtx_result$N_original,

            N_scaffolds_gt1Mb =
                nrow(
                    xtx_result$
                        scaffold_qc
                ),

            N_scaffolds_gt1Mb_gt20SNP =
                sum(
                    xtx_result$
                        scaffold_qc$
                        Eligible_gt20SNP
                ),

            N_local_score_SNPs =
                nrow(
                    xtx_result$scores
                ),

            N_candidate_windows =
                nrow(
                    xtx_result$windows
                )
        ),


        data.table(

            Dataset =
                "all_snps",

            Analysis =
                "C2_Host",

            Representative_seed =
                REP_SEED,

            N_original_SNPs =
                c2_result$N_original,

            N_SNPs_gt1Mb =
                c2_result$N_gt1Mb,

            Percent_SNPs_retained =
                100 *
                c2_result$N_gt1Mb /
                c2_result$N_original,

            N_scaffolds_gt1Mb =
                nrow(
                    c2_result$
                        scaffold_qc
                ),

            N_scaffolds_gt1Mb_gt20SNP =
                sum(
                    c2_result$
                        scaffold_qc$
                        Eligible_gt20SNP
                ),

            N_local_score_SNPs =
                nrow(
                    c2_result$scores
                ),

            N_candidate_windows =
                nrow(
                    c2_result$windows
                )
        )
    ),

    use.names = TRUE,

    fill = TRUE
)


fwrite(
    RUN_SUMMARY,
    file.path(
        OUTBASE,
        "local_score_gt1Mb_run_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
)


############################################################
## 11. SETTINGS
############################################################

SETTINGS <- data.table(

    Parameter = c(
        "Dataset",
        "Representative_seed_XtX",
        "Representative_seed_C2",
        "Minimum_scaffold_length_bp",
        "Scaffold_length_rule",
        "min_maf",
        "min_nsnp",
        "BayPass_min_nsnp_rule",
        "xi_XtX",
        "xi_C2",
        "alpha_XtX",
        "alpha_C2"
    ),

    Value = c(
        "all_snps",
        REP_SEED,
        REP_SEED,
        MIN_SCAFFOLD_LENGTH_BP,
        ">1000000_bp",
        MIN_MAF,
        MIN_NSNP,
        ">20_SNPs_per_scaffold",
        XI_XTX,
        XI_C2,
        ALPHA_XTX,
        ALPHA_C2
    )
)


fwrite(
    SETTINGS,
    file.path(
        OUTBASE,
        "local_score_gt1Mb_settings.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
)


############################################################
## 12. FINAL SUMMARY
############################################################

cat(
    "\n============================================\n"
)

cat(
    "FINAL LOCAL-SCORE RUN SUMMARY\n"
)

cat(
    "============================================\n"
)

print(
    RUN_SUMMARY
)

cat(
    "\nSETTINGS\n"
)

print(
    SETTINGS
)

cat(
    "\nOutput directory:\n",
    OUTBASE,
    "\n"
)

cat(
    "\nNOTE:\n"
)

cat(
    paste0(
        "Local score was applied only to XtX and C2 ",
        "on scaffolds >1 Mb. GEA candidate detection ",
        "was based directly on median BF(dB) >20.\n"
    )
)

cat(
    "\nDONE\n"
)
