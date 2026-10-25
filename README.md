# CactusMothSelection
Scripts for paper: Contemporary ecological heterogeneity shapes candidate adaptive genomic variation across the native range of an invasive herbivore

# CactusMothSelection

Code associated with the population- and landscape-genomic analyses of the native-range *Cactoblastis cactorum* dataset presented in:

**Contemporary ecological heterogeneity shapes candidate adaptive genomic variation across the native range of an invasive herbivore**

This repository contains the final scripts used for BayPass genome scans, local-score candidate detection, pseudo-observed dataset (POD) calibration, candidate-SNP overlap analyses, host-use sensitivity analyses, candidate-gene annotation, variant-effect prediction, and Gene Ontology enrichment.

## Study overview

The analyses use ddRADseq data from 136 *Cactoblastis cactorum* individuals sampled from 28 populations across Argentina.

The final all-SNP dataset contains 132,283 SNPs and was used for the genome-wide selection and genotype–environment association analyses.

Three complementary BayPass analyses were performed:

- **XtX:** genome-wide population differentiation after accounting for the population covariance structure.
- **C2:** differentiation associated with predominant host use, contrasting populations sampled predominantly from native *Opuntia* with populations sampled predominantly from introduced *O. ficus-indica*.
- **GEA:** genotype–environment association analyses for TemperaturePC1, PrecipitationPC1, aridity index, elevation, latitude, and longitude.

The repository contains only the final analyses retained in the manuscript. Exploratory analyses and obsolete intermediate versions are not included.

---

# Repository structure

```text
CactusMothSelection/
├── README.md
└── scripts/
    ├── 01_baypass_core/
    │   ├── run_00_baypass_all_snps.sh
    │   └── run_01_prepare_baypass_results.R
    │
    ├── 02_xtx_c2/
    │   ├── run_02_local_score_XtX_C2.R
    │   └── run_03_annotate_XtX_C2_candidates.R
    │
    ├── 03_gea/
    │   └── run_04_GEA_BF20_candidates.R
    │
    ├── 04_pod/
    │   ├── run_05a_simulate_neutral_POD.R
    │   ├── run_05b_run_POD_BayPass.sh
    │   └── run_05c_calculate_POD_thresholds.R
    │
    ├── 05_candidate_overlap/
    │   └── run_06_candidate_SNP_overlap_hypergeometric.R
    │
    ├── 06_c2_geographic_sensitivity/
    │   ├── run_07a_prepare_geographic_C2_subset.R
    │   ├── run_07b_run_geographic_C2_sensitivity.sh
    │   └── run_07c_local_score_C2_geographic_sensitivity.R
    │
    ├── 07_c2_peak_followup/
    │   ├── run_08a_C2_peakSNP_host_frequency_profiles.R
    │   ├── run_08b_C2_peakSNP_lineage_consistency.R
    │   └── run_08c_C2_peakSNP_summary_figure_tables.R
    │
    └── 08_genes_analyses/
        ├── 01_candidate_gene_annotation/
        ├── 02_BLAST_InterPro/
        ├── 03_VEP/
        └── 04_GO_enrichment/
```

---

# Analysis workflow

## 1. BayPass core analysis

### `run_00_baypass_all_snps.sh`

Runs BayPass v3.1 on the final all-SNP population allele-count matrix.

Three independent MCMC chains were run with seeds:

```text
12345
67890
13579
```

The analyses included:

- estimation of the population covariance matrix Ω;
- XtX;
- the C2 host-use contrast;
- genotype–environment associations.

The final BayPass settings were:

```text
npilot      = 20
pilotlength = 500
burnin      = 10000
nval        = 2500
thin        = 20
```

The C2 contrast coded populations predominantly sampled from native *Opuntia* as −1 and those predominantly sampled from *O. ficus-indica* as +1.

### `run_01_prepare_baypass_results.R`

Processes BayPass outputs and evaluates convergence among the three chains.

Convergence diagnostics include:

- Förstner–Moonen distances among Ω matrices;
- posterior Beta-prior parameter stability;
- pairwise correlations of XtX estimates;
- consistency of C2 statistics;
- covariate-specific consistency of Bayes Factors.

For downstream XtX and C2 analyses, seed **12345** was used as the representative converged run.

For GEA analyses, the final statistic for each SNP–covariate combination is the **median BF(dB) across the three runs**.

---

# 2. XtX and C2 local-score analyses

### `run_02_local_score_XtX_C2.R`

Applies the BayPass local-score approach separately to XtX and C2.

Local scores are calculated only on scaffolds **>1 Mb**.

Final parameters:

```text
MAF > 0.05
xi = 2
alpha = 0.01
min.nsnp = 20
```

Local-score candidate regions are the primary region-level definition of XtX and C2 candidates.

GEA statistics are not analysed using the local-score procedure.

### `run_03_annotate_XtX_C2_candidates.R`

Annotates the final XtX and C2 candidate regions and peak SNPs against the *C. cactorum* reference annotation.

Peak locations are classified as:

- coding;
- exonic non-coding;
- intronic;
- gene-proximal, when located within 5 kb of an annotated gene;
- distal intergenic.

When a peak does not overlap a gene, the nearest annotated gene within 5 kb is reported.

---

# 3. Genotype–environment associations

### `run_04_GEA_BF20_candidates.R`

Processes the final GEA results for six environmental or spatial covariates:

```text
TemperaturePC1
PrecipitationPC1
AridityIndex
Elevation
Latitude
Longitude
```

For each SNP–covariate combination, the median BF across the three converged BayPass runs is used.

Strong associations are defined as:

```text
median BF > 20 dB
```

Significant SNPs belonging to the same covariate, located on the same scaffold, and separated by no more than 1 kb are grouped into candidate loci.

The SNP with the highest median BF within each locus is designated as the locus peak.

---

# 4. Neutral POD calibration

The POD analysis provides an independent SNP-level neutral calibration for XtX and C2.

It complements, but does not replace, the local-score candidate-region analysis.

### `run_05a_simulate_neutral_POD.R`

Generates a neutral pseudo-observed dataset of:

```text
100,000 SNPs
```

using the BayPass simulation framework and parameters estimated from the empirical data, including:

- the empirical Ω matrix;
- empirical Beta-prior parameters;
- the empirical population sample-size distribution.

No association with host use is simulated.

### `run_05b_run_POD_BayPass.sh`

Runs BayPass on the neutral POD to obtain null distributions for:

- XtXst;
- C2.

The host contrast is supplied only so that BayPass calculates the null C2 statistic corresponding to the empirical contrast.

### `run_05c_calculate_POD_thresholds.R`

Calculates the observed pooled MAF of the simulated SNPs and restricts the primary neutral calibration to simulated variants satisfying:

```text
MAF >= 0.05
```

The primary threshold is the **99.9th percentile** of the MAF-filtered neutral distribution.

Final POD99.9 thresholds:

```text
XtXst = 57.39
C2    = 12.08
```

POD99.9 outliers located within local-score candidate regions provide complementary SNP-level support for XtX and C2 candidates.

---

# 5. Candidate-SNP convergence

### `run_06_candidate_SNP_overlap_hypergeometric.R`

Tests whether candidate SNPs are shared among analyses more frequently than expected by chance using one-sided hypergeometric tests.

XtX and C2 SNP sets are defined as SNPs satisfying **both** criteria:

```text
POD99.9 outlier
AND
located within a final local-score candidate region
```

GEA candidate sets consist of SNPs with:

```text
median BF > 20 dB
```

Tests include:

- XtX versus C2;
- XtX versus each of the six GEA covariates;
- C2 versus each of the six GEA covariates.

Benjamini–Hochberg correction is applied across the 12 covariate-specific GEA comparisons.

The background universe corresponds to the unique genomic SNP positions represented in the final BayPass dataset.

---

# 6. Geographic sensitivity of the C2 host-use signal

Because host categories are geographically structured across the native range, the C2 analysis was repeated using a geographically restricted subset.

### `run_07a_prepare_geographic_C2_subset.R`

Constructs the restricted dataset using populations from provinces in which both predominant host categories are represented.

The final restricted analysis contains:

```text
16 populations
8 native-Opuntia populations
8 O. ficus-indica populations
```

### `run_07b_run_geographic_C2_sensitivity.sh`

Re-estimates Ω and reruns the C2 contrast independently for the geographically restricted dataset.

### `run_07c_local_score_C2_geographic_sensitivity.R`

Applies the same final local-score framework used in the primary analysis, including restriction to scaffolds >1 Mb.

Concordance between the primary and geographically restricted analyses is evaluated using:

- genome-wide C2 correspondence;
- ranking of the primary C2 peak SNPs;
- physical overlap among candidate regions;
- shared peak-associated genes.

---

# 7. C2 peak-SNP follow-up

These analyses were added to examine the population allele-frequency patterns underlying the final host-associated candidate regions.

### `run_08a_C2_peakSNP_host_frequency_profiles.R`

Extracts the BayPass posterior population allele frequency (`M_P`) for each of the **17 final C2 peak SNPs** across the 28 populations.

The analysis summarizes allele frequencies according to predominant host category:

```text
14 native-Opuntia populations
14 O. ficus-indica populations
```

Raw `M_P` values are retained without reorientation to `p` or `1-p`.

### `run_08b_C2_peakSNP_lineage_consistency.R`

Evaluates whether host-associated allele-frequency differences at the final C2 peaks show the same direction within the three broad geographic lineages recovered by population-structure analyses.

This analysis is used as a descriptive robustness check against broad geographic confounding.

### `run_08c_C2_peakSNP_summary_figure_tables.R`

Generates the final supplementary summaries and graphical outputs for the C2 peak-SNP follow-up.

---

# 8. Candidate genes and functional annotation

## 8.1 Candidate-gene annotation

Scripts in:

```text
scripts/08_genes_analyses/01_candidate_gene_annotation/
```

build the final XtX, C2, and GEA candidate-gene tables and assign peak SNPs to overlapping or nearby genes.

Candidate-gene convergence is evaluated independently of physical overlap between genomic candidate regions.

---

## 8.2 BLAST and InterPro refinement

Scripts in:

```text
scripts/08_genes_analyses/02_BLAST_InterPro/
```

were used to refine the annotation of recurrent candidate genes that remained uncharacterized.

Predicted proteins were queried using BLASTP against NCBI protein databases, and selected proteins were additionally evaluated using InterProScan.

---

## 8.3 Variant Effect Predictor

Scripts in:

```text
scripts/08_genes_analyses/03_VEP/
```

prepare and annotate all unique candidate SNPs identified across the final XtX, C2, and GEA analyses.

The final union contains:

```text
567 unique candidate SNP positions
```

Variants were annotated using Ensembl Variant Effect Predictor (VEP v111).

When multiple transcript consequences were assigned to the same SNP, one primary consequence was retained using the hierarchy:

```text
protein-altering
>
UTR
>
synonymous
>
non-coding exon/transcript
>
intronic or splice-region
>
upstream/downstream
>
intergenic
```

---

## 8.4 Gene Ontology enrichment

Scripts in:

```text
scripts/08_genes_analyses/04_GO_enrichment/
```

perform GO enrichment analyses using `topGO`.

Analyses include:

- XtX;
- C2;
- pooled GEA;
- individual GEA covariates;
- combined candidate-gene set.

Biological Process, Molecular Function, and Cellular Component ontologies are analysed separately using Fisher's exact test with the `weight01` algorithm.

The background consists of annotated genes overlapping or located within 5 kb of SNPs screened in the final all-SNP BayPass dataset.

Benjamini–Hochberg correction is applied within each analysis and ontology, with:

```text
FDR < 0.05
```

used as the significance threshold.

---

# Data availability

Raw sequencing data are available through NCBI under BioProject:

```text
PRJNA666743
```

Processed genomic datasets and BayPass input files associated with the landscape-genomic analyses are available from Figshare:

**Genomic VCF and BayPass input datasets for native-range landscape genomic analyses of *Cactoblastis cactorum***

DOI:

```text
10.6084/m9.figshare.33340938
```

Scripts for the upstream ddRADseq processing workflow are available separately at:

```text
https://github.com/Parany24/Cactus_moth_RADseq
```

This repository focuses on the downstream population-genomic and candidate-locus analyses associated with the manuscript.

---

# Software

Principal software and R packages used in these analyses include:

- BayPass v3.1
- R
- data.table
- GenomicRanges
- IRanges
- topGO
- Ensembl Variant Effect Predictor v111
- NCBI BLAST+ v2.15.0
- InterProScan v5.73-104.0

Additional package dependencies are specified within individual scripts where applicable.

---

# Reproducibility notes

The scripts were originally developed and executed on an HPC environment and therefore contain paths corresponding to the original analysis directory structure.

Users reproducing the workflow should update project, reference-genome, annotation, and output paths as appropriate for their environment.

The repository preserves the analytical decisions used for the final manuscript, including:

- three independent BayPass runs;
- representative seed 12345 for XtX and C2;
- median BF across runs for GEA;
- local-score analysis restricted to scaffolds >1 Mb;
- `xi = 2` and `alpha = 0.01` for XtX and C2;
- MAF >0.05 filtering;
- BF >20 dB for GEA candidates;
- ≤1-kb grouping of significant GEA SNPs;
- POD99.9 neutral calibration;
- geographically restricted C2 sensitivity analysis;
- C2 peak-SNP allele-frequency follow-up;
- final candidate-SNP VEP annotation;
- GO enrichment using the BayPass-screened genomic background.

---

# Citation

If using these scripts or datasets, please cite the associated manuscript and data repository.

Full manuscript citation will be added following publication.

---

# Contact

Daniel Poveda-Martínez

For questions concerning the analyses or repository, please open a GitHub issue.
