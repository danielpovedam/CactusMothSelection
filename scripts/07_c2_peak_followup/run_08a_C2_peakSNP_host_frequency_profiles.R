library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# PATHS
# ============================================================

BASE <- "~/work/micetral/Population_snps/baypass"

PIJ <- file.path(
  BASE,
  "runs_10x/all_snps/seed_13579/",
  "BP_all_snps_seed13579_summary_pij.out"
)

PEAK_FILE <- file.path(
  BASE,
  "sensitivity_C2_geographic/",
  "comparison_primary_vs_restricted_gt1Mb/",
  "Table_C2_final17_peak_retention_and_rank.tsv"
)

POP_ORDER <- "~/work/micetral/Population_snps/baypass_preparation/input/baypass_population_order.txt"

OUT <- file.path(
  BASE,
  "selection_scan/C2_peakSNP_BayPass_pij_host_profiles"
)

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# POPULATION METADATA
# ============================================================

pop_meta <- tribble(
  ~PopulationID, ~Locality, ~Province, ~Host, ~Host_C2, ~Lineage,
  "BAM","Mayor Buratovich","Buenos Aires","O. ficus-indica",1,"Southern",
  "BAP","Tornquist","Buenos Aires","O. penicilligera",-1,"Southern",
  "BAS","Sierra la Ventana","Buenos Aires","O. ficus-indica",1,"Southern",
  "CAQ","Quirós","Catamarca","O. ficus-indica",1,"Central-northern-western",
  "CAS","Chumbicha","Catamarca","O. ficus-indica",1,"Central-northern-western",
  "CBA","El Aranado","Córdoba","O. megapotamica",-1,"Central-northern-western",
  "CBF","El Fortin","Córdoba","O. rioplatense",-1,"Eastern",
  "CBR","Rincón","Córdoba","O. ficus-indica",1,"Eastern",
  "CBS","Centeno","Santa Fe","O. rioplatense",-1,"Eastern",
  "CBV","Villa del Totoral","Córdoba","O. ficus-indica",1,"Central-northern-western",
  "CHP","Puerto Tirol","Chaco","O. ficus-indica",1,"Central-northern-western",
  "ERB","Bovril","Entre Rios","O. elata",-1,"Eastern",
  "ERC","Chajarí","Entre Rios","O. rioplatense",-1,"Eastern",
  "ERD","Federación","Entre Rios","O. ficus-indica",1,"Eastern",
  "ERF","Federal","Entre Rios","O. ficus-indica",1,"Eastern",
  "ERL","Los conquistadores","Entre Rios","O. rioplatense",-1,"Eastern",
  "ERN","Nogoyá","Entre Rios","O. megapotamica",-1,"Eastern",
  "FOG","Gran Guardia","Formosa","O. anacantha",-1,"Central-northern-western",
  "JUJ","El Carmen","Jujuy","O. ficus-indica",1,"Central-northern-western",
  "JUS","Lumbreras","Salta","O. ficus-indica",1,"Central-northern-western",
  "LPA","Algarrobo del Aguila","La Pampa","O. megapotamica",-1,"Southern",
  "LPS","Santa Rosa","La Pampa","O. bonaerensis",-1,"Southern",
  "LPV","La Victoria","La Pampa","O. megapotamica",-1,"Southern",
  "RNR","Choele Choel","Río Negro","O. megapotamica",-1,"Southern",
  "SEB","Brea Puñuna","Santiago del Estero","O. ficus-indica",1,"Central-northern-western",
  "SEL","Lavalle","Santiago del Estero","O. quimilo",-1,"Central-northern-western",
  "SEP","San Lorenzo","Santiago del Estero","O. ficus-indica",1,"Central-northern-western",
  "TUN","La Cocha","Tucumán","O. ficus-indica",1,"Central-northern-western"
) %>%
  mutate(
    Host_category = ifelse(
      Host_C2 == 1,
      "O. ficus-indica",
      "Native Opuntia"
    ),
    Host_category = factor(
      Host_category,
      levels = c("Native Opuntia", "O. ficus-indica")
    ),
    Lineage = factor(
      Lineage,
      levels = c(
        "Southern",
        "Eastern",
        "Central-northern-western"
      )
    )
  )

# ============================================================
# READ POPULATION ORDER
# ============================================================

pop_order <- fread(
  POP_ORDER,
  header = FALSE
)[[1]]

stopifnot(length(pop_order) == 28)

pop_map <- tibble(
  POP = seq_along(pop_order),
  PopulationID = pop_order
)

# ============================================================
# READ FINAL 17 C2 PEAKS
# ============================================================

peaks <- fread(PEAK_FILE)

stopifnot(nrow(peaks) == 17)

peaks <- peaks %>%
  mutate(
    Candidate_label = case_when(
      Primary_Region_ID == "all_snps_C2_R016" ~ "Masc",
      !is.na(Primary_Gene_symbol) &
        Primary_Gene_symbol != "" ~ Primary_Gene_symbol,
      TRUE ~ Scaffold_norm
    )
  )

# ============================================================
# READ BAYPASS POSTERIOR FREQUENCIES
# ============================================================

pij <- fread(PIJ)

cat("Rows in summary_pij:", nrow(pij), "\n")
cat("Expected:", 132283 * 28, "\n")

# retain only the 17 peak markers
peak_pij <- pij %>%
  filter(MRK %in% peaks$BayPass_row) %>%
  left_join(
    pop_map,
    by = "POP"
  ) %>%
  left_join(
    peaks %>%
      select(
        BayPass_row,
        Primary_Region_ID,
        Scaffold_norm,
        Primary_Statistic_peak_position,
        Primary_Statistic_peak_value,
        Candidate_label,
        Restricted_C2,
        Restricted_rank,
        Restricted_percentile,
        Restricted_top_1pct,
        `Restricted_top_0.1pct`,
        Primary_region_recovered_by_physical_overlap
      ),
    by = c("MRK" = "BayPass_row")
  ) %>%
  left_join(
    pop_meta,
    by = "PopulationID"
  )

stopifnot(nrow(peak_pij) == 17 * 28)
stopifnot(n_distinct(peak_pij$MRK) == 17)
stopifnot(n_distinct(peak_pij$PopulationID) == 28)

# ============================================================
# ORIENT FREQUENCIES CONSISTENTLY FOR VISUALIZATION
# ============================================================
#
# BayPass allele orientation is arbitrary across loci.
# For descriptive comparison, orient each SNP so that the allele
# plotted has the higher posterior mean frequency in
# O. ficus-indica-associated populations.
#
# Raw M_P values are preserved.
# ============================================================

orientation <- peak_pij %>%
  group_by(
    MRK,
    Primary_Region_ID,
    Candidate_label
  ) %>%
  summarise(
    Mean_MP_OFI =
      mean(M_P[Host_C2 == 1], na.rm = TRUE),

    Mean_MP_native =
      mean(M_P[Host_C2 == -1], na.rm = TRUE),

    Delta_MP =
      Mean_MP_OFI - Mean_MP_native,

    Flip =
      Delta_MP < 0,

    .groups = "drop"
  )

peak_pij <- peak_pij %>%
  left_join(
    orientation,
    by = c(
      "MRK",
      "Primary_Region_ID",
      "Candidate_label"
    )
  ) %>%
  mutate(
    M_P_oriented =
      ifelse(
        Flip,
        1 - M_P,
        M_P
      )
  )

# ============================================================
# EXPORT RAW + ORIENTED POPULATION FREQUENCIES
# ============================================================

fwrite(
  peak_pij,
  file.path(
    OUT,
    "Table_C2_17_peakSNP_BayPass_posterior_population_frequencies.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

# ============================================================
# BINARY HOST-CATEGORY SUMMARY
# ============================================================

binary_summary <- peak_pij %>%
  group_by(
    MRK,
    Primary_Region_ID,
    Candidate_label,
    Host_category
  ) %>%
  summarise(
    N_populations = n(),
    Mean_MP = mean(M_P_oriented, na.rm = TRUE),
    Median_MP = median(M_P_oriented, na.rm = TRUE),
    SD_MP = sd(M_P_oriented, na.rm = TRUE),
    Min_MP = min(M_P_oriented, na.rm = TRUE),
    Max_MP = max(M_P_oriented, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  binary_summary,
  file.path(
    OUT,
    "Table_C2_17_peakSNP_binary_host_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

# ============================================================
# HOST-SPECIES SUMMARY
# ============================================================

host_species_summary <- peak_pij %>%
  group_by(
    MRK,
    Primary_Region_ID,
    Candidate_label,
    Host
  ) %>%
  summarise(
    N_populations = n(),
    Mean_MP = mean(M_P_oriented, na.rm = TRUE),
    Median_MP = median(M_P_oriented, na.rm = TRUE),
    Min_MP = min(M_P_oriented, na.rm = TRUE),
    Max_MP = max(M_P_oriented, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  host_species_summary,
  file.path(
    OUT,
    "Table_C2_17_peakSNP_host_species_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

# ============================================================
# LINEAGE x HOST SUMMARY
# ============================================================

lineage_host_summary <- peak_pij %>%
  group_by(
    MRK,
    Primary_Region_ID,
    Candidate_label,
    Lineage,
    Host_category
  ) %>%
  summarise(
    N_populations = n(),
    Mean_MP = mean(M_P_oriented, na.rm = TRUE),
    Median_MP = median(M_P_oriented, na.rm = TRUE),
    Min_MP = min(M_P_oriented, na.rm = TRUE),
    Max_MP = max(M_P_oriented, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  lineage_host_summary,
  file.path(
    OUT,
    "Table_C2_17_peakSNP_lineage_host_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

# ============================================================
# FIGURE 1:
# NATIVE vs O. FICUS-INDICA FOR EACH PEAK
# ============================================================

candidate_order <- peaks %>%
  arrange(BayPass_row) %>%
  pull(Candidate_label)

peak_pij <- peak_pij %>%
  mutate(
    Candidate_label = factor(
      Candidate_label,
      levels = candidate_order
    )
  )

p_binary <- ggplot(
  peak_pij,
  aes(
    x = Host_category,
    y = M_P_oriented
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    aes(shape = Lineage),
    width = 0.12,
    size = 1.7,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ Candidate_label,
    ncol = 4,
    scales = "fixed"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    x = NULL,
    y = "BayPass posterior allele frequency",
    shape = "Major lineage"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    OUT,
    "Figure_C2_17_peakSNP_native_vs_OFI.pdf"
  ),
  p_binary,
  width = 11,
  height = 10
)

ggsave(
  file.path(
    OUT,
    "Figure_C2_17_peakSNP_native_vs_OFI.png"
  ),
  p_binary,
  width = 11,
  height = 10,
  dpi = 400
)

# ============================================================
# FIGURE 2:
# ACTUAL HOST SPECIES
# ============================================================

p_species <- ggplot(
  peak_pij,
  aes(
    x = Host,
    y = M_P_oriented,
    shape = Lineage
  )
) +
  geom_point(
    size = 2,
    alpha = 0.85,
    position = position_jitter(
      width = 0.10,
      height = 0
    )
  ) +
  facet_wrap(
    ~ Candidate_label,
    ncol = 4
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    x = "Host species",
    y = "BayPass posterior allele frequency",
    shape = "Major lineage"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1
    ),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    OUT,
    "Figure_C2_17_peakSNP_host_species.pdf"
  ),
  p_species,
  width = 13,
  height = 11
)

# ============================================================
# FIGURE 3:
# SIMPLE POPULATION HEATMAP
# ============================================================

population_order <- pop_meta %>%
  arrange(
    Lineage,
    Host_category,
    Host,
    PopulationID
  ) %>%
  pull(PopulationID)

peak_pij <- peak_pij %>%
  mutate(
    PopulationID = factor(
      PopulationID,
      levels = population_order
    )
  )

p_heat <- ggplot(
  peak_pij,
  aes(
    x = PopulationID,
    y = Candidate_label,
    fill = M_P_oriented
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.25
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Posterior\nallele frequency"
  ) +
  labs(
    x = "Population",
    y = "C2 candidate region"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "italic"
    ),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(
    OUT,
    "Figure_C2_17_peakSNP_population_heatmap.pdf"
  ),
  p_heat,
  width = 12,
  height = 7
)

# ============================================================
# QC
# ============================================================

cat("\n==============================\n")
cat("FINAL QC\n")
cat("==============================\n")

cat("Peak SNPs:", n_distinct(peak_pij$MRK), "\n")
cat("Populations:", n_distinct(peak_pij$PopulationID), "\n")
cat("Rows:", nrow(peak_pij), "\n")

cat("\nNative vs OFI populations:\n")
print(table(pop_meta$Host_category))

cat("\nActual host species:\n")
print(table(pop_meta$Host))

cat("\nLineage x host category:\n")
print(table(pop_meta$Lineage, pop_meta$Host_category))

cat("\nSensitivity support:\n")
print(
  peaks %>%
    summarise(
      N_peaks = n(),
      Top1pct = sum(
        Restricted_top_1pct,
        na.rm = TRUE
      ),
      Top01pct = sum(
        `Restricted_top_0.1pct`,
        na.rm = TRUE
      ),
      Physical_overlap = sum(
        Primary_region_recovered_by_physical_overlap,
        na.rm = TRUE
      )
    )
)

cat("\nOutput directory:\n", OUT, "\n")
