library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# PATHS
# ============================================================

BASE <- "~/work/micetral/Population_snps/baypass"

INFILE <- file.path(
  BASE,
  "selection_scan/C2_peakSNP_BayPass_pij_host_profiles/",
  "Table_C2_17_peakSNP_BayPass_posterior_population_frequencies.tsv"
)

OUT <- file.path(
  BASE,
  "selection_scan/C2_peakSNP_BayPass_pij_host_profiles"
)

dat <- fread(INFILE) %>% as_tibble()

# ============================================================
# 1. BASIC QC
# ============================================================

stopifnot(n_distinct(dat$MRK) == 17)
stopifnot(n_distinct(dat$PopulationID) == 28)
stopifnot(nrow(dat) == 476)

cat("17 peaks x 28 populations confirmed\n\n")

# ============================================================
# 2. WITHIN-LINEAGE HOST CONTRASTS
#
# Descriptive only:
# Delta = mean oriented frequency in OFI populations
#       - mean oriented frequency in native-host populations
# ============================================================

lineage_delta <- dat %>%
  group_by(
    MRK,
    Primary_Region_ID,
    Candidate_label,
    Lineage
  ) %>%
  summarise(
    N_native = sum(Host_C2 == -1),
    N_OFI = sum(Host_C2 == 1),

    Mean_native = mean(
      M_P_oriented[Host_C2 == -1],
      na.rm = TRUE
    ),

    Mean_OFI = mean(
      M_P_oriented[Host_C2 == 1],
      na.rm = TRUE
    ),

    Delta_OFI_minus_native = Mean_OFI - Mean_native,

    .groups = "drop"
  ) %>%
  mutate(
    Direction = case_when(
      Delta_OFI_minus_native > 0 ~ "OFI higher",
      Delta_OFI_minus_native < 0 ~ "Native higher",
      TRUE ~ "Equal"
    )
  )

fwrite(
  lineage_delta,
  file.path(
    OUT,
    "Table_C2_17_peakSNP_within_lineage_host_differences.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 3. CONSISTENCY ACROSS THE THREE LINEAGES
# ============================================================

consistency <- lineage_delta %>%
  group_by(
    MRK,
    Primary_Region_ID,
    Candidate_label
  ) %>%
  summarise(
    Delta_Southern =
      Delta_OFI_minus_native[Lineage == "Southern"],

    Delta_Eastern =
      Delta_OFI_minus_native[Lineage == "Eastern"],

    Delta_CNW =
      Delta_OFI_minus_native[
        Lineage == "Central-northern-western"
      ],

    N_lineages_OFI_higher =
      sum(Delta_OFI_minus_native > 0),

    N_lineages_native_higher =
      sum(Delta_OFI_minus_native < 0),

    .groups = "drop"
  ) %>%
  mutate(
    Consistency = case_when(
      N_lineages_OFI_higher == 3 ~
        "OFI higher in all 3 lineages",

      N_lineages_OFI_higher == 2 ~
        "OFI higher in 2 of 3 lineages",

      N_lineages_OFI_higher == 1 ~
        "OFI higher in 1 of 3 lineages",

      TRUE ~
        "No lineage with OFI higher"
    )
  )

fwrite(
  consistency,
  file.path(
    OUT,
    "Table_C2_17_peakSNP_lineage_consistency_summary.tsv"
  ),
  sep = "\t"
)

cat("=== CONSISTENCY ACROSS K=3 LINEAGES ===\n")
print(
  consistency %>%
    select(
      Candidate_label,
      Delta_Southern,
      Delta_Eastern,
      Delta_CNW,
      N_lineages_OFI_higher,
      Consistency
    )
)

cat("\nOverall counts:\n")
print(table(consistency$Consistency))

# ============================================================
# 4. ADD SENSITIVITY SUPPORT
# ============================================================

candidate_info <- dat %>%
  distinct(
    MRK,
    Candidate_label,
    Restricted_percentile,
    Restricted_top_1pct,
    `Restricted_top_0.1pct`,
    Primary_region_recovered_by_physical_overlap
  ) %>%
  left_join(
    consistency,
    by = c(
      "MRK",
      "Candidate_label"
    )
  )

fwrite(
  candidate_info,
  file.path(
    OUT,
    "Table_C2_17_peakSNP_integrated_support_summary.tsv"
  ),
  sep = "\t"
)

# ============================================================
# 5. POPULATION ORDER
#
# Geographic lineage first.
# Within lineage: Native populations then OFI.
# Within native category: host species.
# ============================================================

lineage_levels <- c(
  "Southern",
  "Eastern",
  "Central-northern-western"
)

dat <- dat %>%
  mutate(
    Lineage = factor(
      Lineage,
      levels = lineage_levels
    ),

    Host_category = factor(
      Host_category,
      levels = c(
        "Native Opuntia",
        "O. ficus-indica"
      )
    )
  )

population_order <- dat %>%
  distinct(
    PopulationID,
    Lineage,
    Host_category,
    Host
  ) %>%
  arrange(
    Lineage,
    Host_category,
    Host,
    PopulationID
  ) %>%
  pull(PopulationID)

# Candidate order = genomic/BayPass order
candidate_order <- dat %>%
  distinct(
    MRK,
    Candidate_label
  ) %>%
  arrange(MRK) %>%
  pull(Candidate_label)

dat <- dat %>%
  mutate(
    PopulationID = factor(
      PopulationID,
      levels = population_order
    ),

    Candidate_label = factor(
      Candidate_label,
      levels = rev(candidate_order)
    )
  )

# ============================================================
# 6. MAIN HEATMAP
# ============================================================

p_heat <- ggplot(
  dat,
  aes(
    x = PopulationID,
    y = Candidate_label,
    fill = M_P_oriented
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Posterior allele\nfrequency"
  ) +
  labs(
    x = NULL,
    y = "C2 candidate peak"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.y = element_text(
      size = 8
    ),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

ggsave(
  file.path(
    OUT,
    "Figure_C2_17_peakSNP_population_heatmap_FINAL.pdf"
  ),
  p_heat,
  width = 12,
  height = 7.5
)

ggsave(
  file.path(
    OUT,
    "Figure_C2_17_peakSNP_population_heatmap_FINAL.png"
  ),
  p_heat,
  width = 12,
  height = 7.5,
  dpi = 500
)

# ============================================================
# 7. ANNOTATION DATA FOR POPULATIONS
# ============================================================

pop_annotation <- dat %>%
  distinct(
    PopulationID,
    Lineage,
    Host_category,
    Host
  ) %>%
  arrange(
    match(PopulationID, population_order)
  ) %>%
  mutate(
    PopulationID = factor(
      PopulationID,
      levels = population_order
    )
  )

# ------------------------------------------------------------
# Lineage annotation
# ------------------------------------------------------------

p_lineage <- ggplot(
  pop_annotation,
  aes(
    x = PopulationID,
    y = "Lineage",
    fill = Lineage
  )
) +
  geom_tile() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Major lineage"
  ) +
  theme_void(base_size = 9) +
  theme(
    legend.position = "right"
  )

# ------------------------------------------------------------
# Binary host annotation
# ------------------------------------------------------------

p_hostcat <- ggplot(
  pop_annotation,
  aes(
    x = PopulationID,
    y = "Host category",
    fill = Host_category
  )
) +
  geom_tile() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Host category"
  ) +
  theme_void(base_size = 9) +
  theme(
    legend.position = "right"
  )

# ------------------------------------------------------------
# Actual host species annotation
# ------------------------------------------------------------

p_hostspecies <- ggplot(
  pop_annotation,
  aes(
    x = PopulationID,
    y = "Host species",
    fill = Host
  )
) +
  geom_tile() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Host species"
  ) +
  theme_void(base_size = 9) +
  theme(
    legend.position = "right"
  )

# Save annotations individually too
ggsave(
  file.path(OUT, "annotation_lineage.pdf"),
  p_lineage,
  width = 12,
  height = 1.3
)

ggsave(
  file.path(OUT, "annotation_host_category.pdf"),
  p_hostcat,
  width = 12,
  height = 1.3
)

ggsave(
  file.path(OUT, "annotation_host_species.pdf"),
  p_hostspecies,
  width = 12,
  height = 1.3
)

# ============================================================
# 8. LINEAGE-CONTRAST HEATMAP
#
# This is particularly useful for Reviewer Comment 3.
# Rows = candidate peaks
# Columns = K=3 lineages
# Value = OFI mean - native mean
# ============================================================

lineage_delta <- lineage_delta %>%
  mutate(
    Candidate_label = factor(
      Candidate_label,
      levels = rev(candidate_order)
    ),

    Lineage = factor(
      Lineage,
      levels = lineage_levels
    )
  )

max_abs_delta <- max(
  abs(lineage_delta$Delta_OFI_minus_native),
  na.rm = TRUE
)

p_delta <- ggplot(
  lineage_delta,
  aes(
    x = Lineage,
    y = Candidate_label,
    fill = Delta_OFI_minus_native
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Delta_OFI_minus_native
      )
    ),
    size = 2.7
  ) +
  scale_fill_gradient2(
    low = "white",
    mid = "grey90",
    high = "black",
    midpoint = 0,
    limits = c(
      -max_abs_delta,
      max_abs_delta
    ),
    name = expression(
      Delta*" allele frequency"
    )
  ) +
  labs(
    x = "Major genetic lineage",
    y = "C2 candidate peak"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 8
    ),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(
    OUT,
    "Figure_C2_17_peakSNP_within_lineage_host_deltas.pdf"
  ),
  p_delta,
  width = 7,
  height = 7.5
)

ggsave(
  file.path(
    OUT,
    "Figure_C2_17_peakSNP_within_lineage_host_deltas.png"
  ),
  p_delta,
  width = 7,
  height = 7.5,
  dpi = 500
)

# ============================================================
# 9. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("REVIEWER COMMENT 3 - FINAL SUMMARY\n")
cat("========================================\n\n")

cat("Number of candidate peaks: ",
    n_distinct(dat$MRK), "\n")

cat("Number of populations: ",
    n_distinct(dat$PopulationID), "\n\n")

cat("Within-lineage consistency:\n")
print(
  consistency %>%
    select(
      Candidate_label,
      N_lineages_OFI_higher,
      Consistency
    )
)

cat("\nSensitivity + lineage consistency:\n")
print(
  candidate_info %>%
    select(
      Candidate_label,
      Restricted_percentile,
      Restricted_top_1pct,
      `Restricted_top_0.1pct`,
      Primary_region_recovered_by_physical_overlap,
      N_lineages_OFI_higher,
      Consistency
    )
)

cat("\nDone.\n")
