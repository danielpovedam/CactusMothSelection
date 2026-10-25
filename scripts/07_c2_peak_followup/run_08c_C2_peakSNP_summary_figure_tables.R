# ============================================================
# 61_FINAL_RAW_FigureS9_and_tables.R
#
# FINAL SCIENTIFICALLY TRANSPARENT VERSION
#
# IMPORTANT:
# - Uses ONLY raw BayPass posterior population allele
#   frequencies (M_P).
# - NO allele orientation.
# - NO p -> 1-p transformation.
#
# Figure S9
#   A. Native Opuntia vs O. ficus-indica boxplots
#   B. Population-level heatmap
#
# Supplementary tables
#   1. Candidate-level raw-frequency summary
#   2. Complete 17 x 28 population-level raw-frequency table
# ============================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# ============================================================
# 1. PATHS
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

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. READ DATA
# ============================================================

dat <- fread(INFILE) %>%
  as_tibble()

stopifnot(
  nrow(dat) == 476,
  n_distinct(dat$MRK) == 17,
  n_distinct(dat$PopulationID) == 28
)

stopifnot("M_P" %in% names(dat))

# ============================================================
# 3. CLEAN CANDIDATE LABELS
# ============================================================

dat <- dat %>%
  mutate(
    Candidate = case_when(
      Candidate_label == "ptg000007l_1_1" ~ "ptg000007",
      Candidate_label == "ptg000105l_1_1" ~ "ptg000105",
      Candidate_label == "ptg000171l_1_1" ~ "ptg000171",
      TRUE ~ Candidate_label
    )
  )

# ============================================================
# 4. FINAL CANDIDATE ORDER
# ============================================================

candidate_order <- c(
  "Hspg2",
  "Masc",
  "TBC1D13",
  "ptg000105",
  "CCHa1-R",
  "Orai",
  "Vinc",
  "ARFGEF3",
  "Usp8",
  "Dnajc9",
  "OXA1L",
  "Ccac_g10875",
  "Ccac_g3342",
  "Ccac_g1877",
  "ptg000007",
  "Ccac_g4727",
  "ptg000171"
)

stopifnot(
  setequal(
    unique(dat$Candidate),
    candidate_order
  )
)

# ============================================================
# 5. PARSED DISPLAY LABELS
# ============================================================

candidate_display <- c(
  "Hspg2"        = "italic('Hspg2')",
  "Masc"         = "italic('Masc')",
  "TBC1D13"      = "italic('TBC1D13')",
  "ptg000105"    = "plain('ptg000105')",
  "CCHa1-R"      = "italic('CCHa1-R')",
  "Orai"         = "italic('Orai')",
  "Vinc"         = "italic('Vinc')",
  "ARFGEF3"      = "italic('ARFGEF3')",
  "Usp8"         = "italic('Usp8')",
  "Dnajc9"       = "italic('Dnajc9')",
  "OXA1L"        = "italic('OXA1L')",
  "Ccac_g10875"  = "italic('Ccac_g10875')",
  "Ccac_g3342"   = "italic('Ccac_g3342')",
  "Ccac_g1877"   = "italic('Ccac_g1877')",
  "ptg000007"    = "plain('ptg000007')",
  "Ccac_g4727"   = "italic('Ccac_g4727')",
  "ptg000171"    = "plain('ptg000171')"
)

display_order <- unname(
  candidate_display[candidate_order]
)

# ============================================================
# 6. FACTOR LEVELS
# ============================================================

lineage_levels <- c(
  "Southern",
  "Eastern",
  "Central-northern-western"
)

host_levels <- c(
  "Native Opuntia",
  "O. ficus-indica"
)

dat <- dat %>%
  mutate(
    Candidate = factor(
      Candidate,
      levels = candidate_order
    ),

    Candidate_display =
      candidate_display[
        as.character(Candidate)
      ],

    Candidate_display = factor(
      Candidate_display,
      levels = display_order
    ),

    Host_category = factor(
      as.character(Host_category),
      levels = host_levels
    ),

    Lineage = factor(
      as.character(Lineage),
      levels = lineage_levels
    )
  )

# ============================================================
# 7. POPULATION ORDER FOR HEATMAP
#
# Host category -> lineage -> alphabetical population ID
#
# Completely independent of allele frequency.
# ============================================================

pop_order_df <- dat %>%
  distinct(
    PopulationID,
    Host_category,
    Lineage,
    Host
  ) %>%
  arrange(
    Host_category,
    Lineage,
    PopulationID
  ) %>%
  mutate(
    Population_order = row_number()
  )

pop_order <- pop_order_df$PopulationID

# ============================================================
# 8. RAW OVERALL HOST CONTRAST
#
# Delta_raw = mean(M_P_OFI) - mean(M_P_native)
#
# Sign may be positive OR negative.
# ============================================================

overall_raw <- dat %>%
  group_by(
    MRK,
    Candidate
  ) %>%
  summarise(
    N_native_overall =
      sum(Host_C2 == -1),

    N_OFI_overall =
      sum(Host_C2 == 1),

    Mean_native_overall =
      mean(
        M_P[Host_C2 == -1],
        na.rm = TRUE
      ),

    Mean_OFI_overall =
      mean(
        M_P[Host_C2 == 1],
        na.rm = TRUE
      ),

    Delta_overall_raw =
      Mean_OFI_overall -
      Mean_native_overall,

    Overall_direction =
      case_when(
        Delta_overall_raw > 0 ~ "OFI higher",
        Delta_overall_raw < 0 ~ "Native higher",
        TRUE ~ "Equal"
      ),

    .groups = "drop"
  )

# ============================================================
# 9. RAW WITHIN-LINEAGE CONTRASTS
# ============================================================

lineage_raw <- dat %>%
  group_by(
    MRK,
    Candidate,
    Lineage
  ) %>%
  summarise(
    N_native =
      sum(Host_C2 == -1),

    N_OFI =
      sum(Host_C2 == 1),

    Mean_native =
      mean(
        M_P[Host_C2 == -1],
        na.rm = TRUE
      ),

    Mean_OFI =
      mean(
        M_P[Host_C2 == 1],
        na.rm = TRUE
      ),

    Delta_raw =
      Mean_OFI -
      Mean_native,

    .groups = "drop"
  ) %>%
  left_join(
    overall_raw %>%
      select(
        MRK,
        Candidate,
        Delta_overall_raw
      ),
    by = c(
      "MRK",
      "Candidate"
    )
  ) %>%
  mutate(
    Same_direction_as_overall =
      case_when(
        Delta_overall_raw > 0 & Delta_raw > 0 ~ TRUE,
        Delta_overall_raw < 0 & Delta_raw < 0 ~ TRUE,
        Delta_overall_raw == 0 & Delta_raw == 0 ~ TRUE,
        TRUE ~ FALSE
      )
  )

# ============================================================
# 10. CONSISTENCY ACROSS LINEAGES
# ============================================================

consistency <- lineage_raw %>%
  group_by(
    MRK,
    Candidate
  ) %>%
  summarise(
    N_lineages_same_direction =
      sum(
        Same_direction_as_overall,
        na.rm = TRUE
      ),

    Consistent_all_3 =
      N_lineages_same_direction == 3,

    .groups = "drop"
  )

# ============================================================
# 11. PANEL A DATA
# ============================================================

boxdat <- dat %>%
  group_by(
    Candidate_display,
    Host_category
  ) %>%
  arrange(
    Lineage,
    PopulationID,
    .by_group = TRUE
  ) %>%
  mutate(
    point_offset = seq(
      -0.17,
      0.17,
      length.out = n()
    ),

    x_numeric =
      as.numeric(Host_category) +
      point_offset
  ) %>%
  ungroup()

# ============================================================
# 12. PANEL A
# RAW M_P VALUES
# ============================================================

p_box <- ggplot(boxdat) +

  geom_boxplot(
    aes(
      x = as.numeric(Host_category),
      y = M_P,
      group = Host_category,
      fill = Host_category
    ),
    width = 0.60,
    colour = "black",
    linewidth = 0.72,
    outlier.shape = NA,
    alpha = 0.80
  ) +

  geom_point(
    aes(
      x = x_numeric,
      y = M_P,
      shape = Lineage
    ),
    size = 2.35,
    stroke = 0.50,
    colour = "black",
    alpha = 0.86
  ) +

  facet_wrap(
    ~ Candidate_display,
    ncol = 5,
    labeller = label_parsed
  ) +

  scale_shape_manual(
    values = c(
      "Southern" = 16,
      "Eastern" = 17,
      "Central-northern-western" = 15
    ),
    name = "Major lineage"
  ) +

  scale_fill_manual(
    values = c(
      "Native Opuntia" = "#CFCFCF",
      "O. ficus-indica" = "#E78A32"
    ),
    guide = "none"
  ) +

  scale_x_continuous(
    breaks = c(1, 2),

    labels = list(
      expression(
        "Native "*italic(Opuntia)*" spp."
      ),
      expression(
        italic("O. ficus-indica")
      )
    ),

    limits = c(
      0.60,
      2.40
    )
  ) +

  scale_y_continuous(
    limits = c(0, 1),

    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),

    labels = c(
      "0.00",
      "0.25",
      "0.50",
      "0.75",
      "1.00"
    ),

    expand = expansion(
      mult = c(
        0.01,
        0.035
      )
    )
  ) +

  labs(
    x = NULL,
    y = "BayPass posterior allele frequency"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    strip.background =
      element_blank(),

    strip.text =
      element_text(
        size = 12.0,
        margin = margin(
          b = 5
        )
      ),

    axis.text.x =
      element_text(
        angle = 32,
        hjust = 1,
        vjust = 1,
        size = 10.0
      ),

    axis.text.y =
      element_text(
        size = 10.2
      ),

    axis.title.y =
      element_text(
        size = 11.7,
        margin = margin(
          r = 10
        )
      ),

    axis.ticks =
      element_line(
        linewidth = 0.48
      ),

    axis.line =
      element_line(
        linewidth = 0.55
      ),

    panel.spacing.x =
      unit(
        0.65,
        "cm"
      ),

    panel.spacing.y =
      unit(
        0.75,
        "cm"
      ),

    legend.position =
      "bottom",

    legend.direction =
      "horizontal",

    legend.title =
      element_text(
        size = 10.5,
        face = "bold"
      ),

    legend.text =
      element_text(
        size = 10.0
      ),

    legend.key.width =
      unit(
        0.62,
        "cm"
      ),

    legend.spacing.x =
      unit(
        0.24,
        "cm"
      ),

    plot.margin =
      margin(
        4,
        6,
        2,
        6
      )
  )

# ============================================================
# 13. PANEL B DATA
# RAW M_P VALUES
# ============================================================

heatdat <- dat %>%

  left_join(
    pop_order_df %>%
      select(
        PopulationID,
        Population_order
      ),
    by = "PopulationID"
  ) %>%

  mutate(
    PopulationID = factor(
      PopulationID,
      levels = pop_order
    ),

    Candidate_heat = factor(
      as.character(
        Candidate_display
      ),
      levels = rev(
        display_order
      )
    ),

    Host_group = case_when(
      Host_category ==
        "Native Opuntia" ~
        "Native~italic(Opuntia)~'spp.'",

      Host_category ==
        "O. ficus-indica" ~
        "italic('O. ficus-indica')"
    ),

    Host_group = factor(
      Host_group,
      levels = c(
        "Native~italic(Opuntia)~'spp.'",
        "italic('O. ficus-indica')"
      )
    ),

    Lineage_group = factor(
      as.character(Lineage),
      levels = lineage_levels
    )
  )

# ============================================================
# 14. PANEL B
# RAW M_P HEATMAP
# ============================================================

p_heat <- ggplot(
  heatdat,
  aes(
    x = PopulationID,
    y = Candidate_heat,
    fill = M_P
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.32
  ) +

  facet_grid(
    cols = vars(
      Host_group,
      Lineage_group
    ),
    scales = "free_x",
    space = "free_x",

    labeller = labeller(
      Host_group = label_parsed,
      Lineage_group = label_value
    )
  ) +

  scale_fill_gradientn(
    colours = c(
      "#FFFFCC",
      "#C7E9B4",
      "#7FCDBB",
      "#41B6C4",
      "#2C7FB8",
      "#253494"
    ),

    values = c(
      0.00,
      0.20,
      0.40,
      0.60,
      0.80,
      1.00
    ),

    limits = c(
      0,
      1
    ),

    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),

    labels = c(
      "0.00",
      "0.25",
      "0.50",
      "0.75",
      "1.00"
    ),

    name =
      "Posterior allele\nfrequency"
  ) +

  scale_x_discrete(
    drop = TRUE,
    expand = c(
      0,
      0
    )
  ) +

  scale_y_discrete(
    expand = c(
      0,
      0
    ),

    labels = function(x) {
      parse(
        text = x
      )
    }
  ) +

  labs(
    x = "Population",
    y = NULL
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    axis.text.x =
      element_text(
        angle = 55,
        hjust = 1,
        vjust = 1,
        size = 9.2
      ),

    axis.text.y =
      element_text(
        size = 10.2
      ),

    axis.title.x =
      element_text(
        size = 11.0,
        margin = margin(
          t = 7
        )
      ),

    axis.ticks =
      element_blank(),

    axis.line =
      element_blank(),

    panel.spacing.x =
      unit(
        0.11,
        "cm"
      ),

    strip.background =
      element_rect(
        fill = "white",
        colour = NA
      ),

    strip.text.x =
      element_text(
        size = 9.6,
        face = "bold",
        margin = margin(
          t = 3,
          b = 4
        )
      ),

    legend.position =
      "right",

    legend.title =
      element_text(
        size = 9.7,
        face = "bold",
        lineheight = 0.95
      ),

    legend.text =
      element_text(
        size = 9.0
      ),

    legend.margin =
      margin(
        l = 4
      ),

    plot.margin =
      margin(
        3,
        6,
        4,
        6
      )
  ) +

  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(
        3.7,
        "cm"
      ),
      barwidth = unit(
        0.42,
        "cm"
      ),
      ticks = TRUE,
      frame.colour = NA
    )
  )

# ============================================================
# 15. FINAL COMPOSITE
# ============================================================

g_box <- ggplotGrob(
  p_box
)

g_heat <- ggplotGrob(
  p_heat
)

draw_final <- function() {

  grid.newpage()

  lay <- grid.layout(
    nrow = 4,
    ncol = 1,

    heights = unit.c(
      unit(0.48, "in"),
      unit(8.35, "in"),
      unit(0.48, "in"),
      unit(5.00, "in")
    )
  )

  pushViewport(
    viewport(
      layout = lay
    )
  )

  grid.text(
    "A   Host-associated allele-frequency profiles",
    x = unit(0.006, "npc"),
    y = unit(0.52, "npc"),
    just = c("left", "center"),
    gp = gpar(
      fontsize = 14.0,
      fontface = "bold"
    ),
    vp = viewport(
      layout.pos.row = 1
    )
  )

  grid.draw(
    editGrob(
      g_box,
      vp = viewport(
        layout.pos.row = 2
      )
    )
  )

  grid.text(
    "B   Population distribution of candidate alleles",
    x = unit(0.006, "npc"),
    y = unit(0.52, "npc"),
    just = c("left", "center"),
    gp = gpar(
      fontsize = 14.0,
      fontface = "bold"
    ),
    vp = viewport(
      layout.pos.row = 3
    )
  )

  grid.draw(
    editGrob(
      g_heat,
      vp = viewport(
        layout.pos.row = 4
      )
    )
  )

  popViewport()
}

# ============================================================
# 16. EXPORT FIGURE
# ============================================================

FIG_WIDTH <- 16.0
FIG_HEIGHT <- 14.7

pdf(
  file.path(
    OUT,
    "FigureS9_C2_RAW_posterior_frequencies_FINAL.pdf"
  ),
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  useDingbats = FALSE
)

draw_final()

dev.off()

svg(
  filename = file.path(
    OUT,
    "FigureS9_C2_RAW_posterior_frequencies_FINAL.svg"
  ),
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  pointsize = 11,
  onefile = TRUE,
  bg = "white"
)

draw_final()

dev.off()

tiff(
  filename = file.path(
    OUT,
    "FigureS9_C2_RAW_posterior_frequencies_FINAL.tiff"
  ),
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  units = "in",
  res = 600,
  compression = "lzw",
  bg = "white"
)

draw_final()

dev.off()

png(
  filename = file.path(
    OUT,
    "FigureS9_C2_RAW_posterior_frequencies_FINAL.png"
  ),
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  units = "in",
  res = 400,
  bg = "white"
)

draw_final()

dev.off()

# ============================================================
# 17. FINAL CANDIDATE SUMMARY TABLE
# ============================================================

lineage_wide <- lineage_raw %>%
  mutate(
    Lineage_short = case_when(
      Lineage == "Southern" ~ "Southern",
      Lineage == "Eastern" ~ "Eastern",
      Lineage ==
        "Central-northern-western" ~ "CNW"
    )
  ) %>%
  select(
    MRK,
    Candidate,
    Lineage_short,
    N_native,
    N_OFI,
    Mean_native,
    Mean_OFI,
    Delta_raw,
    Same_direction_as_overall
  ) %>%
  pivot_wider(
    names_from = Lineage_short,
    values_from = c(
      N_native,
      N_OFI,
      Mean_native,
      Mean_OFI,
      Delta_raw,
      Same_direction_as_overall
    ),
    names_glue =
      "{.value}_{Lineage_short}"
  )

candidate_meta <- dat %>%
  distinct(
    MRK,
    Candidate,
    Primary_Region_ID,
    Scaffold_norm,
    Primary_Statistic_peak_position,
    Primary_Statistic_peak_value,
    Restricted_C2,
    Restricted_rank,
    Restricted_percentile,
    Restricted_top_1pct,
    `Restricted_top_0.1pct`,
    Primary_region_recovered_by_physical_overlap
  )

table1 <- candidate_meta %>%

  left_join(
    overall_raw,
    by = c(
      "MRK",
      "Candidate"
    )
  ) %>%

  left_join(
    lineage_wide,
    by = c(
      "MRK",
      "Candidate"
    )
  ) %>%

  left_join(
    consistency,
    by = c(
      "MRK",
      "Candidate"
    )
  ) %>%

  mutate(
    Candidate_order =
      match(
        as.character(Candidate),
        candidate_order
      )
  ) %>%

  arrange(
    Candidate_order
  ) %>%

  transmute(
    Candidate =
      as.character(Candidate),

    Region_ID =
      Primary_Region_ID,

    Scaffold =
      Scaffold_norm,

    Peak_position =
      Primary_Statistic_peak_position,

    C2_peak =
      Primary_Statistic_peak_value,

    N_native_overall,
    N_OFI_overall,

    Mean_native_overall,
    Mean_OFI_overall,

    Delta_overall_raw,

    Overall_direction,

    N_native_Southern,
    N_OFI_Southern,
    Mean_native_Southern,
    Mean_OFI_Southern,
    Delta_raw_Southern,
    Same_direction_as_overall_Southern,

    N_native_Eastern,
    N_OFI_Eastern,
    Mean_native_Eastern,
    Mean_OFI_Eastern,
    Delta_raw_Eastern,
    Same_direction_as_overall_Eastern,

    N_native_CNW,
    N_OFI_CNW,
    Mean_native_CNW,
    Mean_OFI_CNW,
    Delta_raw_CNW,
    Same_direction_as_overall_CNW,

    N_lineages_same_direction,

    Consistent_all_3,

    Restricted_C2,

    Restricted_percentile,

    Restricted_top_1pct =
      ifelse(
        Restricted_top_1pct,
        "Yes",
        "No"
      ),

    Restricted_top_0.1pct =
      ifelse(
        `Restricted_top_0.1pct`,
        "Yes",
        "No"
      ),

    Restricted_region_overlap =
      ifelse(
        Primary_region_recovered_by_physical_overlap,
        "Yes",
        "No"
      )
  ) %>%

  mutate(
    across(
      where(is.numeric),
      ~ round(
        .x,
        4
      )
    )
  )

fwrite(
  table1,
  file.path(
    OUT,
    "TableS_C2_peakSNP_RAW_frequency_summary_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  table1,
  file.path(
    OUT,
    "TableS_C2_peakSNP_RAW_frequency_summary_FINAL.csv"
  ),
  quote = TRUE,
  na = "NA"
)

# ============================================================
# 18. FINAL POPULATION-LEVEL TABLE
# ============================================================

table2 <- dat %>%

  mutate(
    Candidate_order =
      match(
        as.character(Candidate),
        candidate_order
      )
  ) %>%

  arrange(
    Candidate_order,
    Host_category,
    Lineage,
    PopulationID
  ) %>%

  transmute(
    Candidate =
      as.character(Candidate),

    Region_ID =
      Primary_Region_ID,

    Scaffold =
      Scaffold_norm,

    Peak_position =
      Primary_Statistic_peak_position,

    C2_peak =
      Primary_Statistic_peak_value,

    BayPass_marker =
      MRK,

    Population =
      PopulationID,

    Locality,

    Province,

    Host_species =
      Host,

    Host_category =
      as.character(Host_category),

    Major_lineage =
      as.character(Lineage),

    BayPass_posterior_frequency =
      M_P,

    BayPass_posterior_SD =
      SD_P,

    Restricted_C2,

    Restricted_percentile,

    Restricted_top_1pct =
      ifelse(
        Restricted_top_1pct,
        "Yes",
        "No"
      ),

    Restricted_top_0.1pct =
      ifelse(
        `Restricted_top_0.1pct`,
        "Yes",
        "No"
      ),

    Restricted_region_overlap =
      ifelse(
        Primary_region_recovered_by_physical_overlap,
        "Yes",
        "No"
      )
  ) %>%

  mutate(
    across(
      c(
        C2_peak,
        BayPass_posterior_frequency,
        BayPass_posterior_SD,
        Restricted_C2,
        Restricted_percentile
      ),
      ~ round(
        .x,
        4
      )
    )
  )

fwrite(
  table2,
  file.path(
    OUT,
    "TableS_C2_peakSNP_population_RAW_frequencies_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  table2,
  file.path(
    OUT,
    "TableS_C2_peakSNP_population_RAW_frequencies_FINAL.csv"
  ),
  quote = TRUE,
  na = "NA"
)

# ============================================================
# 19. EXPORT POPULATION ORDER
# ============================================================

fwrite(
  pop_order_df,
  file.path(
    OUT,
    "Table_FigureS9_RAW_population_order.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

# ============================================================
# 20. FINAL SCIENTIFIC QC
# ============================================================

cat("\n")
cat("=================================================\n")
cat("FIGURE S9 RAW-FREQUENCY FINAL QC\n")
cat("=================================================\n")

cat(
  "Candidate regions: ",
  n_distinct(dat$MRK),
  "\n",
  sep = ""
)

cat(
  "Populations: ",
  n_distinct(dat$PopulationID),
  "\n",
  sep = ""
)

cat(
  "Observations: ",
  nrow(dat),
  "\n\n",
  sep = ""
)

cat("RAW OVERALL HOST CONTRASTS:\n")

print(
  overall_raw %>%
    arrange(
      match(
        as.character(Candidate),
        candidate_order
      )
    ),
  n = Inf
)

cat("\nRAW WITHIN-LINEAGE CONTRASTS:\n")

print(
  lineage_raw %>%
    select(
      Candidate,
      Lineage,
      Mean_native,
      Mean_OFI,
      Delta_raw,
      Same_direction_as_overall
    ) %>%
    arrange(
      match(
        as.character(Candidate),
        candidate_order
      ),
      Lineage
    ),
  n = Inf
)

cat("\nDIRECTIONAL CONSISTENCY SUMMARY:\n")

print(
  consistency %>%
    arrange(
      match(
        as.character(Candidate),
        candidate_order
      )
    ),
  n = Inf
)

cat("\nNumber consistent in all 3 lineages: ")

cat(
  sum(
    consistency$Consistent_all_3
  ),
  "/17\n",
  sep = ""
)

cat("\nCounts of lineages matching overall direction:\n")

print(
  table(
    consistency$N_lineages_same_direction
  )
)

cat("\nFinal table sizes:\n")

cat(
  "Candidate summary: ",
  nrow(table1),
  " rows\n",
  sep = ""
)

cat(
  "Population table: ",
  nrow(table2),
  " rows\n",
  sep = ""
)

cat(
  "\nOutputs written to:\n",
  OUT,
  "\n"
)
