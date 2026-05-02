# calculate FIS
##############################
library(dplyr)
library(hierfstat)

# Fis within populations / lines
# using the same wild populations as the Nei analysis
# and the same domesticated lines
# calculate:
#   1. Fis for all loci except LG 9 and 10
#   2. Fis for all non-colour loci except LG 9 and 10
# output both in the same table
# with sample_type indicating wild or domesticated
##############################

# -------------------------
# User options
# -------------------------
dom_species_to_keep <- "majus"

custom_order_wild <- c(
  "BAR****",
  "BES****",
  "CADI****",
  "CAL****",
  "CLA****",
  "D27**",
  "D_***_R*",
  "ELS****",
  "GARL****",
  "GBOU****",
  "GCAM****",
  "GCHI****",
  "GDOS****",
  "GFAB****",
  "GPER****",
  "GPUI****",
  "LU****",
  "M406*",
  "MLYS****",
  "MP113*",
  "MRT****",
  "PAL****",
  "POM****",
  "YP044*",
  "PLAYF",
  "PLAYN",
  "PLAHZ",
  "PLAMN",
  "PLAMF",
  "SUE****",
  "USS****",
  "VIE****",
  "VIL****"
)

# -------------------------
# Helper functions
# -------------------------
dosage_to_hierfstat <- function(x) {
  x[x %in% c(-9, -10)] <- NA
  out <- rep(NA_integer_, length(x))
  out[x == 0] <- 11L
  out[x == 1] <- 12L
  out[x == 2] <- 22L
  out
}

calc_fis_for_marker_set <- function(marker_info_use, fis_samples, group_order, suffix_label) {
  marker_info_use <- marker_info_use %>%
    mutate(
      LG_num = suppressWarnings(as.numeric(gsub("[^0-9.]", "", LG))),
      cM_num = suppressWarnings(as.numeric(cM))
    ) %>%
    filter(!(LG_num %in% c(9, 10))) %>%
    arrange(LG_num, cM_num)

  geno_cols <- marker_info_use$LocusName
  geno_cols <- geno_cols[geno_cols %in% names(fis_samples)]

  cat("\n", suffix_label, "\n", sep = "")
  cat("Loci before monomorphic filter:", length(geno_cols), "\n")

  geno_use <- fis_samples[, geno_cols, drop = FALSE]

  keep_loci <- sapply(geno_use, function(x) {
    x[x %in% c(-9, -10)] <- NA
    length(unique(x[!is.na(x)])) > 1
  })

  geno_cols_use <- names(keep_loci)[keep_loci]

  cat("Loci retained after monomorphic filter:", length(geno_cols_use), "\n")

  samples_use <- fis_samples %>%
    select(PlantID, group, all_of(geno_cols_use))

  hf_geno <- samples_use[, geno_cols_use, drop = FALSE]
  for (j in seq_along(hf_geno)) {
    hf_geno[[j]] <- dosage_to_hierfstat(hf_geno[[j]])
  }

  pop_factor <- factor(samples_use$group, levels = group_order)

  hf_data <- data.frame(
    pop = as.integer(pop_factor),
    hf_geno,
    check.names = FALSE
  )

  bs <- basic.stats(hf_data)

  Ho_per_pop <- colMeans(bs$Ho, na.rm = TRUE)
  Hs_per_pop <- colMeans(bs$Hs, na.rm = TRUE)
  Fis_per_pop <- 1 - (Ho_per_pop / Hs_per_pop)

  stopifnot(length(group_order) == length(Ho_per_pop))

  data.frame(
    group = group_order,
    Ho_mean = as.numeric(Ho_per_pop),
    Hs_mean = as.numeric(Hs_per_pop),
    Fis = as.numeric(Fis_per_pop),
    row.names = NULL
  )
}

# -------------------------
# 1. Build wild sample set
# -------------------------
marker_info_all_tmp <- marker_info_filtered %>%
  mutate(LG_num = suppressWarnings(as.numeric(gsub("[^0-9.]", "", LG)))) %>%
  filter(!(LG_num %in% c(9, 10))) %>%
  arrange(LG_num)

geno_cols_all <- marker_info_all_tmp$LocusName
geno_cols_all <- geno_cols_all[geno_cols_all %in% names(data_filtered)]

wild_fis <- wild_plot_samples %>%
  select(PlantID, location_pattern) %>%
  left_join(
    data_filtered %>%
      select(PlantID, all_of(geno_cols_all)),
    by = "PlantID"
  ) %>%
  mutate(group = location_pattern) %>%
  select(PlantID, group, all_of(geno_cols_all))

# -------------------------
# 2. Build domesticated sample set
# -------------------------
dom_fis <- data_filtered %>%
  filter(domesticated == "domesticated", !is.na(line))

if (!is.null(dom_species_to_keep)) {
  dom_fis <- dom_fis %>%
    filter(species == dom_species_to_keep)
}

dom_fis <- dom_fis %>%
  mutate(group = line) %>%
  select(PlantID, group, all_of(geno_cols_all))

# -------------------------
# 3. Combine
# -------------------------
fis_samples <- bind_rows(wild_fis, dom_fis)

cat("Individuals used:", nrow(fis_samples), "\n")
cat("Groups used:", length(unique(fis_samples$group)), "\n")

group_sizes <- fis_samples %>%
  count(group, sort = FALSE)

print(group_sizes)

# -------------------------
# 4. Build group order
# -------------------------
dom_line_order <- data_filtered %>%
  filter(domesticated == "domesticated", !is.na(line)) %>%
  mutate(
    species = factor(
      species,
      levels = c("majus", "majus_ssp_tortuosum", "molle", "siculum")
    ),
    line_num = suppressWarnings(as.numeric(sub("^L", "", line)))
  ) %>%
  distinct(line, species, line_num) %>%
  arrange(species, line_num) %>%
  pull(line)

if (!is.null(dom_species_to_keep)) {
  dom_line_order <- data_filtered %>%
    filter(
      domesticated == "domesticated",
      !is.na(line),
      species == dom_species_to_keep
    ) %>%
    mutate(
      species = factor(
        species,
        levels = c("majus", "majus_ssp_tortuosum", "molle", "siculum")
      ),
      line_num = suppressWarnings(as.numeric(sub("^L", "", line)))
    ) %>%
    distinct(line, species, line_num) %>%
    arrange(species, line_num) %>%
    pull(line)
}

group_order <- c(
  custom_order_wild[custom_order_wild %in% unique(wild_fis$group)],
  dom_line_order[dom_line_order %in% unique(dom_fis$group)]
)

# -------------------------
# 5. Calculate Fis for all loci
# -------------------------
fis_all <- calc_fis_for_marker_set(
  marker_info_use = marker_info_filtered,
  fis_samples = fis_samples,
  group_order = group_order,
  suffix_label = "All loci except LG 9 and 10"
) %>%
  rename(
    Ho_mean_all = Ho_mean,
    Hs_mean_all = Hs_mean,
    Fis_all = Fis
  )

# -------------------------
# 6. Calculate Fis for non-colour loci
# -------------------------
fis_noncolour <- calc_fis_for_marker_set(
  marker_info_use = marker_info_filtered %>% filter(is.na(effect) | effect != "colour"),
  fis_samples = fis_samples,
  group_order = group_order,
  suffix_label = "Non-colour loci except LG 9 and 10"
) %>%
  rename(
    Ho_mean_noncolour = Ho_mean,
    Hs_mean_noncolour = Hs_mean,
    Fis_noncolour = Fis
  )

# -------------------------
# 7. Merge results
# -------------------------
fis_results <- data.frame(
  group = group_order,
  sample_type = ifelse(group_order %in% unique(wild_fis$group), "wild", "domesticated"),
  n_individuals = group_sizes$n[match(group_order, group_sizes$group)],
  row.names = NULL
) %>%
  left_join(fis_all, by = "group") %>%
  left_join(fis_noncolour, by = "group")

print(fis_results)

# -------------------------
# 8. Save output
# -------------------------
write.csv(
  fis_results,
  file = "Fis_by_population_all_and_noncolour_excluding_LG9_10.csv",
  row.names = FALSE
)
