########## compute the per-sample heterozygosity and compare the means between magenta and yellow populations

library(dplyr)

# -------------------------
# 1. define locus sets
# -------------------------

marker_info_use <- marker_info_filtered %>%
  mutate(
    LG_num = suppressWarnings(as.numeric(gsub("[^0-9.]", "", LG)))
  )

# non-colour loci
non_colour_loci <- marker_info_use %>%
  filter(LocusName %in% names(data_filtered)) %>%
  filter(effect != "colour" | is.na(effect)) %>%
  pull(LocusName)

# all loci used in the main genotype plots
# assuming "all 103 loci" = exclude LG 10 only
all103_loci <- marker_info_use %>%
  filter(LocusName %in% names(data_filtered)) %>%
  filter(LG_num != 10) %>%
  pull(LocusName)

# if instead you really want to exclude both LG 9 and LG 10, use this:
# all103_loci <- marker_info_use %>%
#   filter(LocusName %in% names(data_filtered)) %>%
#   filter(!(LG_num %in% c(9, 10))) %>%
#   pull(LocusName)

cat("Number of non-colour loci:", length(non_colour_loci), "\n")
cat("Number of all103 loci:", length(all103_loci), "\n")

# -------------------------
# 2. define wild samples
# -------------------------

wild_het <- wild_plot_samples %>%
  mutate(
    group = as.character(location_pattern),
    group_type = "wild",
    sample_id = PlantID
  ) %>%
  select(group_type, group, sample_id, everything())

# -------------------------
# 3. define domesticated samples
#    majus only
# -------------------------

dom_het <- data_filtered %>%
  filter(
    domesticated == "domesticated",
    !is.na(line),
    species == "majus"
  ) %>%
  mutate(
    group = as.character(line),
    group_type = "domesticated",
    sample_id = paste(line, replicate_number, sep = "_")
  ) %>%
  select(group_type, group, sample_id, everything())

# combine
het_data <- bind_rows(
  wild_het,
  dom_het
)

# -------------------------
# 4. helper to compute per-sample heterozygosity
# -------------------------

compute_prop_het <- function(df, loci) {
  geno_mat <- as.matrix(df[, loci, drop = FALSE])
  geno_mat[geno_mat %in% c(-9, -10)] <- NA
  storage.mode(geno_mat) <- "numeric"
  
  apply(geno_mat, 1, function(x) {
    mean(x == 1, na.rm = TRUE)
  })
}

# -------------------------
# 5. compute per-sample heterozygosity for both locus sets
# -------------------------

het_data$prop_het_non_colour <- compute_prop_het(het_data, non_colour_loci)
het_data$prop_het_all103 <- compute_prop_het(het_data, all103_loci)

cat("Samples with NA prop_het_non_colour:", sum(is.na(het_data$prop_het_non_colour)), "\n")
cat("Samples with NA prop_het_all103:", sum(is.na(het_data$prop_het_all103)), "\n")

# -------------------------
# 6. summarise by wild population / domesticated line
# -------------------------

het_summary <- het_data %>%
  group_by(group_type, group) %>%
  summarise(
    mean_prop_het_non_colour = mean(prop_het_non_colour, na.rm = TRUE),
    sd_prop_het_non_colour = sd(prop_het_non_colour, na.rm = TRUE),
    mean_prop_het_all103 = mean(prop_het_all103, na.rm = TRUE),
    sd_prop_het_all103 = sd(prop_het_all103, na.rm = TRUE),
    n_individuals = sum(!is.na(prop_het_non_colour) | !is.na(prop_het_all103)),
    .groups = "drop"
  ) %>%
  mutate(
    wild_colour_group = case_when(
      group_type == "wild" & group %in% red_pops ~ "magenta",
      group_type == "wild" & group %in% yellow_pops ~ "yellow",
      group_type == "wild" ~ "other",
      TRUE ~ NA_character_
    )
  )

# optional: make colour group an ordered factor
het_summary$wild_colour_group <- factor(
  het_summary$wild_colour_group,
  levels = c("magenta", "yellow", "other")
)

# -------------------------
# 7. order output
# -------------------------

custom_order_wild <- c(
  "BAR****", "BES****", "CAL****", "CLA****", "D_***_R*", "GARL****",
  "GCHI****", "GPER****", "GPUI****", "M406*", "MP113*", "MRT****",
  "PLAMN", "PLAMF", "SUE****", "VIE****", "D27**", "PLAHZ", "CADI****",
  "ELS****", "GBOU****", "GCAM****", "GDOS****", "GFAB****", "LU****",
  "MLYS****", "PAL****", "POM****", "YP044*", "PLAYF", "PLAYN",
  "USS****", "VIL****"
)

dom_order <- het_summary %>%
  filter(group_type == "domesticated") %>%
  pull(group) %>%
  unique()

dom_order <- dom_order[order(as.numeric(sub("^L", "", dom_order)))]

het_summary <- het_summary %>%
  mutate(
    plot_order = case_when(
      group_type == "wild" ~ match(group, custom_order_wild),
      group_type == "domesticated" ~ match(group, dom_order) + 1000
    )
  ) %>%
  arrange(plot_order) %>%
  select(-plot_order)

print(het_summary)

# -------------------------
# 8. write output
# -------------------------

write.csv(
  het_summary,
  file = "heterozygosity_summary_by_group.csv",
  row.names = FALSE
)

cat("Wrote: heterozygosity_summary_by_group.csv\n")


# -------------------------
# 9. Mann–Whitney U tests (wild populations only)
# -------------------------

test_data <- het_summary %>%
  filter(
    group_type == "wild",
    wild_colour_group %in% c("magenta", "yellow")
  )

cat("Number of magenta populations:", sum(test_data$wild_colour_group == "magenta"), "\n")
cat("Number of yellow populations:", sum(test_data$wild_colour_group == "yellow"), "\n")

# -------------------------
# 9a. non-colour loci
# -------------------------

wilcox_non_colour <- wilcox.test(
  mean_prop_het_non_colour ~ wild_colour_group,
  data = test_data,
  exact = FALSE
)

print(wilcox_non_colour)

# -------------------------
# 9b. all 103 loci
# -------------------------

wilcox_all103 <- wilcox.test(
  mean_prop_het_all103 ~ wild_colour_group,
  data = test_data,
  exact = FALSE
)

print(wilcox_all103)

# -------------------------
# 10. report summary
# -------------------------

cat("\n--- Summary ---\n")

cat("Non-colour loci:\n")
cat("  W =", wilcox_non_colour$statistic,
    " p =", wilcox_non_colour$p.value, "\n")

cat("All loci (103):\n")
cat("  W =", wilcox_all103$statistic,
    " p =", wilcox_all103$p.value, "\n")