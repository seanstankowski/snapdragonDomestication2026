# =========================================================
# Weir & Cockerham FST between magenta vs yellow wild samples
# using hierfstat
#
# Uses:
#   - wild_plot_samples
#   - red_pops as magenta group
#   - yellow_pops as yellow group
#   - all loci except LG 9 and LG 10
#
# Assumes these objects already exist:
#   wild_plot_samples
#   marker_info_filtered
#   data_filtered
#   red_pops
#   yellow_pops
# =========================================================

library(dplyr)
library(hierfstat)

# -------------------------
# 1. define loci: exclude LG 9 and 10
# -------------------------

marker_info_use <- marker_info_filtered %>%
  mutate(
    LG_num = suppressWarnings(as.numeric(gsub("[^0-9.]", "", LG)))
  )

fst_loci <- marker_info_use %>%
  filter(LocusName %in% names(data_filtered)) %>%
  filter(!(LG_num %in% c(9, 10))) %>%
  pull(LocusName)

cat("Number of loci used:", length(fst_loci), "\n")

# -------------------------
# 2. define samples
#    pop 1 = magenta
#    pop 2 = yellow
# -------------------------

fst_data <- wild_plot_samples %>%
  filter(location_pattern %in% c(red_pops, yellow_pops)) %>%
  mutate(
    pop = case_when(
      location_pattern %in% red_pops ~ 1L,
      location_pattern %in% yellow_pops ~ 2L
    )
  ) %>%
  select(pop, all_of(fst_loci))

cat("Individuals in magenta:", sum(fst_data$pop == 1), "\n")
cat("Individuals in yellow:", sum(fst_data$pop == 2), "\n")

# -------------------------
# 3. convert 0/1/2 genotypes to hierfstat coding
#
# 0 -> 11
# 1 -> 12
# 2 -> 22
# -9/-10/NA -> NA
# -------------------------

convert_to_hierfstat <- function(x) {
  x[x %in% c(-9, -10)] <- NA
  ifelse(
    x == 0, 11,
    ifelse(
      x == 1, 12,
      ifelse(x == 2, 22, NA)
    )
  )
}

fst_geno <- fst_data

for (locus in fst_loci) {
  fst_geno[[locus]] <- convert_to_hierfstat(fst_geno[[locus]])
}

fst_geno <- as.data.frame(fst_geno)

# -------------------------
# 4. quick checks
# -------------------------

cat("Input dimensions:", nrow(fst_geno), "individuals x", ncol(fst_geno), "columns\n")
cat("Number of missing genotypes:", sum(is.na(fst_geno[, fst_loci])), "\n")

print(table(unlist(fst_geno[, fst_loci]), useNA = "ifany"))

# -------------------------
# 5. compute Weir & Cockerham FST
# -------------------------

fst_res <- hierfstat::wc(fst_geno)

# -------------------------
# 6. extract per-locus FST
# -------------------------

fst_per_locus <- fst_res$per.loc[, "FST"]

fst_table <- data.frame(
  locus = fst_loci,
  FST = as.numeric(fst_per_locus)
) %>%
  left_join(
    marker_info_use %>%
      select(LocusName, LG, cM, effect, gene_id),
    by = c("locus" = "LocusName")
  )

# -------------------------
# 7. overall summary
# -------------------------

cat("Mean per-locus FST:", mean(fst_table$FST, na.rm = TRUE), "\n")
cat("Median per-locus FST:", median(fst_table$FST, na.rm = TRUE), "\n")
cat("Overall FST from wc():", fst_res$FST, "\n")

# -------------------------
# 8. write output
# -------------------------

write.csv(
  fst_table,
  "fst_magenta_vs_yellow_per_locus.csv",
  row.names = FALSE
)

cat("Wrote: fst_magenta_vs_yellow_per_locus.csv\n")

# -------------------------
# colour vector: colour genes = red
# -------------------------

bar_cols <- ifelse(
  fst_table$effect == "colour",
  "firebrick",
  "grey70"
)

# -------------------------
# plot
# -------------------------

plot(
  fst_table$FST,
  ylim = c(0, 1),
  pch = 20,
  type = "h",
  lwd = 4,
  col = bar_cols
)

# optional: add legend
legend(
  "topright",
  legend = c("colour gene", "other"),
  col = c("firebrick", "grey70"),
  lwd = 4,
  bty = "n"
)