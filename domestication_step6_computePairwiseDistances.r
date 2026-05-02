#################### compute distance matrix of pairwise distances

library(dplyr)
library(ggplot2)
library(ape)

# install ggtree if needed
if (!requireNamespace("ggtree", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("ggtree")
}
library(ggtree)

# -------------------------
# 1. define sample sets
# -------------------------

dom_dip <- data_filtered %>%
  filter(
    domesticated == "domesticated",
    !is.na(line),
    species == "majus"
  ) %>%
  mutate(
    group = "domesticated",
    unit = line,
    sample_id = paste(line, replicate_number, sep = "_")
  )

wild_dip <- wild_plot_samples %>%
  mutate(
    group = "wild",
    unit = as.character(location_pattern),
    sample_id = PlantID
  )

# -------------------------
# 2. choose loci to use
#    keep loci that are present in both datasets
#    and are not annotated as colour
#    (retain loci with effect == NA)
# -------------------------

non_colour_loci <- marker_info_filtered %>%
  filter(LocusName %in% names(dom_dip)) %>%
  filter(LocusName %in% names(wild_dip)) %>%
  filter(effect != "colour" | is.na(effect)) %>%
  pull(LocusName)

cat("Number of non-colour loci used:", length(non_colour_loci), "\n")

# -------------------------
# 3. combine diploid datasets
# -------------------------

all_dip <- bind_rows(
  dom_dip %>% select(group, unit, sample_id, all_of(non_colour_loci)),
  wild_dip %>% select(group, unit, sample_id, all_of(non_colour_loci))
)

cat("Number of domesticated individuals:", sum(all_dip$group == "domesticated"), "\n")
cat("Number of wild individuals:", sum(all_dip$group == "wild"), "\n")

# -------------------------
# 4. build diploid genotype matrix
# -------------------------

dip_mat <- as.matrix(all_dip[, non_colour_loci, drop = FALSE])
dip_mat[dip_mat %in% c(-9, -10)] <- NA
storage.mode(dip_mat) <- "numeric"
rownames(dip_mat) <- all_dip$sample_id

# -------------------------
# 5. pairwise diploid allele-level distance
#
# For genotypes g_i and g_j coded 0/1/2:
# p_i = g_i / 2
# p_j = g_j / 2
#
# per-locus distance =
#   p_i(1 - p_j) + (1 - p_i)p_j
#
# So:
#   0 vs 0 = 0
#   0 vs 1 = 0.5
#   0 vs 2 = 1
#   1 vs 1 = 0.5
#   1 vs 2 = 0.5
#   2 vs 2 = 0
# -------------------------

pairwise_allele_distance <- function(mat) {
  n <- nrow(mat)
  out <- matrix(NA_real_, n, n)
  rownames(out) <- rownames(mat)
  colnames(out) <- rownames(mat)
  
  for (i in seq_len(n)) {
    out[i, i] <- 0
    
    if (i < n) {
      xi <- mat[i, ]
      
      for (j in (i + 1):n) {
        xj <- mat[j, ]
        keep <- !is.na(xi) & !is.na(xj)
        
        if (sum(keep) == 0) {
          d <- NA_real_
        } else {
          p_i <- xi[keep] / 2
          p_j <- xj[keep] / 2
          d_locus <- p_i * (1 - p_j) + (1 - p_i) * p_j
          d <- mean(d_locus)
        }
        
        out[i, j] <- d
        out[j, i] <- d
      }
    }
  }
  
  out
}

dist_mat_dip <- pairwise_allele_distance(dip_mat)

cat("Distance matrix dimensions:", nrow(dist_mat_dip), "x", ncol(dist_mat_dip), "\n")
cat("Symmetric:", isTRUE(all.equal(dist_mat_dip, t(dist_mat_dip), check.attributes = FALSE)), "\n")
cat("Diagonal all zero:", all(diag(dist_mat_dip) == 0, na.rm = TRUE), "\n")

# -------------------------
# 6. build pairwise dataframe
# -------------------------

pair_idx <- which(upper.tri(dist_mat_dip), arr.ind = TRUE)

pairs_df <- data.frame(
  i = pair_idx[, 1],
  j = pair_idx[, 2],
  distance = dist_mat_dip[pair_idx]
)

pairs_df$group_i <- all_dip$group[pairs_df$i]
pairs_df$group_j <- all_dip$group[pairs_df$j]
pairs_df$unit_i  <- all_dip$unit[pairs_df$i]
pairs_df$unit_j  <- all_dip$unit[pairs_df$j]

# -------------------------
# 7. classify comparisons
# -------------------------

pairs_df$comparison <- NA_character_

pairs_df$comparison[
  pairs_df$group_i == "domesticated" &
    pairs_df$group_j == "domesticated" &
    pairs_df$unit_i == pairs_df$unit_j
] <- "Within domesticated lines"

pairs_df$comparison[
  pairs_df$group_i == "wild" &
    pairs_df$group_j == "wild" &
    pairs_df$unit_i == pairs_df$unit_j
] <- "Within wild populations"

pairs_df$comparison[
  pairs_df$group_i == "domesticated" &
    pairs_df$group_j == "domesticated" &
    pairs_df$unit_i != pairs_df$unit_j
] <- "Between domesticated lines"

pairs_df$comparison[
  pairs_df$group_i == "wild" &
    pairs_df$group_j == "wild" &
    pairs_df$unit_i != pairs_df$unit_j
] <- "Between wild populations"

pairs_df$comparison[
  (pairs_df$group_i == "wild" & pairs_df$group_j == "domesticated") |
    (pairs_df$group_i == "domesticated" & pairs_df$group_j == "wild")
] <- "Between wild and domesticated"

plot_df_dip <- pairs_df %>%
  filter(!is.na(comparison), !is.na(distance))

plot_df_dip$comparison <- factor(
  plot_df_dip$comparison,
  levels = c(
    "Within domesticated lines",
    "Within wild populations",
    "Between domesticated lines",
    "Between wild populations",
    "Between wild and domesticated"
  )
)

# -------------------------
# 8. summary table
# -------------------------

distance_summary_dip <- plot_df_dip %>%
  group_by(comparison) %>%
  summarise(
    mean_distance = mean(distance, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  )

print(distance_summary_dip)
print(table(plot_df_dip$comparison))

