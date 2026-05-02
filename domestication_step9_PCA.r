# =========================================================
# PCA of wild + domesticated samples
# using non-colour loci only
# Missing genotypes imputed by random draws from the same
# wild population or domesticated line
# =========================================================

library(dplyr)
library(ggplot2)

# -------------------------
# 1. Define non-colour loci
# -------------------------

non_colour_loci <- marker_info_filtered %>%
  filter(LocusName %in% names(data_filtered)) %>%
  filter(effect != "colour" | is.na(effect)) %>%
  pull(LocusName)

cat("Number of non-colour loci:", length(non_colour_loci), "\n")

# -------------------------
# 2. Build PCA input data
# -------------------------

pca_wild <- wild_plot_samples %>%
  mutate(
    group = "wild",
    unit = as.character(location_pattern),
    sample_id = PlantID
  ) %>%
  select(sample_id, group, unit, all_of(non_colour_loci))

pca_dom <- dom_data %>%
  filter(species == "majus") %>%
  mutate(
    group = "domesticated",
    unit = as.character(line),
    sample_id = paste(line, replicate_number, sep = "_")
  ) %>%
  select(sample_id, group, unit, all_of(non_colour_loci))

pca_data <- bind_rows(pca_wild, pca_dom)

# -------------------------
# 3. Build genotype matrix
# -------------------------

pca_mat <- as.matrix(pca_data[, non_colour_loci, drop = FALSE])
pca_mat[pca_mat %in% c(-9, -10)] <- NA
storage.mode(pca_mat) <- "numeric"
rownames(pca_mat) <- pca_data$sample_id

# -------------------------
# 4. Impute missing genotypes
#    by random genotype drawn from same population/line
# -------------------------

set.seed(792)

impute_group_random <- function(mat, groups) {
  out <- mat
  
  for (g in unique(groups)) {
    idx <- which(groups == g)
    
    for (j in seq_len(ncol(out))) {
      missing_idx <- idx[is.na(out[idx, j])]
      
      if (length(missing_idx) > 0) {
        donor_values <- out[idx, j]
        donor_values <- donor_values[!is.na(donor_values)]
        
        if (length(donor_values) > 0) {
          out[missing_idx, j] <- sample(
            donor_values,
            size = length(missing_idx),
            replace = TRUE
          )
        }
      }
    }
  }
  
  out
}

pca_mat <- impute_group_random(
  mat = pca_mat,
  groups = pca_data$unit
)

# fallback: if an entire population/line is missing at a locus,
# impute from the global locus distribution
for (j in seq_len(ncol(pca_mat))) {
  missing <- is.na(pca_mat[, j])
  
  if (any(missing)) {
    donor_values <- pca_mat[, j]
    donor_values <- donor_values[!is.na(donor_values)]
    
    if (length(donor_values) > 0) {
      pca_mat[missing, j] <- sample(
        donor_values,
        size = sum(missing),
        replace = TRUE
      )
    }
  }
}

# drop loci still all missing, if any
keep_loci <- colSums(is.na(pca_mat)) == 0
pca_mat <- pca_mat[, keep_loci, drop = FALSE]

cat("Number of samples:", nrow(pca_mat), "\n")
cat("Number of loci used in PCA:", ncol(pca_mat), "\n")
cat("Remaining missing genotypes:", sum(is.na(pca_mat)), "\n")

# -------------------------
# 5. Run PCA
# -------------------------

pca <- prcomp(
  pca_mat,
  center = TRUE,
  scale. = FALSE
)

var_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_df <- data.frame(
  sample_id = pca_data$sample_id,
  group = pca_data$group,
  unit = pca_data$unit,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3]
)

# -------------------------
# 6. Plot PC1 vs PC2
# -------------------------

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(colour = group),
    size = 2,
    alpha = 0.8
  ) +
  scale_colour_manual(
    values = c(
      "wild" = "steelblue",
      "domesticated" = "firebrick"
    )
  ) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    colour = "Group",
    title = "PCA of wild and domesticated samples using non-colour loci"
  )