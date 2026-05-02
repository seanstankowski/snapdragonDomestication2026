########################################################################################################################
# Combined 0/1/2 plots, polarised by delta p between red and yellow populations:
#   1. wild populations
#   2. domesticated lines
# Polarisation rule:
#   genotype 2 = allele more common in red than yellow populations
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(stringr)
library(patchwork)

# -------------------------
# 0. Helper functions
# -------------------------
dms_to_decimal <- function(deg, min, sec, hemi) {
  x <- deg + min / 60 + sec / 3600
  if (hemi %in% c("S", "W")) x <- -x
  x
}

haversine_m <- function(lat1, lon1, lat2, lon2) {
  R <- 6371000
  to_rad <- pi / 180
  phi1 <- lat1 * to_rad
  phi2 <- lat2 * to_rad
  dphi <- (lat2 - lat1) * to_rad
  dlambda <- (lon2 - lon1) * to_rad
  a <- sin(dphi / 2)^2 +
    cos(phi1) * cos(phi2) * sin(dlambda / 2)^2
  2 * R * atan2(sqrt(a), sqrt(1 - a))
}

wildcard_to_regex <- function(x) {
  x <- gsub("([.|()\\^{}+$?])", "\\\\\\1", x)
  x <- gsub("\\*", ".", x)
  paste0("^", x, "$")
}

match_pattern <- function(id, regex_vec) {
  hits <- names(regex_vec)[str_detect(id, regex_vec)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

flip_012_cols <- function(df, loci) {
  for (m in loci) {
    x <- df[[m]]
    df[[m]] <- ifelse(x == 0, 2,
               ifelse(x == 2, 0, x))
  }
  df
}

# -------------------------
# 1. Wild population definitions
# -------------------------
location_patterns <- c(
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
  "SUE****",
  "USS****",
  "VIE****",
  "VIL****",
  "YP044*"
)

location_regex <- setNames(
  vapply(location_patterns, wildcard_to_regex, character(1)),
  location_patterns
)

custom_order_wild <- c(
  "BAR****",
  "BES****",
  "CAL****",
  "CLA****",
  "D_***_R*",
  "GARL****",
  "GCHI****",
  "GPER****",
  "GPUI****",
  "M406*",
  "MP113*",
  "MRT****",
  "PLAMN",
  "PLAMF",
  "SUE****",
  "VIE****",
  "D27**",
  "PLAHZ",
  "CADI****",
  "ELS****",
  "GBOU****",
  "GCAM****",
  "GDOS****",
  "GFAB****",
  "LU****",
  "MLYS****",
  "PAL****",
  "POM****",
  "YP044*",
  "PLAYF",
  "PLAYN",
  "USS****",
  "VIL****"
)

red_pops <- c(
  "BAR****",
  "BES****",
  "CAL****",
  "CLA****",
  "D_***_R*",
  "GARL****",
  "GCHI****",
  "GPER****",
  "GPUI****",
  "M406*",
  "MP113*",
  "MRT****",
  "PLAMN",
  "PLAMF",
  "SUE****",
  "VIE****"
)

yellow_pops <- c(
  "ELS****",
  "GBOU****",
  "GCAM****",
  "GDOS****",
  "GFAB****",
  "LU****",
  "MLYS****",
  "PAL****",
  "POM****",
  "YP044*",
  "PLAYF",
  "PLAYN",
  "USS****",
  "VIL****",
  "CADI****",
  "D27**"
)

hybrid_pops <- c("PLAHZ")

# -------------------------
# 2. Extra focal wild sites
# -------------------------
site_targets <- tibble::tribble(
  ~site,    ~lat,                                               ~lon,                                               ~n_target,
  "PLAYF",  dms_to_decimal(42, 19, 23.46, "N"),                 dms_to_decimal(2,  2, 50.77, "E"),                  10,
  "PLAYN",  dms_to_decimal(42, 19, 27.89, "N"),                 dms_to_decimal(2,  3, 53.45, "E"),                  10,
  "PLAMN",  dms_to_decimal(42, 19, 19.70, "N"),                 dms_to_decimal(2,  5, 15.80, "E"),                  10,
  "PLAMF",  dms_to_decimal(42, 19, 12.79, "N"),                 dms_to_decimal(2,  5, 49.83, "E"),                  10,
  "PLAHZ",  dms_to_decimal(42, 19, 21.43, "N"),                 dms_to_decimal(2,  4, 26.79, "E"),                  11
)

radius_m <- 100

# -------------------------
# 3. Build wild_plot_samples
# -------------------------
wild_all <- data_filtered %>%
  filter(
    domesticated == "wild",
    !is.na(PlantID),
    !is.na(CorrectedLatitude),
    !is.na(CorrectedLongitude)
  ) %>%
  mutate(
    location_pattern = vapply(
      PlantID,
      match_pattern,
      character(1),
      regex_vec = location_regex
    )
  )

set.seed(792)

wild_named_sampled <- wild_all %>%
  filter(!is.na(location_pattern)) %>%
  group_by(location_pattern) %>%
  group_modify(~ slice_sample(.x, n = min(10, nrow(.x)))) %>%
  ungroup()

wild_remaining <- wild_all %>%
  filter(!(PlantID %in% wild_named_sampled$PlantID))

site_samples_list <- vector("list", nrow(site_targets))
used_ids <- character(0)

for (i in seq_len(nrow(site_targets))) {
  candidates <- wild_remaining %>%
    filter(!(PlantID %in% used_ids)) %>%
    mutate(
      dist_m = haversine_m(
        CorrectedLatitude, CorrectedLongitude,
        site_targets$lat[i], site_targets$lon[i]
      )
    ) %>%
    filter(dist_m <= radius_m)
  
  if (nrow(candidates) < site_targets$n_target[i]) {
    stop(
      paste0(
        "Not enough candidates for ", site_targets$site[i],
        ": need ", site_targets$n_target[i],
        ", found ", nrow(candidates)
      )
    )
  }
  
  picked <- candidates %>%
    slice_sample(n = site_targets$n_target[i]) %>%
    mutate(location_pattern = site_targets$site[i])
  
  site_samples_list[[i]] <- picked
  used_ids <- c(used_ids, picked$PlantID)
}

wild_site_sampled <- bind_rows(site_samples_list)

wild_plot_samples <- bind_rows(wild_named_sampled, wild_site_sampled)

cat("Final wild plot sample size:", nrow(wild_plot_samples), "\n")
print(table(wild_plot_samples$location_pattern))
stopifnot(all(wild_plot_samples$domesticated == "wild"))
stopifnot(!anyDuplicated(wild_plot_samples$PlantID))

# -------------------------
# 4. Marker order and annotation
# -------------------------
marker_info_plot <- marker_info_filtered %>%
  mutate(
    LG_num = suppressWarnings(as.numeric(gsub("[^0-9.]", "", LG))),
    cM_num = suppressWarnings(as.numeric(cM))
  ) %>%
  filter(LG_num != 10) %>%
  arrange(LG_num, cM_num)

geno_cols <- marker_info_plot$LocusName
geno_cols <- geno_cols[geno_cols %in% names(data_filtered)]

marker_info_plot <- marker_info_plot %>%
  filter(LocusName %in% geno_cols) %>%
  mutate(
    marker = factor(LocusName, levels = geno_cols),
    marker_index = seq_along(geno_cols),
    gene_strip = ifelse(effect == "colour" & !is.na(gene_id) & gene_id != "", gene_id, "none")
  )

chr_breaks <- marker_info_plot %>%
  group_by(LG) %>%
  summarise(
    start = min(marker_index),
    end = max(marker_index),
    .groups = "drop"
  )

vline_positions <- chr_breaks$end[-nrow(chr_breaks)] + 0.5

# -------------------------
# 5. Build domesticated dataset
# -------------------------
dom_data <- data_filtered %>%
  filter(domesticated == "domesticated", !is.na(line)) %>%
  select(line, replicate_number, species, all_of(geno_cols))

# -------------------------
# 6. Polarise genotypes using delta p between red and yellow wild populations
#    so genotype 2 = allele more common in red than yellow
# -------------------------
red_subset <- wild_plot_samples %>%
  filter(location_pattern %in% red_pops)

yellow_subset <- wild_plot_samples %>%
  filter(location_pattern %in% yellow_pops)

geno_red <- as.matrix(red_subset[, geno_cols, drop = FALSE])
geno_yellow <- as.matrix(yellow_subset[, geno_cols, drop = FALSE])

geno_red[geno_red %in% c(-9, -10)] <- NA
geno_yellow[geno_yellow %in% c(-9, -10)] <- NA

p_red <- colMeans(geno_red / 2, na.rm = TRUE)
p_yellow <- colMeans(geno_yellow / 2, na.rm = TRUE)
delta_p <- p_red - p_yellow

flip_loci <- names(delta_p)[!is.na(delta_p) & delta_p < 0]

cat("Number of loci flipped:", length(flip_loci), "\n")

polarisation_summary <- data.frame(
  marker = names(delta_p),
  p_red = as.numeric(p_red),
  p_yellow = as.numeric(p_yellow),
  delta_p = as.numeric(delta_p),
  flipped = names(delta_p) %in% flip_loci
)

wild_plot_samples <- flip_012_cols(wild_plot_samples, flip_loci)
dom_data <- flip_012_cols(dom_data, flip_loci)

# optional sanity check
red_subset_check <- wild_plot_samples %>%
  filter(location_pattern %in% red_pops)

yellow_subset_check <- wild_plot_samples %>%
  filter(location_pattern %in% yellow_pops)

geno_red_check <- as.matrix(red_subset_check[, geno_cols, drop = FALSE])
geno_yellow_check <- as.matrix(yellow_subset_check[, geno_cols, drop = FALSE])

geno_red_check[geno_red_check %in% c(-9, -10)] <- NA
geno_yellow_check[geno_yellow_check %in% c(-9, -10)] <- NA

delta_p_check <- colMeans(geno_red_check / 2, na.rm = TRUE) -
  colMeans(geno_yellow_check / 2, na.rm = TRUE)

cat("Minimum delta p after flipping:", min(delta_p_check, na.rm = TRUE), "\n")

# -------------------------
# 7. Wild plotting data
# -------------------------
plot_df_wild <- wild_plot_samples %>%
  mutate(
    location_pattern = factor(location_pattern, levels = custom_order_wild)
  ) %>%
  arrange(location_pattern, PlantID) %>%
  select(PlantID, location_pattern, all_of(geno_cols)) %>%
  mutate(individual_id = PlantID)

plot_df_wild$individual_id <- factor(
  plot_df_wild$individual_id,
  levels = rev(plot_df_wild$individual_id)
)

geno_long_wild <- plot_df_wild %>%
  pivot_longer(
    cols = all_of(geno_cols),
    names_to = "marker",
    values_to = "genotype"
  )

geno_long_wild$marker <- factor(geno_long_wild$marker, levels = geno_cols)
geno_long_wild$genotype[geno_long_wild$genotype %in% c(-9, -10)] <- NA
geno_long_wild$genotype <- factor(geno_long_wild$genotype, levels = c(0, 1, 2))

group_runs_wild <- rle(as.character(plot_df_wild$location_pattern))
hline_positions_wild <- cumsum(group_runs_wild$lengths) + 0.5
hline_positions_wild <- hline_positions_wild[-length(hline_positions_wild)]

# -------------------------
# 8. Domesticated plotting data
# -------------------------
plot_df_dom <- dom_data %>%
  mutate(
    species = factor(
      species,
      levels = c("majus", "majus_ssp_tortuosum", "molle", "siculum")
    ),
    line_num = suppressWarnings(as.numeric(sub("^L", "", line)))
  ) %>%
  arrange(species, line_num, suppressWarnings(as.numeric(replicate_number))) %>%
  mutate(individual_id = paste(line, replicate_number, sep = "_")) %>%
  select(line, individual_id, all_of(geno_cols))

plot_df_dom$individual_id <- factor(
  plot_df_dom$individual_id,
  levels = rev(plot_df_dom$individual_id)
)

geno_long_dom <- plot_df_dom %>%
  pivot_longer(
    cols = all_of(geno_cols),
    names_to = "marker",
    values_to = "genotype"
  )

geno_long_dom$marker <- factor(geno_long_dom$marker, levels = geno_cols)
geno_long_dom$genotype[geno_long_dom$genotype %in% c(-9, -10)] <- NA
geno_long_dom$genotype <- factor(geno_long_dom$genotype, levels = c(0, 1, 2))

line_sizes <- plot_df_dom %>%
  count(line)

hline_positions_dom <- cumsum(line_sizes$n) + 0.5
hline_positions_dom <- hline_positions_dom[-length(hline_positions_dom)]

# -------------------------
# 9. Annotation strip
# -------------------------
annot_strip <- marker_info_plot %>%
  transmute(
    marker,
    strip_y = "Gene strip",
    gene_strip
  )

# -------------------------
# 10. Wild plot
# -------------------------
p_wild <- ggplot() +
  geom_tile(
    data = geno_long_wild,
    aes(x = marker, y = individual_id, fill = genotype)
  ) +
  geom_tile(
    data = annot_strip,
    aes(x = marker, y = strip_y, fill = gene_strip),
    height = 4.4
  ) +
  geom_vline(xintercept = vline_positions, linewidth = 0, colour = "white") +
  geom_hline(yintercept = hline_positions_wild, linewidth = 0.6, colour = "white") +
  scale_fill_manual(
    values = c(
      "0" = "steelblue",
      "1" = "gold",
      "2" = "firebrick",
      "none" = "grey85",
      "cre" = "#F8766D",
      "def" = "#D89000",
      "dich" = "#A3A500",
      "el" = "#39B600",
      "fla" = "#00BF7D",
      "mixta" = "#00BFC4",
      "ros" = "#619CFF",
      "sulf" = "#DB72FB",
      "ven" = "#FF61C3"
    ),
    breaks = c("0", "1", "2", "none", "cre", "def", "dich", "el", "fla", "mixta", "ros", "sulf", "ven"),
    labels = c("0", "1", "2", "No gene", "cre", "def", "dich", "el", "fla", "mixta", "ros", "sulf", "ven")
  ) +
  theme_classic() +
  labs(
    x = "Marker",
    y = "Wild individuals",
    fill = "Annotation / genotype"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.3, "cm")
  )

# -------------------------
# 11. Domesticated plot
# -------------------------
p_dom <- ggplot() +
  geom_tile(
    data = geno_long_dom,
    aes(x = marker, y = individual_id, fill = genotype)
  ) +
  geom_tile(
    data = annot_strip,
    aes(x = marker, y = strip_y, fill = gene_strip),
    height = 4.4
  ) +
  geom_vline(xintercept = vline_positions, linewidth = 0, colour = "white") +
  geom_hline(yintercept = hline_positions_dom, linewidth = 0.4, colour = "white") +
  scale_fill_manual(
    values = c(
      "0" = "steelblue",
      "1" = "gold",
      "2" = "firebrick",
      "none" = "grey85",
      "cre" = "#F8766D",
      "def" = "#D89000",
      "dich" = "#A3A500",
      "el" = "#39B600",
      "fla" = "#00BF7D",
      "mixta" = "#00BFC4",
      "ros" = "#619CFF",
      "sulf" = "#DB72FB",
      "ven" = "#FF61C3"
    ),
    breaks = c("0", "1", "2", "none", "cre", "def", "dich", "el", "fla", "mixta", "ros", "sulf", "ven"),
    labels = c("0", "1", "2", "No gene", "cre", "def", "dich", "el", "fla", "mixta", "ros", "sulf", "ven")
  ) +
  theme_classic() +
  labs(
    x = "Marker",
    y = "Domesticated individuals",
    fill = "Annotation / genotype"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.3, "cm")
  )

# -------------------------
# 12. Combine on one page
# -------------------------
p_wild / p_dom
