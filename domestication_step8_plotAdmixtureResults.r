######################## admixture plotting script -- assumes output of admiture is in working directory

library(dplyr)
library(ggplot2)
library(tidyr)

# -------------------------
# 1. user settings
# -------------------------

prefix <- "snapdragon_admixture"
K <- 3                             ###################################### change value of K here

custom_order_wild <- c(
  "BAR****", "BES****", "CAL****", "CLA****", "D_***_R*", "GARL****",
  "GCHI****", "GPER****", "GPUI****", "M406*", "MP113*", "MRT****",
  "PLAMN", "PLAMF", "SUE****", "VIE****", "D27**", "PLAHZ", "CADI****",
  "ELS****", "GBOU****", "GCAM****", "GDOS****", "GFAB****", "LU****",
  "MLYS****", "PAL****", "POM****", "YP044*", "PLAYF", "PLAYN",
  "USS****", "VIL****"
)

# -------------------------
# 2. read ADMIXTURE output
# -------------------------

Q_file <- paste0(prefix, ".", K, ".Q")
fam_file <- paste0(prefix, ".fam")

Q <- read.table(Q_file, header = FALSE)
colnames(Q) <- paste0("Cluster", seq_len(ncol(Q)))

fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

admix_df <- bind_cols(fam, Q)

# -------------------------
# 3. rebuild metadata from original objects
# -------------------------

wild_meta <- wild_plot_samples %>%
  transmute(
    IID = PlantID,
    group = "wild",
    unit = as.character(location_pattern),
    unit_order = match(unit, custom_order_wild)
  )

dom_meta <- data_filtered %>%
  filter(
    domesticated == "domesticated",
    !is.na(line),
    species == "majus"
  ) %>%
  transmute(
    IID = paste(line, replicate_number, sep = "_"),
    group = "domesticated",
    unit = as.character(line),
    unit_order = suppressWarnings(as.numeric(sub("^L", "", unit)))
  )

sample_meta <- bind_rows(wild_meta, dom_meta)

admix_df <- admix_df %>%
  left_join(sample_meta, by = "IID") %>%
  filter(!is.na(group))

cat("Matched samples:", nrow(admix_df), "\n")

# -------------------------
# 4. order samples exactly once
# -------------------------

cluster_cols <- paste0("Cluster", seq_len(K))

admix_df <- admix_df %>%
  mutate(
    block_order = ifelse(group == "wild", 1L, 2L)
  ) %>%
  arrange(
    block_order,
    unit_order,
    unit,
    dplyr::desc(.data[[cluster_cols[1]]]),
    if (K >= 2) dplyr::desc(.data[[cluster_cols[2]]]) else row_number()
  )

# for K > 2, refine ordering with remaining clusters
if (K > 2) {
  for (cl in cluster_cols[3:K]) {
    admix_df <- admix_df %>% arrange(block_order, unit_order, unit, dplyr::desc(.data[[cluster_cols[1]]]),
                                     dplyr::desc(.data[[cluster_cols[2]]]), dplyr::desc(.data[[cl]]))
  }
}

admix_df <- admix_df %>%
  mutate(plot_pos = seq_len(n()))

# sanity check
print(head(admix_df[, c("IID", "group", "unit", "unit_order", "plot_pos")]))

# -------------------------
# 5. build unit positions from the same ordered table
# -------------------------

unit_sizes <- admix_df %>%
  group_by(block_order, group, unit, unit_order) %>%
  summarise(
    start = min(plot_pos),
    end = max(plot_pos),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(block_order, unit_order, unit) %>%
  mutate(
    mid = (start + end) / 2
  )

separator_positions <- unit_sizes$end[-nrow(unit_sizes)] + 0.5

block_separator <- admix_df %>%
  group_by(block_order, group) %>%
  summarise(end = max(plot_pos), .groups = "drop") %>%
  arrange(block_order)

block_separator <- if (nrow(block_separator) > 1) block_separator$end[1] + 0.5 else numeric(0)

# -------------------------
# 6. long format for plotting
# -------------------------

plot_df <- admix_df[, c("IID", "unit", "plot_pos", cluster_cols), drop = FALSE] %>%
  pivot_longer(
    cols = all_of(cluster_cols),
    names_to = "cluster",
    values_to = "q"
  )

plot_df$cluster <- factor(plot_df$cluster, levels = cluster_cols)

# sanity check
print(head(plot_df))

# -------------------------
# 7. colours
# -------------------------

cluster_palette <- setNames(
  grDevices::hcl.colors(K, palette = "Set 2"),
  cluster_cols
)

# -------------------------
# 8. plot
# -------------------------

ggplot(plot_df, aes(x = plot_pos, y = q, fill = cluster)) +
  geom_col(width = 1) +
  geom_vline(
    xintercept = separator_positions,
    colour = "white",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = block_separator,
    colour = "white",
    linewidth = 1.2
  ) +
  geom_text(
    data = unit_sizes,
    aes(x = mid, y = 1.02, label = unit),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 0,
    vjust = 0.5,
    size = 2.6
  ) +
  scale_fill_manual(values = cluster_palette, drop = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.08), clip = "off") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(10, 10, 30, 10),
    legend.position = "right"
  ) +
  labs(
    x = "Individuals",
    y = "Ancestry proportion",
    fill = "Cluster",
    title = paste("ADMIXTURE results (K =", K, ")")
  )