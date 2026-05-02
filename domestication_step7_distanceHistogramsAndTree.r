################################### plot histograms of pariwise distances and tree
######## run relevant sections seprately

# ---------------------------------------------------------------------------
# 9. histograms
# ---------------------------------------------------------------------------

ggplot(plot_df_dip, aes(x = distance)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  geom_vline(
    data = distance_summary_dip,
    aes(xintercept = mean_distance),
    colour = "blue",
    linewidth = 0.8
  ) +
  facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Pairwise diploid allele-level distance",
    y = "Count"
  )





# ---------------------------------------------------------------------------
# 10. build NJ tree
# ---------------------------------------------------------------------------

dist_obj_dip <- as.dist(dist_mat_dip)
tree_nj_dip <- nj(dist_obj_dip)

# optional: try bionj instead
# tree_nj_dip <- bionj(dist_obj_dip)

tree_plot <- ladderize(tree_nj_dip, right = TRUE)

# fix negative branch lengths
n_neg <- sum(tree_plot$edge.length < 0, na.rm = TRUE)
cat("Number of negative branches:", n_neg, "\n")
tree_plot$edge.length[tree_plot$edge.length < 0] <- 0

# -------------------------
# 11. aesthetics
# -------------------------

cb_palette <- c(
  "gray50",
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#0072B2",
  "#D55E00",
  "#F0E442"
)

wild_cols <- list(
  red = "magenta",
  yellow = "yellow",
  hybrid = "orange",
  other = "grey70"
)

wild_shapes <- list(
  red = 21,
  yellow = 22,
  hybrid = 24,
  other = 21
)

dom_lines <- unique(all_dip$unit[all_dip$group == "domesticated"])
dom_line_num <- suppressWarnings(as.numeric(sub("^L", "", dom_lines)))
dom_lines <- dom_lines[order(dom_line_num)]

shape_pool <- c(21, 22, 23, 24, 25)

dom_line_shape <- setNames(
  rep(shape_pool, length.out = length(dom_lines)),
  dom_lines
)

dom_line_fill <- setNames(
  rep(cb_palette, length.out = length(dom_lines)),
  dom_lines
)

# -------------------------
# 12. build tip metadata
# -------------------------

tip_meta <- all_dip %>%
  select(sample_id, group, unit) %>%
  mutate(
    category = case_when(
      group == "domesticated" ~ "domesticated",
      unit %in% red_pops ~ "wild_red",
      unit %in% yellow_pops ~ "wild_yellow",
      unit %in% hybrid_pops ~ "wild_hybrid",
      TRUE ~ "wild_other"
    ),
    shape_code = case_when(
      group == "domesticated" ~ unname(dom_line_shape[unit]),
      category == "wild_red" ~ wild_shapes$red,
      category == "wild_yellow" ~ wild_shapes$yellow,
      category == "wild_hybrid" ~ wild_shapes$hybrid,
      TRUE ~ wild_shapes$other
    ),
    fill_code = case_when(
      group == "domesticated" ~ unname(dom_line_fill[unit]),
      category == "wild_red" ~ wild_cols$red,
      category == "wild_yellow" ~ wild_cols$yellow,
      category == "wild_hybrid" ~ wild_cols$hybrid,
      TRUE ~ wild_cols$other
    )
  )

# -------------------------
# 13. plot with ggtree
# -------------------------

p <- ggtree(tree_plot, layout = "equal_angle", open.angle = 200, color = "grey75", size = 0.25) %<+% tip_meta +
  geom_tippoint(
    aes(fill = fill_code, shape = factor(shape_code)),
    size = 2.2,
    stroke = 0.35,
    colour = "black"
  ) +
  scale_shape_manual(
    values = sort(unique(tip_meta$shape_code)),
    breaks = as.character(sort(unique(tip_meta$shape_code))),
    guide = "none"
  ) +
  scale_fill_identity() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Neighbour-joining tree of diploid individuals")

print(p)

# -------------------------
# 14. quick checks
# -------------------------

cat("Number of loci used:", length(non_colour_loci), "\n")
cat("Number of individuals:", nrow(all_dip), "\n")
cat("Number of domesticated lines:", length(dom_lines), "\n")
