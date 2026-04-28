############################### Genotype phenotype association plots for each gene

library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)

# -------------------------
# 1. Marker/phenotype combinations
# -------------------------
plot_specs <- tibble::tribble(
  ~marker,                ~phenotype,
  "ros_assembly_567004",  "Red",
  "ros_assembly_715001",  "Red",
  "s1187_290152",         "Yellow",
  "s316_93292",           "Yellow",
  "s91_181717",           "Yellow"
)

set.seed(123)

# -------------------------
# 2. Define samples once
#    wild sample size matched to domesticated majus
# -------------------------
dom_data <- data_filtered %>%
  filter(
    domesticated == "domesticated",
    species == "majus"
  )

n_dom <- nrow(dom_data)

wild_data <- data_filtered %>%
  filter(domesticated == "wild") %>%
  slice_sample(n = n_dom)

base_plot_data <- bind_rows(
  wild_data %>% mutate(sample_type = "wild"),
  dom_data %>% mutate(sample_type = "domesticated_majus")
) %>%
  mutate(
    sample_type = factor(
      sample_type,
      levels = c("wild", "domesticated_majus")
    )
  )

# -------------------------
# 3. Plotting function
# -------------------------
make_marker_plot <- function(marker_to_plot, phenotype_to_plot, df) {
  plot_data <- df %>%
    filter(!is.na(.data[[phenotype_to_plot]])) %>%
    mutate(
      genotype = .data[[marker_to_plot]],
      genotype = ifelse(genotype %in% c(-9, -10), NA, genotype),
      genotype = ifelse(genotype == 0, 2,
                 ifelse(genotype == 2, 0, genotype)),
      genotype = factor(genotype, levels = c(0, 1, 2))
    ) %>%
    filter(!is.na(genotype))
  mean_data <- plot_data %>%
    group_by(sample_type, genotype) %>%
    summarise(
      mean_value = mean(.data[[phenotype_to_plot]], na.rm = TRUE),
      .groups = "drop"
    )
  ggplot(plot_data, aes(x = genotype, y = .data[[phenotype_to_plot]])) +
    ggbeeswarm::geom_quasirandom(
      width = 0.18,
      alpha = 0.75,
      size = 1.2,
      colour = "grey35"
    ) +
    geom_line(
      data = mean_data,
      aes(x = genotype, y = mean_value, group = 1),
      inherit.aes = FALSE,
      linewidth = 0.8,
      colour = "red"
    ) +
    geom_point(
      data = mean_data,
      aes(x = genotype, y = mean_value),
      inherit.aes = FALSE,
      size = 2,
      colour = "red"
    ) +
    facet_wrap(~ sample_type, nrow = 1) +
    theme_classic() +
    labs(
      x = "",
      y = phenotype_to_plot
    ) +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      strip.text = element_blank(),
      strip.background = element_blank()
    )
}

# -------------------------
# 4. Build all five plots
# -------------------------
plot_list <- lapply(seq_len(nrow(plot_specs)), function(i) {
  make_marker_plot(
    marker_to_plot = plot_specs$marker[i],
    phenotype_to_plot = plot_specs$phenotype[i],
    df = base_plot_data
  )
})

# -------------------------
# 5. Arrange in one column
# -------------------------
wrap_plots(plot_list, ncol = 1)



############# spearman correlations for the same marker/phenotype combinations

library(dplyr)

# -------------------------
# 1. Marker/phenotype combinations
# -------------------------
test_specs <- tibble::tribble(
  ~marker,                ~phenotype,
  "ros_assembly_567004",  "Red",
  "ros_assembly_715001",  "Red",
  "s1187_290152",         "Yellow",
  "s316_93292",           "Yellow",
  "s91_181717",           "Yellow"
)

set.seed(123)

# -------------------------
# 2. Define samples once
#    wild sample size matched to domesticated majus
# -------------------------
dom_data <- data_filtered %>%
  filter(
    domesticated == "domesticated",
    species == "majus"
  )

n_dom <- nrow(dom_data)

wild_data <- data_filtered %>%
  filter(domesticated == "wild") %>%
  slice_sample(n = n_dom)

base_plot_data <- bind_rows(
  wild_data %>% mutate(sample_type = "wild"),
  dom_data %>% mutate(sample_type = "domesticated_majus")
) %>%
  mutate(
    sample_type = factor(
      sample_type,
      levels = c("wild", "domesticated_majus")
    )
  )

# -------------------------
# 3. Function to calculate Spearman correlation
# -------------------------
calc_spearman_one <- function(df, marker_to_test, phenotype_to_test) {
  dat <- df %>%
    filter(!is.na(.data[[phenotype_to_test]])) %>%
    mutate(
      genotype = .data[[marker_to_test]],
      genotype = ifelse(genotype %in% c(-9, -10), NA, genotype),
      # flip homozygotes to match your plotting code
      genotype = ifelse(genotype == 0, 2,
                 ifelse(genotype == 2, 0, genotype))
    ) %>%
    filter(!is.na(genotype))
  x <- dat$genotype
  y <- dat[[phenotype_to_test]]
  if (length(unique(x)) < 2 || length(unique(y)) < 2 || nrow(dat) < 3) {
    return(data.frame(
      n = nrow(dat),
      rho = NA_real_,
      p_value = NA_real_
    ))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(
    n = nrow(dat),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}

# -------------------------
# 4. Run across all markers and both sample types
# -------------------------
results_list <- vector("list", nrow(test_specs) * 2)
k <- 1

for (i in seq_len(nrow(test_specs))) {
  this_marker <- test_specs$marker[i]
  this_pheno  <- test_specs$phenotype[i]
  for (this_type in levels(base_plot_data$sample_type)) {
    this_df <- base_plot_data %>%
      filter(sample_type == this_type)
    res <- calc_spearman_one(
      df = this_df,
      marker_to_test = this_marker,
      phenotype_to_test = this_pheno
    )
    results_list[[k]] <- data.frame(
      sample_type = this_type,
      marker = this_marker,
      phenotype = this_pheno,
      res
    )
    k <- k + 1
  }
}

correlation_results <- bind_rows(results_list)

print(correlation_results)

# optional save
write.csv(
  correlation_results,
  "spearman_correlations_selected_markers.csv",
  row.names = FALSE
)
