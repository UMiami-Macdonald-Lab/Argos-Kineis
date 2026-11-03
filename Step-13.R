# ============================================================
# INSPECT FITTED MOVEMENT DISTRIBUTIONS FOR EACH CLUSTER
# ============================================================

library(bayesmove)
library(dplyr)
library(ggplot2)

# Extract fitted probabilities for each movement-variable bin.
fitted_bins <- bayesmove::get_behav_hist(
  dat = clustering_fit,
  nburn = lda_burn,
  ngibbs = lda_iterations,
  nmaxclust = max_clusters,
  var.names = c(
    "Step length",
    "Turning angle"
  )
)

# Preserve the original cluster numbering.
fitted_bins <- fitted_bins %>%
  mutate(
    cluster = factor(
      as.character(behav),
      levels = as.character(seq_len(max_clusters)),
      labels = paste("Cluster", seq_len(max_clusters))
    )
  )

stopifnot(!anyNA(fitted_bins$cluster))

cluster_distribution_plot <- ggplot(
  fitted_bins,
  aes(
    x = bin,
    y = prop
  )
) +
  geom_col(
    fill = "#3B8C95",
    width = 0.8
  ) +
  facet_grid(
    cluster ~ var,
    scales = "free_x"
  ) +
  scale_x_continuous(
    breaks = 1:8
  ) +
  labs(
    x = "Movement-variable bin",
    y = "Fitted bin probability"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

print(cluster_distribution_plot)

# Reprint the boundaries to help interpret the plot.
print(step_bin_table, n = Inf)
print(angle_bin_table, n = Inf)