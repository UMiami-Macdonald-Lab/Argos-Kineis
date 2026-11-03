# ============================================================
# EXTRACT CLUSTER CONTRIBUTIONS
# ============================================================

library(bayesmove)
library(dplyr)
library(tibble)

# Posterior mean cluster proportions for each segment.
theta <- bayesmove::extract_prop(
  res = clustering_fit,
  ngibbs = lda_iterations,
  nburn = lda_burn,
  nmaxclust = max_clusters
)

theta <- as.matrix(theta)

stopifnot(
  nrow(theta) == nrow(segment_bin_counts),
  ncol(theta) == max_clusters,
  all(is.finite(theta)),
  all(theta >= 0),
  all(abs(rowSums(theta) - 1) < 1e-5)
)

# Preserve the model's original cluster numbers.
colnames(theta) <- paste0(
  "cluster_",
  seq_len(max_clusters)
)

cluster_contributions <- tibble(
  cluster = colnames(theta),
  mean_segment_prop = colMeans(theta)
) %>%
  arrange(desc(mean_segment_prop)) %>%
  mutate(
    cumulative_prop = cumsum(mean_segment_prop),
    
    # Include the cluster that crosses the 90% threshold.
    # Use full-precision values for this calculation.
    retained_over_90_rule =
      lag(cumulative_prop, default = 0) <= 0.90
  )

print(
  cluster_contributions,
  n = Inf,
  width = Inf
)

cat(
  "\nNumber of clusters retained by the >90% rule:",
  sum(cluster_contributions$retained_over_90_rule),
  "\n"
)