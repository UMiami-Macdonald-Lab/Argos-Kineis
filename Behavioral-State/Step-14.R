# ============================================================
# CLUSTER COMPOSITION FOR EACH SHARK AND TRACKING SYSTEM
# ============================================================

library(dplyr)
library(tidyr)

cluster_names <- paste0(
  "cluster_",
  seq_len(max_clusters)
)

stopifnot(
  nrow(theta) == nrow(segment_bin_counts),
  ncol(theta) == length(cluster_names)
)

colnames(theta) <- cluster_names

# Attach the segment keys to their fitted proportions.
segment_proportions <- bind_cols(
  as_tibble(
    segment_bin_counts[, c("id", "tseg")]
  ),
  as_tibble(theta)
)

stopifnot(
  anyDuplicated(
    segment_proportions[c("id", "tseg")]
  ) == 0L
)

# Repeat each segment's composition over its movement rows.
# These remain segment-level estimates, not independently
# estimated probabilities for individual fixes.
movement_with_clusters <- tracks_segmented %>%
  left_join(
    segment_proportions,
    by = c("id", "tseg")
  )

proportion_matrix <- as.matrix(
  movement_with_clusters[, cluster_names]
)

stopifnot(
  nrow(movement_with_clusters) == nrow(tracks_segmented),
  all(is.finite(proportion_matrix)),
  all(proportion_matrix >= 0),
  all(abs(rowSums(proportion_matrix) - 1) < 1e-5)
)

# Average over daily movement rows within each shark/system.
# Longer segments therefore contribute more than short ones.
individual_budgets <- movement_with_clusters %>%
  group_by(track_id, system) %>%
  summarise(
    n_movement_rows = n(),
    
    across(
      all_of(cluster_names),
      ~ mean(.x)
    ),
    
    .groups = "drop"
  )

# Full results: all seven clusters retained.
print(
  individual_budgets,
  n = Inf,
  width = Inf
)

# Compact display of the two major clusters and the remainder.
major_budget_summary <- individual_budgets %>%
  transmute(
    track_id,
    system,
    n_movement_rows,
    cluster_1_percent = 100 * cluster_1,
    cluster_2_percent = 100 * cluster_2,
    
    remaining_clusters_percent =
      100 * (1 - cluster_1 - cluster_2)
  )

print(
  major_budget_summary,
  n = Inf,
  width = Inf
)