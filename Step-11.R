# ============================================================
# BAYESMOVE STAGE 2: JOINTLY CLUSTER THE TRACK SEGMENTS
# ============================================================

library(bayesmove)
library(dplyr)

# ------------------------------------------------------------
# 1. SUMMARIZE BIN COUNTS WITHIN EACH SEGMENT
# ------------------------------------------------------------

segment_input <- tracks_segmented %>%
  select(
    id,
    tseg,
    SL,
    TA
  ) %>%
  as.data.frame()

segment_bin_counts <- bayesmove::summarize_tsegs(
  dat = segment_input,
  nbins = c(5L, 8L)
)

# Check that both movement variables account for every row.
sl_columns <- grep(
  "^y1\\.",
  names(segment_bin_counts),
  value = TRUE
)

ta_columns <- grep(
  "^y2\\.",
  names(segment_bin_counts),
  value = TRUE
)

stopifnot(
  length(sl_columns) == 5L,
  length(ta_columns) == 8L,
  nrow(segment_bin_counts) == nrow(segment_lengths),
  
  all(
    rowSums(segment_bin_counts[, sl_columns, drop = FALSE]) ==
      rowSums(segment_bin_counts[, ta_columns, drop = FALSE])
  ),
  
  sum(
    as.matrix(
      segment_bin_counts[, sl_columns, drop = FALSE]
    )
  ) == nrow(tracks_segmented)
)

cat(
  "Segments entering clustering:",
  nrow(segment_bin_counts),
  "\n"
)

# ------------------------------------------------------------
# 2. SET CLUSTERING PARAMETERS
# ------------------------------------------------------------

lda_seed <- 20261104L

lda_iterations <- 10000L
lda_burn <- 5000L

max_clusters <- 7L

lda_alpha <- 0.1
lda_gamma <- 0.1

# ------------------------------------------------------------
# 3. FIT ONE JOINT CLUSTERING MODEL
# ------------------------------------------------------------

set.seed(lda_seed)

clustering_fit <- bayesmove::cluster_segments(
  dat = segment_bin_counts,
  
  gamma1 = lda_gamma,
  alpha = lda_alpha,
  
  ngibbs = lda_iterations,
  nburn = lda_burn,
  
  nmaxclust = max_clusters,
  
  # Two data streams: step length and turning angle.
  ndata.types = 2L
)

# Save a checkpoint.
saveRDS(
  list(
    fit = clustering_fit,
    segment_bin_counts = segment_bin_counts,
    tracks_segmented = tracks_segmented,
    seed = lda_seed,
    iterations = lda_iterations,
    burn_in = lda_burn,
    max_clusters = max_clusters,
    alpha = lda_alpha,
    gamma1 = lda_gamma,
    step_limits = step_limits,
    angle_limits = angle_limits
  ),
  file = "Bayesmove_clustering_checkpoint.rds"
)

# ------------------------------------------------------------
# 4. INSPECT THE CLUSTERING LIKELIHOOD TRACE
# ------------------------------------------------------------

plot(
  clustering_fit$loglikel,
  type = "l",
  xlab = "Iteration",
  ylab = "Log likelihood",
  main = "Joint segment-clustering model"
)

abline(
  v = lda_burn,
  col = "red",
  lty = 2
)