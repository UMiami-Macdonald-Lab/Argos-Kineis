# ============================================================
# CLUSTERING STABILITY: SECOND INDEPENDENT RUN
# ============================================================

library(bayesmove)

# Different seed; otherwise use the same settings.
lda_seed_repeat <- 20261105L

set.seed(lda_seed_repeat)

clustering_fit_repeat <- bayesmove::cluster_segments(
  dat = segment_bin_counts,
  
  gamma1 = lda_gamma,
  alpha = lda_alpha,
  
  ngibbs = lda_iterations,
  nburn = lda_burn,
  
  nmaxclust = max_clusters,
  ndata.types = 2L
)

# Preserve the second run separately.
saveRDS(
  list(
    fit = clustering_fit_repeat,
    seed = lda_seed_repeat,
    iterations = lda_iterations,
    burn_in = lda_burn,
    max_clusters = max_clusters
  ),
  file = "Bayesmove_clustering_repeat_checkpoint.rds"
)

# Inspect its likelihood trace.
plot(
  clustering_fit_repeat$loglikel,
  type = "l",
  xlab = "Iteration",
  ylab = "Log likelihood",
  main = "Joint clustering: second independent run"
)

abline(
  v = lda_burn,
  col = "red",
  lty = 2
)

# Confirm the repeat fit exists.
stopifnot(exists("clustering_fit_repeat"))

# Check whether the stored likelihood traces are identical.
cat(
  "Are the two likelihood traces exactly identical? ",
  identical(
    clustering_fit$loglikel,
    clustering_fit_repeat$loglikel
  ),
  "\n"
)

# Explicitly plot the second run.
plot(
  clustering_fit_repeat$loglikel,
  type = "l",
  xlab = "Iteration",
  ylab = "Log likelihood",
  main = "SECOND RUN — seed 20261105"
)

abline(
  v = lda_burn,
  col = "red",
  lty = 2
)