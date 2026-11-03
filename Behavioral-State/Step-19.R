# ============================================================
# LONGER CLUSTERING RUNS TO CHECK STABILITY
# ============================================================

library(bayesmove)

long_iterations <- 50000L
long_burn <- 25000L

# Keep the same data, priors, and maximum cluster count.
fit_long_chain <- function(seed_value) {
  
  set.seed(seed_value)
  
  bayesmove::cluster_segments(
    dat = segment_bin_counts,
    gamma1 = lda_gamma,
    alpha = lda_alpha,
    ngibbs = long_iterations,
    nburn = long_burn,
    nmaxclust = max_clusters,
    ndata.types = 2L
  )
}

# Two separate runs.
clustering_long_a <- fit_long_chain(
  seed_value = 20261104L
)

clustering_long_b <- fit_long_chain(
  seed_value = 20261105L
)

# Save both fits without replacing the earlier checkpoints.
saveRDS(
  list(
    run_a = clustering_long_a,
    run_b = clustering_long_b,
    seeds = c(20261104L, 20261105L),
    iterations = long_iterations,
    burn_in = long_burn,
    max_clusters = max_clusters,
    alpha = lda_alpha,
    gamma1 = lda_gamma
  ),
  file = "Bayesmove_long_clustering_checkpoint.rds"
)

# ------------------------------------------------------------
# COMPARE COMPLETE LIKELIHOOD TRACES
# ------------------------------------------------------------

shared_ylim <- range(
  c(
    clustering_long_a$loglikel,
    clustering_long_b$loglikel
  ),
  finite = TRUE
)

old_par <- par(
  mfrow = c(1, 2)
)

plot(
  clustering_long_a$loglikel,
  type = "l",
  ylim = shared_ylim,
  xlab = "Iteration",
  ylab = "Log likelihood",
  main = "Long run A"
)

abline(
  v = long_burn,
  col = "red",
  lty = 2
)

plot(
  clustering_long_b$loglikel,
  type = "l",
  ylim = shared_ylim,
  xlab = "Iteration",
  ylab = "Log likelihood",
  main = "Long run B"
)

abline(
  v = long_burn,
  col = "red",
  lty = 2
)

par(old_par)