# ============================================================
# COMPARE CLUSTER DISTRIBUTIONS AND CONTRIBUTIONS BETWEEN RUNS
# ============================================================

library(bayesmove)
library(dplyr)
library(tibble)
library(ggplot2)

# Number of analyzed daily movement rows in each segment.
sl_columns <- grep(
  "^y1\\.",
  names(segment_bin_counts),
  value = TRUE
)

segment_weights <- rowSums(
  as.matrix(
    segment_bin_counts[, sl_columns, drop = FALSE]
  )
)

stopifnot(
  all(segment_weights > 0),
  sum(segment_weights) == nrow(tracks_segmented)
)

summarize_long_run <- function(fit, run_label) {
  
  theta_run <- bayesmove::extract_prop(
    res = fit,
    ngibbs = long_iterations,
    nburn = long_burn,
    nmaxclust = max_clusters
  ) %>%
    as.matrix()
  
  stopifnot(
    nrow(theta_run) == length(segment_weights),
    ncol(theta_run) == max_clusters,
    all(is.finite(theta_run)),
    all(theta_run >= 0),
    all(abs(rowSums(theta_run) - 1) < 1e-5)
  )
  
  # Weight segment compositions by their movement-row counts.
  weighted_theta <- sweep(
    theta_run,
    MARGIN = 1,
    STATS = segment_weights,
    FUN = "*"
  )
  
  contributions <- tibble(
    run = run_label,
    
    cluster = paste0(
      "cluster_",
      seq_len(max_clusters)
    ),
    
    mean_segment_prop = colMeans(theta_run),
    
    mean_daily_row_prop =
      colSums(weighted_theta) / sum(segment_weights)
  ) %>%
    arrange(desc(mean_segment_prop))
  
  bin_distributions <- bayesmove::get_behav_hist(
    dat = fit,
    nburn = long_burn,
    ngibbs = long_iterations,
    nmaxclust = max_clusters,
    var.names = c(
      "Step length",
      "Turning angle"
    )
  ) %>%
    mutate(
      run = run_label,
      cluster = factor(
        as.character(behav),
        levels = as.character(seq_len(max_clusters)),
        labels = paste(
          "Cluster",
          seq_len(max_clusters)
        )
      )
    )
  
  stopifnot(
    !anyNA(bin_distributions$cluster)
  )
  
  list(
    theta = theta_run,
    contributions = contributions,
    bins = bin_distributions
  )
}

# Summarize each run separately.
long_summary_a <- summarize_long_run(
  clustering_long_a,
  "A"
)

long_summary_b <- summarize_long_run(
  clustering_long_b,
  "B"
)

# ------------------------------------------------------------
# 1. PRINT CONTRIBUTIONS
# ------------------------------------------------------------

long_run_contributions <- bind_rows(
  long_summary_a$contributions,
  long_summary_b$contributions
)

print(
  long_run_contributions,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 2. PLOT EACH RUN'S FITTED MOVEMENT DISTRIBUTIONS
# ------------------------------------------------------------

plot_run_distributions <- function(summary_object, title) {
  
  ggplot(
    summary_object$bins,
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
      title = title,
      x = "Movement-variable bin",
      y = "Fitted bin probability"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank()
    )
}

print(
  plot_run_distributions(
    long_summary_a,
    "Long run A"
  )
)

print(
  plot_run_distributions(
    long_summary_b,
    "Long run B"
  )
)


print(
  long_run_contributions,
  n = Inf,
  width = Inf
)


# ============================================================
# CHECK POSTERIOR STORAGE BEFORE WITHIN-RUN DIAGNOSTICS
# ============================================================

cat("\nbayesmove version:\n")
print(packageVersion("bayesmove"))

cat("\nRUN A: phi structure\n")
str(
  clustering_long_a$phi,
  max.level = 2,
  list.len = 2,
  vec.len = 3
)

cat("\nRUN B: phi structure\n")
str(
  clustering_long_b$phi,
  max.level = 2,
  list.len = 2,
  vec.len = 3
)

cat("\nRUN B: theta structure\n")
str(
  clustering_long_b$theta,
  max.level = 2,
  list.len = 2,
  vec.len = 3
)