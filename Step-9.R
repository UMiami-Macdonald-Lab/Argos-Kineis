# ============================================================
# BAYESMOVE STAGE 1: SEGMENT THE TRACKS
# ============================================================

library(bayesmove)
library(dplyr)
library(future)

# ------------------------------------------------------------
# 1. PREPARE THE SIX ORDERED TRACK SERIES
# ------------------------------------------------------------

tracks_for_segmentation <- tracks_discrete %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(
    obs = row_number(),
    time1 = row_number()
  ) %>%
  ungroup()

# Confirm that the retained movement rows remain daily.
segmentation_spacing <- tracks_for_segmentation %>%
  group_by(id) %>%
  summarise(
    n_rows = n(),
    
    all_daily = all(
      abs(
        as.numeric(diff(date), units = "secs") - 86400
      ) < 1
    ),
    
    .groups = "drop"
  )

print(
  segmentation_spacing,
  n = Inf,
  width = Inf
)

stopifnot(
  nrow(segmentation_spacing) == 6L,
  all(segmentation_spacing$all_daily),
  !anyNA(tracks_for_segmentation$SL),
  !anyNA(tracks_for_segmentation$TA)
)

# Keep a full-data list for later breakpoint assignment.
tracks_list <- bayesmove::df_to_list(
  dat = as.data.frame(tracks_for_segmentation),
  ind = "id"
)

# The segmentation model receives ONLY:
# id, discretized step length, discretized turning angle.
segmentation_input <- lapply(
  tracks_list,
  function(z) {
    as.data.frame(
      z[, c("id", "SL", "TA")]
    )
  }
)

# ------------------------------------------------------------
# 2. SET MODEL PARAMETERS
# ------------------------------------------------------------

seg_seed <- 20261103L

seg_iterations <- 10000L
seg_burn <- 5000L

# Order corresponds to SL, then TA.
nbins <- c(5L, 8L)

# Tutorial starting value for the Dirichlet prior.
seg_alpha <- 1

# Run sequentially for a straightforward reproducible setup.
future::plan(future::sequential)

set.seed(seg_seed)

# ------------------------------------------------------------
# 3. FIT SEGMENTATION MODELS
# ------------------------------------------------------------

segmentation_fit <- bayesmove::segment_behavior(
  data = segmentation_input,
  ngibbs = seg_iterations,
  nbins = nbins,
  alpha = seg_alpha
)

# Save the fitted models and inputs so they can be reused.
saveRDS(
  list(
    fit = segmentation_fit,
    tracks_list = tracks_list,
    seed = seg_seed,
    iterations = seg_iterations,
    burn_in = seg_burn,
    nbins = nbins,
    alpha = seg_alpha,
    step_limits = step_limits,
    angle_limits = angle_limits
  ),
  file = "Bayesmove_segmentation_checkpoint.rds"
)

# ------------------------------------------------------------
# 4. INSPECT SEGMENTATION DIAGNOSTICS
# ------------------------------------------------------------

# Log marginal likelihood through the MCMC iterations.
lml_plot <- bayesmove::traceplot(
  data = segmentation_fit,
  type = "LML"
)

print(lml_plot)

# Number of breakpoints through the MCMC iterations.
breakpoint_trace_plot <- bayesmove::traceplot(
  data = segmentation_fit,
  type = "nbrks"
)

print(breakpoint_trace_plot)


# Log marginal likelihood traces for all six tracks.
bayesmove::traceplot(
  data = segmentation_fit,
  type = "LML"
)