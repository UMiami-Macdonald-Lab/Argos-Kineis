# ============================================================
# PROPOSE SHARED BIN BOUNDARIES — NO MODEL FITTING YET
# ============================================================

library(dplyr)
library(tibble)

# Print the numerical movement summary from the previous step.
print(
  movement_summary,
  n = Inf,
  width = Inf
)

# Pool the usable movement rows across both systems.
# step is in METERS; angle is in RADIANS.
stopifnot(
  all(is.finite(movement_for_bayesmove$step)),
  all(is.finite(movement_for_bayesmove$angle))
)

# Five candidate step-length bins.
step_quantiles <- quantile(
  movement_for_bayesmove$step,
  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.90,
    1
  ),
  names = FALSE
)

# Check for repeated quantiles before discretizing.
step_limits <- unique(step_quantiles)

if (length(step_limits) < length(step_quantiles)) {
  warning(
    "Repeated step quantiles: fewer than five distinct ",
    "step-length bins are available."
  )
}

if (length(step_limits) < 3L) {
  stop("Too few distinct step boundaries. Review before proceeding.")
}

# Eight candidate turning-angle bins.
angle_limits <- seq(
  from = -pi,
  to = pi,
  length.out = 9
)

# Display step boundaries in kilometers.
step_bin_table <- tibble(
  bin = seq_len(length(step_limits) - 1L),
  lower_km = head(step_limits, -1) / 1000,
  upper_km = tail(step_limits, -1) / 1000
)

# Display turning-angle boundaries in degrees.
angle_bin_table <- tibble(
  bin = seq_len(length(angle_limits) - 1L),
  lower_degrees = head(angle_limits, -1) * 180 / pi,
  upper_degrees = tail(angle_limits, -1) * 180 / pi
)

print(
  step_bin_table,
  n = Inf,
  width = Inf
)

print(
  angle_bin_table,
  n = Inf,
  width = Inf
)

# Inspect the largest daily displacements.
largest_steps <- movement_for_bayesmove %>%
  select(
    track_id,
    system,
    date,
    step_km
  ) %>%
  arrange(desc(step_km)) %>%
  slice_head(n = 10)

print(
  largest_steps,
  n = Inf,
  width = Inf
)

# Proposed shared step-length boundaries.
print(
  step_bin_table,
  n = Inf,
  width = Inf
)

# Largest daily displacements and their timestamps.
print(
  largest_steps,
  n = Inf,
  width = Inf
)

# ============================================================
# DISCRETIZE MOVEMENT AND CHECK BIN OCCUPANCY
# ============================================================

library(dplyr)
library(tidyr)

# Use the boundaries calculated in the previous step.
stopifnot(
  length(step_limits) == 6L,
  length(angle_limits) == 9L,
  all(diff(step_limits) > 0),
  all(diff(angle_limits) > 0)
)

# SL = step-length bin
# TA = turning-angle bin
#
# include.lowest = TRUE retains the minimum value.
# With right = TRUE, internal boundaries belong to
# the interval ending at that boundary.

tracks_discrete <- movement_for_bayesmove %>%
  mutate(
    SL = as.integer(
      cut(
        step,
        breaks = step_limits,
        include.lowest = TRUE,
        right = TRUE
      )
    ),
    
    TA = as.integer(
      cut(
        angle,
        breaks = angle_limits,
        include.lowest = TRUE,
        right = TRUE
      )
    )
  ) %>%
  arrange(id, date)

# Confirm every movement row received valid bins.
stopifnot(
  nrow(tracks_discrete) == 330L,
  !anyNA(tracks_discrete$SL),
  !anyNA(tracks_discrete$TA),
  all(tracks_discrete$SL %in% 1:5),
  all(tracks_discrete$TA %in% 1:8)
)

# ------------------------------------------------------------
# 1. STEP-LENGTH BIN OCCUPANCY BY SYSTEM
# ------------------------------------------------------------

step_occupancy <- tracks_discrete %>%
  count(system, SL, name = "n") %>%
  complete(
    system,
    SL = 1:5,
    fill = list(n = 0L)
  ) %>%
  group_by(system) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

print(
  step_occupancy,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 2. TURNING-ANGLE BIN OCCUPANCY BY SYSTEM
# ------------------------------------------------------------

angle_occupancy <- tracks_discrete %>%
  count(system, TA, name = "n") %>%
  complete(
    system,
    TA = 1:8,
    fill = list(n = 0L)
  ) %>%
  group_by(system) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

print(
  angle_occupancy,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 3. CHECK THE POOLED STEP-LENGTH PROPORTIONS
# ------------------------------------------------------------

pooled_step_occupancy <- tracks_discrete %>%
  count(SL, name = "n") %>%
  mutate(
    percent = 100 * n / sum(n)
  )

print(
  pooled_step_occupancy,
  n = Inf,
  width = Inf
)