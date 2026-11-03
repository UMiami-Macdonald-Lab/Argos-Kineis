# ============================================================
# EXTRACT BREAKPOINTS AND INSPECT THE RESULTING SEGMENTS
# ============================================================

library(bayesmove)
library(dplyr)

# ------------------------------------------------------------
# 1. SELECT MAP ITERATIONS AFTER BURN-IN
# ------------------------------------------------------------

map_iterations <- bayesmove::get_MAP(
  dat = segmentation_fit$LML,
  nburn = seg_burn
)

print(map_iterations)

# ------------------------------------------------------------
# 2. EXTRACT THE CORRESPONDING BREAKPOINTS
# ------------------------------------------------------------

breakpoints <- bayesmove::get_breakpts(
  dat = segmentation_fit$brkpts,
  MAP.est = map_iterations
)

print(breakpoints)

# ------------------------------------------------------------
# 3. ASSIGN SEGMENT NUMBERS TO MOVEMENT ROWS
# ------------------------------------------------------------

tracks_segmented <- bayesmove::assign_tseg(
  dat = tracks_list,
  brkpts = breakpoints
)

# Confirm all movement rows were retained and assigned.
stopifnot(
  nrow(tracks_segmented) ==
    nrow(tracks_for_segmentation),
  !anyNA(tracks_segmented$tseg)
)

# ------------------------------------------------------------
# 4. SUMMARIZE SEGMENT LENGTHS
# ------------------------------------------------------------

segment_lengths <- tracks_segmented %>%
  group_by(
    track_id,
    system,
    id,
    tseg
  ) %>%
  summarise(
    n_movement_rows = n(),
    first_movement_time = min(date),
    last_movement_time = max(date),
    .groups = "drop"
  )

segmentation_summary <- segment_lengths %>%
  group_by(
    track_id,
    system,
    id
  ) %>%
  summarise(
    n_segments = n(),
    n_breakpoints = n() - 1L,
    
    shortest_segment_rows = min(n_movement_rows),
    
    median_segment_rows = median(n_movement_rows),
    
    longest_segment_rows = max(n_movement_rows),
    
    total_movement_rows = sum(n_movement_rows),
    
    .groups = "drop"
  )

print(
  segmentation_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 5. PLOT BREAKPOINTS OVER THE MOVEMENT VARIABLES
# ------------------------------------------------------------

bayesmove::plot_breakpoints(
  data = tracks_list,
  as_date = TRUE,
  var_names = c(
    "step_km",
    "angle"
  ),
  var_labels = c(
    "Daily step length (km)",
    "Turning angle (rad)"
  ),
  brkpts = breakpoints
)