# ================================================================
# UPDATED SUPPLEMENTAL 6: FIGURES AND TABLES FROM SAVED RESULTS
#
# Reference = long run A
# Sensitivity = long run B
#
# No models are refitted.
# No manuscript or Google Doc is edited.
#
# Required files in data_dir:
#   Bayesmove_segmentation_checkpoint.rds
#   Bayesmove_clustering_checkpoint.rds
#   Bayesmove_long_clustering_checkpoint.rds
#   Argos-SSM_daily_aligned.csv
#   Doppler-SSM_daily_aligned.csv
#
# Optional:
#   SSM_fit_status.csv
#
# Outputs:
#   PNG/PDF figures
#   Full-precision CSV tables
#   Figure/table captions
#   Interpretation notes
#
# Date labels:
#   Tracks spanning <=35 days: every 5 days
#   Longer tracks: every 14 days
# ================================================================


# ------------------------------------------------------------
# 0. SETTINGS AND PACKAGES
# ------------------------------------------------------------

# Change this to your analysis folder.
data_dir <- "."

out_dir <- file.path(
  data_dir,
  "Supplemental_6_Updated_Outputs"
)

figure_dpi <- 600

working_ids <- c(
  "34885",
  "34890",
  "45270"
)

# Install once if needed:
#
# install.packages(c(
#   "bayesmove",
#   "dplyr",
#   "tidyr",
#   "tibble",
#   "readr",
#   "ggplot2",
#   "patchwork"
# ))

packages <- c(
  "bayesmove",
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "ggplot2",
  "patchwork"
)

missing_packages <- packages[
  !vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages)) {
  stop(
    "Install these packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------

write_table <- function(x, filename) {
  readr::write_csv(
    x,
    file.path(out_dir, filename)
  )
}

read_checkpoint <- function(filename) {
  
  path <- file.path(
    data_dir,
    filename
  )
  
  if (!file.exists(path)) {
    stop("Missing checkpoint: ", path)
  }
  
  readRDS(path)
}

save_figure <- function(
    p,
    filename,
    width,
    height
) {
  
  ggplot2::ggsave(
    file.path(
      out_dir,
      paste0(filename, ".png")
    ),
    plot = p,
    width = width,
    height = height,
    dpi = figure_dpi,
    bg = "white"
  )
  
  ggplot2::ggsave(
    file.path(
      out_dir,
      paste0(filename, ".pdf")
    ),
    plot = p,
    width = width,
    height = height,
    bg = "white"
  )
}

percent_label <- function(x) {
  paste0(
    round(100 * x),
    "%"
  )
}

# More date labels for short tracks such as 45270.
adaptive_date_breaks <- function(limits) {
  
  first_day <- as.POSIXct(
    as.Date(
      limits[1],
      tz = "UTC"
    ),
    tz = "UTC"
  )
  
  last_day <- as.POSIXct(
    as.Date(
      limits[2],
      tz = "UTC"
    ),
    tz = "UTC"
  )
  
  span_days <- as.numeric(
    difftime(
      last_day,
      first_day,
      units = "days"
    )
  )
  
  seq(
    from = first_day,
    to = last_day,
    by = if (span_days <= 35) {
      "5 days"
    } else {
      "14 days"
    }
  )
}


# ------------------------------------------------------------
# 1. LOAD AND VALIDATE
# ------------------------------------------------------------

seg_checkpoint <- read_checkpoint(
  "Bayesmove_segmentation_checkpoint.rds"
)

base_checkpoint <- read_checkpoint(
  "Bayesmove_clustering_checkpoint.rds"
)

long_checkpoint <- read_checkpoint(
  "Bayesmove_long_clustering_checkpoint.rds"
)

tracks <- as_tibble(
  base_checkpoint$tracks_segmented
) %>%
  mutate(
    track_id = as.character(track_id),
    id = as.character(id)
  ) %>%
  arrange(id, date)

counts <- as.data.frame(
  base_checkpoint$segment_bin_counts
)

counts$id <- as.character(counts$id)

step_limits <- base_checkpoint$step_limits
angle_limits <- base_checkpoint$angle_limits

n_iterations <- long_checkpoint$iterations
n_burn <- long_checkpoint$burn_in
n_clusters <- long_checkpoint$max_clusters

cluster_cols <- paste0(
  "cluster_",
  seq_len(n_clusters)
)

stopifnot(
  n_clusters == 7L,
  n_iterations > n_burn,
  
  all(
    c(
      "id",
      "track_id",
      "system",
      "date",
      "tseg",
      "step",
      "angle"
    ) %in% names(tracks)
  ),
  
  all(working_ids %in% tracks$track_id),
  
  all(is.finite(tracks$step)),
  all(tracks$step > 0),
  
  all(is.finite(tracks$angle)),
  !anyNA(tracks$date),
  
  anyDuplicated(
    counts[c("id", "tseg")]
  ) == 0L,
  
  anyDuplicated(
    tracks[c("track_id", "system", "date")]
  ) == 0L
)

if (!inherits(tracks$date, "POSIXct")) {
  stop("Checkpoint date must be POSIXct.")
}

sl_columns <- grep(
  "^y1\\.",
  names(counts),
  value = TRUE
)

ta_columns <- grep(
  "^y2\\.",
  names(counts),
  value = TRUE
)

weights <- rowSums(
  as.matrix(
    counts[, sl_columns, drop = FALSE]
  )
)

stopifnot(
  length(sl_columns) == length(step_limits) - 1L,
  length(ta_columns) == length(angle_limits) - 1L,
  
  all(weights > 0),
  
  sum(weights) == nrow(tracks),
  
  all(
    weights ==
      rowSums(
        as.matrix(
          counts[, ta_columns, drop = FALSE]
        )
      )
  )
)

# Verify segment keys and counts.
actual_counts <- tracks %>%
  count(
    id,
    tseg,
    name = "actual_n"
  )

count_check <- as_tibble(
  counts[, c("id", "tseg")]
) %>%
  mutate(
    model_n = weights
  ) %>%
  full_join(
    actual_counts,
    by = c("id", "tseg")
  )

stopifnot(
  !anyNA(count_check),
  all(
    count_check$model_n ==
      count_check$actual_n
  )
)


# ------------------------------------------------------------
# EXTRACT EACH LONG RUN
# ------------------------------------------------------------

extract_run <- function(fit, run_name) {
  
  th <- as.matrix(
    bayesmove::extract_prop(
      fit,
      ngibbs = n_iterations,
      nburn = n_burn,
      nmaxclust = n_clusters
    )
  )
  
  stopifnot(
    nrow(th) == nrow(counts),
    ncol(th) == n_clusters,
    all(is.finite(th)),
    all(th >= 0),
    all(abs(rowSums(th) - 1) < 1e-5)
  )
  
  colnames(th) <- cluster_cols
  
  contributions <- tibble(
    run = run_name,
    cluster = cluster_cols,
    
    mean_segment_prop = colMeans(th),
    
    mean_daily_row_prop =
      colSums(
        sweep(
          th,
          1,
          weights,
          "*"
        )
      ) / sum(weights)
  )
  
  segment_rank <- contributions %>%
    arrange(
      desc(mean_segment_prop)
    ) %>%
    mutate(
      segment_rank = row_number(),
      
      cumulative_segment_prop =
        cumsum(mean_segment_prop),
      
      retained_segment_90 =
        lag(
          cumulative_segment_prop,
          default = 0
        ) <= 0.90
    ) %>%
    select(
      run,
      cluster,
      segment_rank,
      cumulative_segment_prop,
      retained_segment_90
    )
  
  daily_rank <- contributions %>%
    arrange(
      desc(mean_daily_row_prop)
    ) %>%
    mutate(
      daily_rank = row_number(),
      
      cumulative_daily_row_prop =
        cumsum(mean_daily_row_prop),
      
      retained_daily_90 =
        lag(
          cumulative_daily_row_prop,
          default = 0
        ) <= 0.90
    ) %>%
    select(
      run,
      cluster,
      daily_rank,
      cumulative_daily_row_prop,
      retained_daily_90
    )
  
  contributions <- contributions %>%
    left_join(
      segment_rank,
      by = c("run", "cluster")
    ) %>%
    left_join(
      daily_rank,
      by = c("run", "cluster")
    ) %>%
    arrange(segment_rank)
  
  bins <- as_tibble(
    bayesmove::get_behav_hist(
      fit,
      nburn = n_burn,
      ngibbs = n_iterations,
      nmaxclust = n_clusters,
      
      var.names = c(
        "Step length",
        "Turning angle"
      )
    )
  ) %>%
    mutate(
      run = run_name,
      
      cluster_number = as.integer(
        as.character(behav)
      )
    )
  
  stopifnot(
    !anyNA(bins$cluster_number),
    
    all(
      bins$cluster_number %in%
        seq_len(n_clusters)
    )
  )
  
  list(
    theta = th,
    contributions = contributions,
    bins = bins
  )
}

run_a <- extract_run(
  long_checkpoint$run_a,
  "A"
)

run_b <- extract_run(
  long_checkpoint$run_b,
  "B"
)

contributions <- bind_rows(
  run_a$contributions,
  run_b$contributions
)

segment_props <- bind_cols(
  as_tibble(
    counts[, c("id", "tseg")]
  ),
  as_tibble(run_a$theta)
)

reference_rows <- tracks %>%
  select(
    -any_of(cluster_cols)
  ) %>%
  left_join(
    segment_props,
    by = c("id", "tseg")
  )

p_matrix <- as.matrix(
  reference_rows[, cluster_cols]
)

stopifnot(
  nrow(reference_rows) == nrow(tracks),
  all(is.finite(p_matrix)),
  all(abs(rowSums(p_matrix) - 1) < 1e-5)
)

reference_rows$other <- rowSums(
  p_matrix[, 4:7, drop = FALSE]
)


# ------------------------------------------------------------
# 2. FIGURE SETTINGS
# ------------------------------------------------------------

state_levels <- c(
  "cluster_1",
  "cluster_2",
  "cluster_3",
  "other"
)

state_labels <- c(
  "Shorter steps, broad turns",
  "Longer steps, mixed turns",
  "Shorter steps, small turns",
  "Remaining clusters"
)

state_colors <- setNames(
  c(
    "#3B8C95",
    "#E3AA42",
    "#8664A5",
    "#BDBDBD"
  ),
  state_labels
)

system_colors <- c(
  "Argos" = "#5FAF89",
  "Kin\u00e9is" = "#62B7C8"
)

label_system <- function(x) {
  factor(
    x,
    levels = c(
      "Argos",
      "Kineis"
    ),
    labels = c(
      "Argos",
      "Kin\u00e9is"
    )
  )
}

paper_theme <- theme_bw(
  base_size = 10
) +
  theme(
    panel.grid.minor = element_blank(),
    
    strip.background = element_rect(
      fill = "grey95",
      color = "grey60"
    ),
    
    legend.position = "bottom",
    legend.title = element_blank()
  )


# ------------------------------------------------------------
# 3. TABLES
# ------------------------------------------------------------

# TABLE 6.1: COVERAGE AND SEGMENTATION

segment_lengths <- tracks %>%
  group_by(
    track_id,
    system,
    id,
    tseg
  ) %>%
  summarise(
    n_movement_rows = n(),
    first_time = min(date),
    last_time = max(date),
    .groups = "drop"
  )

seg_summary <- segment_lengths %>%
  group_by(
    track_id,
    system
  ) %>%
  summarise(
    n_segments = n(),
    n_breakpoints = n() - 1L,
    
    shortest_segment_rows =
      min(n_movement_rows),
    
    median_segment_rows =
      median(n_movement_rows),
    
    longest_segment_rows =
      max(n_movement_rows),
    
    n_single_row_segments =
      sum(n_movement_rows == 1L),
    
    n_complete_movement_rows =
      sum(n_movement_rows),
    
    .groups = "drop"
  )

read_aligned <- function(
    filename,
    system_name
) {
  
  path <- file.path(
    data_dir,
    filename
  )
  
  if (!file.exists(path)) {
    stop(
      "Missing aligned predictions: ",
      path
    )
  }
  
  readr::read_csv(
    path,
    
    col_types = readr::cols(
      track_id = readr::col_character(),
      date = readr::col_datetime()
    ),
    
    show_col_types = FALSE
  ) %>%
    transmute(
      track_id,
      system = system_name,
      date,
      lon,
      lat
    ) %>%
    filter(
      track_id %in% working_ids
    )
}

aligned <- bind_rows(
  read_aligned(
    "Argos-SSM_daily_aligned.csv",
    "Argos"
  ),
  
  read_aligned(
    "Doppler-SSM_daily_aligned.csv",
    "Kineis"
  )
)

coverage <- aligned %>%
  group_by(
    track_id,
    system
  ) %>%
  summarise(
    n_daily_positions = n(),
    first_prediction = min(date),
    last_prediction = max(date),
    .groups = "drop"
  )

species_key <- tibble(
  track_id = working_ids,
  
  species = c(
    "Tiger shark",
    "Tiger shark",
    "Bull shark"
  )
)

table_61 <- coverage %>%
  left_join(
    species_key,
    by = "track_id"
  ) %>%
  left_join(
    seg_summary,
    by = c("track_id", "system")
  ) %>%
  select(
    track_id,
    species,
    system,
    everything()
  ) %>%
  arrange(
    track_id,
    system
  )

write_table(
  table_61,
  "Table_6_1_Coverage_and_segmentation.csv"
)

write_table(
  segment_lengths,
  "Supporting_segment_lengths.csv"
)


# TABLE 6.2: BIN BOUNDARIES

boundaries <- bind_rows(
  tibble(
    var = "Step length",
    
    bin = seq_len(
      length(step_limits) - 1L
    ),
    
    lower = head(step_limits, -1) / 1000,
    upper = tail(step_limits, -1) / 1000,
    units = "km"
  ),
  
  tibble(
    var = "Turning angle",
    
    bin = seq_len(
      length(angle_limits) - 1L
    ),
    
    lower = head(angle_limits, -1) * 180 / pi,
    upper = tail(angle_limits, -1) * 180 / pi,
    units = "degrees"
  )
)

write_table(
  boundaries,
  "Table_6_2_Shared_bin_boundaries.csv"
)


# TABLE 6.3: CLUSTER CONTRIBUTIONS

write_table(
  contributions,
  "Table_6_3_Cluster_contributions_both_runs.csv"
)


# TABLE 6.4: INDIVIDUAL COMPOSITIONS

budgets <- reference_rows %>%
  group_by(
    track_id,
    system
  ) %>%
  summarise(
    n_movement_rows = n(),
    
    across(
      all_of(cluster_cols),
      ~ mean(.x)
    ),
    
    .groups = "drop"
  )

write_table(
  budgets,
  "Supporting_individual_composition_proportions.csv"
)

write_table(
  budgets %>%
    mutate(
      across(
        all_of(cluster_cols),
        ~ 100 * .x
      )
    ),
  
  "Table_6_4_Individual_composition_all_clusters_percent.csv"
)

major_budgets <- budgets %>%
  transmute(
    track_id,
    system,
    n_movement_rows,
    
    cluster_1_percent = 100 * cluster_1,
    cluster_2_percent = 100 * cluster_2,
    cluster_3_percent = 100 * cluster_3,
    
    remaining_percent =
      100 * (
        1 - cluster_1 - cluster_2 - cluster_3
      )
  )

write_table(
  major_budgets,
  "Table_6_4_Individual_composition_display_percent.csv"
)


# TABLE 6.5: PAIRED DIFFERENCES

paired_differences <- budgets %>%
  select(
    track_id,
    system,
    all_of(cluster_cols)
  ) %>%
  pivot_longer(
    all_of(cluster_cols),
    names_to = "cluster",
    values_to = "proportion"
  ) %>%
  pivot_wider(
    names_from = system,
    values_from = proportion
  ) %>%
  mutate(
    Argos_percent = 100 * Argos,
    Kineis_percent = 100 * Kineis,
    
    Kineis_minus_Argos_pp =
      100 * (Kineis - Argos)
  )

stopifnot(
  !anyNA(paired_differences$Argos),
  !anyNA(paired_differences$Kineis)
)

write_table(
  paired_differences,
  "Table_6_5_Paired_composition_differences.csv"
)

paired_summary <- paired_differences %>%
  group_by(cluster) %>%
  summarise(
    n_sharks = n(),
    
    mean_difference_pp =
      mean(Kineis_minus_Argos_pp),
    
    median_difference_pp =
      median(Kineis_minus_Argos_pp),
    
    min_difference_pp =
      min(Kineis_minus_Argos_pp),
    
    max_difference_pp =
      max(Kineis_minus_Argos_pp),
    
    .groups = "drop"
  )

write_table(
  paired_summary,
  "Supporting_equal_shark_paired_difference_summary.csv"
)


# TABLE 6.6: DOMINANT-CLUSTER AGREEMENT
#
# This uses all seven reference-run clusters.
# Ties are left unclassified.
# Mixed compositions remain the main outcome.

row_max <- apply(
  p_matrix,
  1,
  max
)

ties <- rowSums(
  abs(
    sweep(
      p_matrix,
      1,
      row_max,
      "-"
    )
  ) <= 1e-10
) > 1L

reference_rows$dominant_cluster <- cluster_cols[
  max.col(
    p_matrix,
    ties.method = "first"
  )
]

reference_rows$dominant_cluster[ties] <- NA_character_

get_system_rows <- function(system_name) {
  
  reference_rows %>%
    filter(
      system == system_name
    ) %>%
    select(
      track_id,
      date,
      dominant_cluster,
      all_of(cluster_cols)
    )
}

argos_rows <- get_system_rows("Argos")
kineis_rows <- get_system_rows("Kineis")

matched <- inner_join(
  argos_rows,
  kineis_rows,
  by = c("track_id", "date"),
  suffix = c("_Argos", "_Kineis")
)

stopifnot(
  nrow(matched) == nrow(argos_rows),
  nrow(matched) == nrow(kineis_rows)
)

a_mat <- as.matrix(
  matched[
    ,
    paste0(cluster_cols, "_Argos")
  ]
)

k_mat <- as.matrix(
  matched[
    ,
    paste0(cluster_cols, "_Kineis")
  ]
)

# Sum of minimum proportions across all clusters:
# 1 = identical compositions
# 0 = non-overlapping compositions
# This is not a probability of behavioral agreement.

matched$composition_similarity <- rowSums(
  pmin(
    a_mat,
    k_mat
  )
)

matched <- matched %>%
  mutate(
    comparable =
      !is.na(dominant_cluster_Argos) &
      !is.na(dominant_cluster_Kineis),
    
    agreement = if_else(
      comparable,
      dominant_cluster_Argos == dominant_cluster_Kineis,
      NA
    )
  )

agreement_summary <- function(x) {
  
  x %>%
    summarise(
      n_matched_movement_days = n(),
      
      n_comparable = sum(comparable),
      
      n_agreements = sum(
        agreement,
        na.rm = TRUE
      ),
      
      agreement_percent = if (
        sum(comparable) > 0
      ) {
        100 *
          sum(agreement, na.rm = TRUE) /
          sum(comparable)
      } else {
        NA_real_
      },
      
      mean_composition_similarity =
        mean(composition_similarity),
      
      .groups = "drop"
    )
}

agreement_by_shark <- agreement_summary(
  group_by(
    matched,
    track_id
  )
)

pooled_agreement <- agreement_summary(
  matched
)

write_table(
  agreement_by_shark,
  "Table_6_6_Matched_day_agreement_by_shark.csv"
)

write_table(
  pooled_agreement,
  "Supporting_pooled_day_agreement.csv"
)

equal_shark_agreement <- agreement_by_shark %>%
  summarise(
    n_sharks = n(),
    
    mean_shark_agreement_percent = if (
      all(is.na(agreement_percent))
    ) {
      NA_real_
    } else {
      mean(
        agreement_percent,
        na.rm = TRUE
      )
    },
    
    mean_shark_composition_similarity =
      mean(mean_composition_similarity)
  )

write_table(
  equal_shark_agreement,
  "Supporting_equal_shark_agreement.csv"
)

write_table(
  matched,
  "Supporting_matched_day_comparisons.csv"
)

write_table(
  reference_rows,
  "Supporting_reference_movement_and_compositions.csv"
)


# ADDITIONAL MOVEMENT-METRIC SUMMARIES

movement_stats <- tracks %>%
  group_by(
    track_id,
    system
  ) %>%
  summarise(
    n_movement_rows = n(),
    
    min_step_km = min(step) / 1000,
    
    q25_step_km =
      unname(quantile(step, 0.25)) / 1000,
    
    median_step_km =
      median(step) / 1000,
    
    mean_step_km =
      mean(step) / 1000,
    
    q75_step_km =
      unname(quantile(step, 0.75)) / 1000,
    
    q90_step_km =
      unname(quantile(step, 0.90)) / 1000,
    
    max_step_km =
      max(step) / 1000,
    
    median_abs_angle_degrees =
      median(abs(angle)) * 180 / pi,
    
    .groups = "drop"
  )

write_table(
  movement_stats,
  "Supporting_movement_metric_summary.csv"
)

largest_steps <- tracks %>%
  select(
    track_id,
    system,
    date,
    step,
    angle
  ) %>%
  arrange(desc(step)) %>%
  slice_head(n = 10)

write_table(
  largest_steps,
  "Supporting_largest_steps_meters.csv"
)

dominant_metrics <- reference_rows %>%
  filter(
    !is.na(dominant_cluster)
  ) %>%
  group_by(
    dominant_cluster
  ) %>%
  summarise(
    n_movement_rows = n(),
    
    median_step_m = median(step),
    mean_step_m = mean(step),
    
    median_log_step = median(log(step)),
    mean_log_step = mean(log(step)),
    
    median_abs_angle_rad =
      median(abs(angle)),
    
    mean_abs_angle_rad =
      mean(abs(angle)),
    
    .groups = "drop"
  )

write_table(
  dominant_metrics,
  "Supporting_metrics_grouped_by_dominant_cluster.csv"
)

fitted_bins_both_runs <- bind_rows(
  run_a$bins,
  run_b$bins
) %>%
  left_join(
    boundaries,
    by = c("var", "bin")
  )

write_table(
  fitted_bins_both_runs,
  "Supporting_fitted_bin_probabilities_both_runs.csv"
)


# ------------------------------------------------------------
# 4. FIGURE 6.1: DISCRETIZATION
# ------------------------------------------------------------

plot_tracks <- tracks %>%
  mutate(
    system_label = label_system(system),
    step_km = step / 1000,
    angle_degrees = angle * 180 / pi
  )

p_step <- ggplot(
  plot_tracks,
  
  aes(
    step_km,
    y = after_stat(density),
    fill = system_label
  )
) +
  geom_histogram(
    breaks = seq(
      0,
      max(plot_tracks$step_km) * 1.000001,
      length.out = 31
    ),
    
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = step_limits / 1000,
    linetype = 2,
    linewidth = 0.35
  ) +
  facet_wrap(
    ~ system_label,
    ncol = 1
  ) +
  scale_fill_manual(
    values = system_colors,
    guide = "none"
  ) +
  labs(
    x = "Daily step length (km)",
    y = "Density"
  ) +
  paper_theme

p_angle <- ggplot(
  plot_tracks,
  
  aes(
    angle_degrees,
    y = after_stat(density),
    fill = system_label
  )
) +
  geom_histogram(
    breaks = seq(
      -180,
      180,
      length.out = 17
    ),
    
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = angle_limits * 180 / pi,
    linetype = 2,
    linewidth = 0.35
  ) +
  facet_wrap(
    ~ system_label,
    ncol = 1
  ) +
  scale_fill_manual(
    values = system_colors,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(
      -180,
      -90,
      0,
      90,
      180
    )
  ) +
  labs(
    x = "Signed turning angle (degrees)",
    y = "Density"
  ) +
  paper_theme

fig_61 <- (
  p_step | p_angle
) +
  plot_annotation(
    tag_levels = "A"
  )

save_figure(
  fig_61,
  "Figure_6_1_Movement_discretization",
  10,
  6
)


# ------------------------------------------------------------
# 5. FIGURE 6.2: THETA SUMMARY
# ------------------------------------------------------------

theta_long <- as_tibble(
  run_a$theta
) %>%
  pivot_longer(
    everything(),
    names_to = "cluster",
    values_to = "prop"
  ) %>%
  mutate(
    cluster = factor(
      cluster,
      levels = cluster_cols,
      labels = seq_len(n_clusters)
    )
  )

p_theta <- ggplot(
  theta_long,
  aes(
    cluster,
    prop
  )
) +
  geom_boxplot(
    fill = "#9BC6CB",
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 1
    ),
    
    alpha = 0.3,
    size = 1
  ) +
  labs(
    x = "Run A cluster",
    y = "Segment composition"
  ) +
  paper_theme

mean_data <- run_a$contributions %>%
  select(
    cluster,
    mean_segment_prop,
    mean_daily_row_prop
  ) %>%
  pivot_longer(
    -cluster,
    names_to = "weighting",
    values_to = "prop"
  ) %>%
  mutate(
    cluster = factor(
      cluster,
      levels = cluster_cols,
      labels = seq_len(n_clusters)
    ),
    
    weighting = factor(
      weighting,
      
      levels = c(
        "mean_segment_prop",
        "mean_daily_row_prop"
      ),
      
      labels = c(
        "Equal segment",
        "Daily row"
      )
    )
  )

p_means <- ggplot(
  mean_data,
  aes(
    cluster,
    prop,
    fill = weighting
  )
) +
  geom_col(
    position = "dodge"
  ) +
  scale_fill_manual(
    values = c(
      "Equal segment" = "#3B8C95",
      "Daily row" = "#E3AA42"
    )
  ) +
  labs(
    x = "Run A cluster",
    y = "Mean contribution"
  ) +
  paper_theme

cum_data <- bind_rows(
  run_a$contributions %>%
    transmute(
      rank = segment_rank,
      cumulative = cumulative_segment_prop,
      weighting = "Equal segment"
    ),
  
  run_a$contributions %>%
    transmute(
      rank = daily_rank,
      cumulative = cumulative_daily_row_prop,
      weighting = "Daily row"
    )
)

p_cumulative <- ggplot(
  cum_data,
  aes(
    rank,
    cumulative,
    color = weighting
  )
) +
  geom_hline(
    yintercept = 0.9,
    linetype = 2,
    color = "grey40"
  ) +
  geom_line() +
  geom_point() +
  scale_color_manual(
    values = c(
      "Equal segment" = "#3B8C95",
      "Daily row" = "#E3AA42"
    )
  ) +
  scale_x_continuous(
    breaks = seq_len(n_clusters)
  ) +
  scale_y_continuous(
    labels = percent_label
  ) +
  labs(
    x = "Contribution rank within weighting",
    y = "Cumulative contribution"
  ) +
  paper_theme

fig_62 <- (
  p_theta | p_means | p_cumulative
) +
  plot_annotation(
    tag_levels = "A"
  )

save_figure(
  fig_62,
  "Figure_6_2_Reference_cluster_contributions",
  13,
  4.5
)


# ------------------------------------------------------------
# 6. FIGURE 6.3: FITTED DISTRIBUTIONS
# ------------------------------------------------------------

bin_plot <- function(
    bins,
    title,
    clusters_to_show
) {
  
  bins <- bins %>%
    filter(
      cluster_number %in%
        clusters_to_show
    ) %>%
    mutate(
      cluster_label = factor(
        cluster_number,
        levels = clusters_to_show,
        
        labels = paste(
          "Cluster",
          clusters_to_show
        )
      )
    )
  
  ggplot(
    bins,
    aes(
      bin,
      prop
    )
  ) +
    geom_col(
      fill = "#3B8C95",
      width = 0.8
    ) +
    facet_grid(
      cluster_label ~ var,
      scales = "free_x"
    ) +
    scale_x_continuous(
      breaks = 1:8
    ) +
    scale_y_continuous(
      breaks = c(
        0,
        0.5,
        1
      )
    ) +
    labs(
      title = title,
      x = "Movement-variable bin",
      y = "Fitted bin probability"
    ) +
    paper_theme
}

fig_63 <- bin_plot(
  run_a$bins,
  "Reference run A",
  1:7
)

save_figure(
  fig_63,
  "Figure_6_3_Reference_fitted_distributions",
  9,
  10
)


# ------------------------------------------------------------
# 7. FIGURE 6.4: COMPOSITION
# ------------------------------------------------------------

composition_long <- reference_rows %>%
  select(
    track_id,
    system,
    date,
    all_of(state_levels)
  ) %>%
  pivot_longer(
    all_of(state_levels),
    names_to = "cluster",
    values_to = "proportion"
  ) %>%
  mutate(
    state = factor(
      cluster,
      levels = state_levels,
      labels = state_labels
    ),
    
    system_label = label_system(system)
  )

# Improved dates:
# 45270 receives five-day spacing.
# Longer tracks retain fourteen-day spacing.

p_time <- ggplot(
  composition_long,
  
  aes(
    date,
    proportion,
    fill = state
  )
) +
  geom_col(
    width = 86400 * 0.95
  ) +
  facet_wrap(
    vars(
      track_id,
      system_label
    ),
    ncol = 2,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = state_colors,
    drop = FALSE
  ) +
  scale_x_datetime(
    breaks = adaptive_date_breaks,
    date_labels = "%b %d",
    timezone = "UTC"
  ) +
  scale_y_continuous(
    breaks = c(
      0,
      0.5,
      1
    ),
    
    labels = percent_label
  ) +
  labs(
    x = "Date (UTC)",
    y = "Estimated segment composition"
  ) +
  paper_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    )
  )

plot_budgets <- composition_long %>%
  group_by(
    track_id,
    system_label,
    state
  ) %>%
  summarise(
    proportion = mean(proportion),
    .groups = "drop"
  )

p_budget <- ggplot(
  plot_budgets,
  
  aes(
    system_label,
    proportion,
    fill = state
  )
) +
  geom_col(
    width = 0.65
  ) +
  facet_wrap(
    ~ track_id,
    nrow = 1
  ) +
  scale_fill_manual(
    values = state_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = c(
      0,
      0.5,
      1
    ),
    
    labels = percent_label
  ) +
  labs(
    x = NULL,
    y = "Mean daily-row composition"
  ) +
  paper_theme

fig_64 <- (
  p_time / p_budget
) +
  plot_layout(
    heights = c(
      3,
      1.1
    ),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  )

fig_64 <- fig_64 &
  theme(
    legend.position = "bottom"
  ) &
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )

save_figure(
  fig_64,
  "Figure_6_4_Reference_compositions",
  11,
  11
)


# ------------------------------------------------------------
# 8. FIGURE 6.5: RUN SENSITIVITY
# ------------------------------------------------------------

# Cluster numbers are not assumed to match across runs.

trace_data <- bind_rows(
  tibble(
    run = "A",
    
    iteration = seq_along(
      long_checkpoint$run_a$loglikel
    ),
    
    loglik = as.numeric(
      long_checkpoint$run_a$loglikel
    )
  ),
  
  tibble(
    run = "B",
    
    iteration = seq_along(
      long_checkpoint$run_b$loglikel
    ),
    
    loglik = as.numeric(
      long_checkpoint$run_b$loglikel
    )
  )
)

p_trace <- ggplot(
  trace_data,
  aes(
    iteration,
    loglik
  )
) +
  geom_line(
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = n_burn,
    color = "#C83C3C",
    linetype = 2
  ) +
  facet_wrap(
    ~ run,
    nrow = 1
  ) +
  labs(
    x = "Iteration",
    y = "Log likelihood"
  ) +
  paper_theme

p_bins_a <- bin_plot(
  run_a$bins,
  "Run A: clusters 1–4",
  1:4
)

p_bins_b <- bin_plot(
  run_b$bins,
  "Run B: clusters 1–4",
  1:4
)

fig_65 <- (
  p_trace /
    (p_bins_a | p_bins_b)
) +
  plot_layout(
    heights = c(
      1,
      3
    )
  )

save_figure(
  fig_65,
  "Figure_6_5_Clustering_sensitivity",
  12,
  11
)


# ------------------------------------------------------------
# EXTRA DIAGNOSTICS
# ------------------------------------------------------------

save_figure(
  bin_plot(
    run_b$bins,
    "Sensitivity run B: all clusters",
    1:7
  ),
  
  "Supporting_Run_B_all_cluster_distributions",
  9,
  10
)

pdf(
  file.path(
    out_dir,
    "Supporting_segmentation_diagnostics.pdf"
  ),
  width = 11,
  height = 8
)

tryCatch(
  {
    bayesmove::traceplot(
      data = seg_checkpoint$fit,
      type = "LML"
    )
    
    bayesmove::traceplot(
      data = seg_checkpoint$fit,
      type = "nbrks"
    )
  },
  
  finally = dev.off()
)



# ------------------------------------------------------------
# OPTIONAL INITIAL SSM TRACK COUNTS
# ------------------------------------------------------------

status_path <- file.path(
  data_dir,
  "SSM_fit_status.csv"
)

if (file.exists(status_path)) {
  
  initial_status <- readr::read_csv(
    status_path,
    
    col_types = readr::cols(
      track_id = readr::col_character()
    ),
    
    show_col_types = FALSE
  )
  
  initial_counts <- initial_status %>%
    filter(
      status == "fit_returned_review_diagnostics"
    ) %>%
    distinct(
      system,
      track_id
    ) %>%
    count(
      system,
      name = "n_tracks_returning_predictions"
    )
  
  write_table(
    initial_counts,
    "Supporting_initial_SSM_returned_track_counts.csv"
  )
}


# ------------------------------------------------------------
# PRINT MAIN RESULTS
# ------------------------------------------------------------

cat("\nUPDATED INDIVIDUAL COMPOSITION (%)\n")

print(
  major_budgets,
  n = Inf,
  width = Inf
)

cat("\nMATCHED-DAY AGREEMENT\n")

print(
  agreement_by_shark,
  n = Inf,
  width = Inf
)

cat("\nPOOLED AGREEMENT\n")

print(
  pooled_agreement,
  width = Inf
)

message(
  "\nFinished. Outputs: ",
  normalizePath(out_dir)
)


# Display all five figures.
print(fig_61)
print(fig_62)
print(fig_63)
print(fig_64)
print(fig_65)