# ============================================================
# UPDATED BAYESMOVE FIGURES AND RESULTS
# Reference analysis: long run A
# ============================================================

library(bayesmove)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

figure_dir <- "Bayesmove_updated_figures"

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. EXTRACT THE UPDATED REFERENCE RESULTS
# ------------------------------------------------------------

reference_fit <- clustering_long_a

reference_theta <- bayesmove::extract_prop(
  res = reference_fit,
  ngibbs = long_iterations,
  nburn = long_burn,
  nmaxclust = max_clusters
) %>%
  as.matrix()

reference_cluster_names <- paste0(
  "cluster_",
  seq_len(max_clusters)
)

colnames(reference_theta) <- reference_cluster_names

stopifnot(
  nrow(reference_theta) == nrow(segment_bin_counts),
  all(is.finite(reference_theta)),
  all(
    abs(rowSums(reference_theta) - 1) < 1e-5
  )
)

reference_segment_props <- bind_cols(
  as_tibble(
    segment_bin_counts[, c("id", "tseg")]
  ),
  as_tibble(reference_theta)
)

# Start from the segmented data, without previously joined
# cluster proportions.
reference_movement <- tracks_segmented %>%
  left_join(
    reference_segment_props,
    by = c("id", "tseg")
  ) %>%
  arrange(track_id, system, date)

stopifnot(
  nrow(reference_movement) == nrow(tracks_segmented),
  !anyNA(
    reference_movement[, reference_cluster_names]
  )
)

# Combine smaller clusters for DISPLAY only.
other_columns <- paste0(
  "cluster_",
  4:max_clusters
)

reference_movement$other <- rowSums(
  as.matrix(
    reference_movement[, other_columns, drop = FALSE]
  )
)

# ------------------------------------------------------------
# 2. SHARED FIGURE SETTINGS
# ------------------------------------------------------------

display_levels <- c(
  "cluster_1",
  "cluster_2",
  "cluster_3",
  "other"
)

display_labels <- c(
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
  display_labels
)

paper_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = "grey95",
      color = "grey60"
    ),
    legend.position = "bottom",
    legend.title = element_blank()
  )

save_paper_figure <- function(
    plot,
    filename,
    width,
    height
) {
  
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(filename, ".png")
    ),
    plot = plot,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(filename, ".pdf")
    ),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
}

# ------------------------------------------------------------
# 3. FIGURE: FITTED MOVEMENT DISTRIBUTIONS
# ------------------------------------------------------------

reference_bins <- bayesmove::get_behav_hist(
  dat = reference_fit,
  nburn = long_burn,
  ngibbs = long_iterations,
  nmaxclust = max_clusters,
  var.names = c(
    "Step length",
    "Turning angle"
  )
) %>%
  mutate(
    cluster_number = as.integer(
      as.character(behav)
    )
  )

# Show the three focal clusters.
# Their fitted bin probabilities are not renormalized.
distribution_data <- reference_bins %>%
  filter(cluster_number %in% 1:3) %>%
  mutate(
    state = factor(
      cluster_number,
      levels = 1:3,
      labels = display_labels[1:3]
    )
  )

distribution_plot <- ggplot(
  distribution_data,
  aes(
    x = bin,
    y = prop,
    fill = state
  )
) +
  geom_col(
    width = 0.8
  ) +
  facet_grid(
    state ~ var,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = state_colors,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = 1:8
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1)
  ) +
  labs(
    x = "Movement-variable bin",
    y = "Fitted bin probability"
  ) +
  paper_theme

print(distribution_plot)

save_paper_figure(
  distribution_plot,
  "Bayesmove_movement_distributions",
  width = 9,
  height = 6
)

# ------------------------------------------------------------
# 4. FIGURE: COMPOSITION THROUGH TIME
# ------------------------------------------------------------

reference_long <- reference_movement %>%
  select(
    track_id,
    system,
    date,
    all_of(display_levels)
  ) %>%
  pivot_longer(
    cols = all_of(display_levels),
    names_to = "cluster",
    values_to = "proportion"
  ) %>%
  mutate(
    state = factor(
      cluster,
      levels = display_levels,
      labels = display_labels
    ),
    
    system_label = factor(
      system,
      levels = c("Argos", "Kineis"),
      labels = c("Argos", "Kin\u00e9is")
    )
  )

time_plot <- ggplot(
  reference_long,
  aes(
    x = date,
    y = proportion,
    fill = state
  )
) +
  geom_col(
    width = 86400 * 0.95
  ) +
  facet_wrap(
    vars(track_id, system_label),
    ncol = 2,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = state_colors,
    drop = FALSE
  ) +
  scale_x_datetime(
    date_breaks = "2 weeks",
    date_labels = "%b %d"
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0%", "50%", "100%")
  ) +
  labs(
    x = "Date (UTC)",
    y = "Estimated segment composition"
  ) +
  paper_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )

print(time_plot)

save_paper_figure(
  time_plot,
  "Bayesmove_composition_through_time",
  width = 11,
  height = 8
)

# ------------------------------------------------------------
# 5. UPDATED INDIVIDUAL COMPOSITION TABLE
# ------------------------------------------------------------

# Keep ALL seven clusters in the full results.
updated_individual_budgets <- reference_movement %>%
  group_by(track_id, system) %>%
  summarise(
    n_movement_rows = n(),
    
    across(
      all_of(reference_cluster_names),
      ~ mean(.x)
    ),
    
    .groups = "drop"
  )

updated_summary <- updated_individual_budgets %>%
  transmute(
    track_id,
    system,
    n_movement_rows,
    
    shorter_broad_percent =
      100 * cluster_1,
    
    longer_mixed_percent =
      100 * cluster_2,
    
    shorter_small_turns_percent =
      100 * cluster_3,
    
    remaining_percent =
      100 * (
        1 - cluster_1 - cluster_2 - cluster_3
      )
  )

print(
  updated_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 6. FIGURE: COMPOSITION BY SHARK AND SYSTEM
# ------------------------------------------------------------

budget_plot_data <- reference_long %>%
  group_by(
    track_id,
    system_label,
    state
  ) %>%
  summarise(
    proportion = mean(proportion),
    .groups = "drop"
  )

budget_plot <- ggplot(
  budget_plot_data,
  aes(
    x = system_label,
    y = proportion,
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
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c(
      "0%", "25%", "50%", "75%", "100%"
    )
  ) +
  labs(
    x = NULL,
    y = "Mean composition across daily movement rows"
  ) +
  paper_theme +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )

print(budget_plot)

save_paper_figure(
  budget_plot,
  "Bayesmove_individual_composition",
  width = 10,
  height = 5
)

# ------------------------------------------------------------
# 7. EXPORT UPDATED RESULTS AND BIN BOUNDARIES
# ------------------------------------------------------------

readr::write_csv(
  updated_individual_budgets,
  file.path(
    figure_dir,
    "Updated_individual_budgets_all_clusters.csv"
  )
)

readr::write_csv(
  updated_summary,
  file.path(
    figure_dir,
    "Updated_individual_composition_percent.csv"
  )
)

readr::write_csv(
  reference_bins,
  file.path(
    figure_dir,
    "Updated_fitted_bin_probabilities.csv"
  )
)

readr::write_csv(
  step_bin_table,
  file.path(
    figure_dir,
    "Step_length_bin_boundaries_km.csv"
  )
)

readr::write_csv(
  angle_bin_table,
  file.path(
    figure_dir,
    "Turning_angle_bin_boundaries_degrees.csv"
  )
)

message(
  "Figures and tables saved in: ",
  normalizePath(figure_dir)
)