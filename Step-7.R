# ============================================================
# INSPECT MOVEMENT DISTRIBUTIONS BEFORE DEFINING BINS
# ============================================================

library(dplyr)
library(ggplot2)

# Keep rows with complete, usable movement variables.
movement_for_bayesmove <- movement_metrics %>%
  filter(usable_movement) %>%
  mutate(
    step_km = step / 1000,
    
    system_label = factor(
      system,
      levels = c("Argos", "Kineis"),
      labels = c("Argos", "Kin\u00e9is")
    )
  ) %>%
  arrange(id, date)

# ------------------------------------------------------------
# 1. NUMERICAL SUMMARY BY SHARK AND SYSTEM
# ------------------------------------------------------------

movement_summary <- movement_for_bayesmove %>%
  group_by(track_id, system) %>%
  summarise(
    n_movement_rows = n(),
    
    min_step_km = min(step_km),
    
    q25_step_km = unname(
      quantile(step_km, 0.25)
    ),
    
    median_step_km = median(step_km),
    
    q75_step_km = unname(
      quantile(step_km, 0.75)
    ),
    
    q90_step_km = unname(
      quantile(step_km, 0.90)
    ),
    
    max_step_km = max(step_km),
    
    median_abs_angle_deg = median(
      abs(angle) * 180 / pi
    ),
    
    .groups = "drop"
  )

print(
  movement_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 2. STEP-LENGTH DISTRIBUTIONS
# ------------------------------------------------------------

# Common histogram boundaries allow comparison across panels.
step_breaks <- seq(
  from = 0,
  to = max(movement_for_bayesmove$step_km),
  length.out = 26
)

step_plot <- ggplot(
  movement_for_bayesmove,
  aes(
    x = step_km,
    y = after_stat(density),
    fill = system_label
  )
) +
  geom_histogram(
    breaks = step_breaks,
    color = "white",
    linewidth = 0.2
  ) +
  facet_grid(
    track_id ~ system_label
  ) +
  scale_fill_manual(
    values = c(
      "Argos" = "#5FAF89",
      "Kin\u00e9is" = "#62B7C8"
    ),
    guide = "none"
  ) +
  labs(
    x = "Daily step length (km)",
    y = "Density"
  ) +
  theme_bw(base_size = 11)

print(step_plot)

# ------------------------------------------------------------
# 3. TURNING-ANGLE DISTRIBUTIONS
# ------------------------------------------------------------

angle_plot <- ggplot(
  movement_for_bayesmove,
  aes(
    x = angle,
    y = after_stat(density),
    fill = system_label
  )
) +
  geom_histogram(
    breaks = seq(
      -pi,
      pi,
      length.out = 17
    ),
    color = "white",
    linewidth = 0.2
  ) +
  facet_grid(
    track_id ~ system_label
  ) +
  scale_fill_manual(
    values = c(
      "Argos" = "#5FAF89",
      "Kin\u00e9is" = "#62B7C8"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(
      -pi,
      -pi / 2,
      0,
      pi / 2,
      pi
    ),
    labels = c(
      "-180",
      "-90",
      "0",
      "90",
      "180"
    )
  ) +
  labs(
    x = "Turning angle (degrees)",
    y = "Density"
  ) +
  theme_bw(base_size = 11)

print(angle_plot)