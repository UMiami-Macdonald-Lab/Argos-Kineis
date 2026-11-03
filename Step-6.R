# ============================================================
# CALCULATE MOVEMENT VARIABLES FOR THE THREE PAIRED SHARKS
# ============================================================

library(dplyr)
library(readr)
library(lubridate)
library(sf)
library(bayesmove)

# Change if your aligned CSVs are in another folder.
data_dir <- "."

working_ids <- c(
  "34885",
  "34890",
  "45270"
)

read_aligned <- function(filename, system_name) {
  
  read_csv(
    file.path(data_dir, filename),
    col_types = cols(
      track_id = col_character(),
      date = col_datetime()
    ),
    show_col_types = FALSE
  ) %>%
    transmute(
      track_id,
      system = system_name,
      date,
      lon,
      lat
    )
}

movement_positions <- bind_rows(
  read_aligned(
    "Argos-SSM_daily_aligned.csv",
    "Argos"
  ),
  
  read_aligned(
    "Doppler-SSM_daily_aligned.csv",
    "Kineis"
  )
) %>%
  filter(track_id %in% working_ids) %>%
  mutate(
    # Separate shark × system series.
    id = paste(track_id, system, sep = "__")
  ) %>%
  arrange(id, date)

stopifnot(
  !anyNA(movement_positions$date),
  all(is.finite(movement_positions$lon)),
  all(is.finite(movement_positions$lat)),
  anyDuplicated(
    movement_positions[c("id", "date")]
  ) == 0L
)

# ------------------------------------------------------------
# CONFIRM UNINTERRUPTED DAILY SPACING
# ------------------------------------------------------------

spacing_check <- movement_positions %>%
  arrange(id, date) %>%
  group_by(id) %>%
  summarise(
    all_daily = all(
      abs(
        as.numeric(diff(date), units = "secs") - 86400
      ) < 1
    ),
    .groups = "drop"
  )

print(
  spacing_check,
  n = Inf,
  width = Inf
)

stopifnot(
  all(spacing_check$all_daily)
)

# ------------------------------------------------------------
# PROJECT COORDINATES INTO METERS: UTM ZONE 17N
# ------------------------------------------------------------

# Review projection suitability if tracks extend far outside
# South Florida and adjacent waters.

projected <- movement_positions %>%
  st_as_sf(
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  ) %>%
  st_transform(
    crs = 32617
  )

xy <- st_coordinates(projected)

movement_positions$x <- xy[, 1]
movement_positions$y <- xy[, 2]

# ------------------------------------------------------------
# CALCULATE MOVEMENT FOR EACH SHARK × SYSTEM
# ------------------------------------------------------------

movement_metrics <- bind_rows(
  lapply(
    split(
      as.data.frame(movement_positions),
      movement_positions$id
    ),
    
    function(z) {
      
      bayesmove::prep_data(
        dat = z,
        coord.names = c("x", "y"),
        id = "id"
      )
    }
  )
) %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(
    previous_step = lag(step),
    
    # A turning angle needs both adjacent steps to have
    # a defined direction.
    usable_movement =
      is.finite(step) &
      is.finite(angle) &
      is.finite(dt) &
      abs(dt - 86400) < 1 &
      step > 0 &
      !is.na(previous_step) &
      previous_step > 0
  ) %>%
  ungroup()

# ------------------------------------------------------------
# SUMMARIZE AVAILABILITY BEFORE REMOVING ANY ROWS
# ------------------------------------------------------------

movement_check <- movement_metrics %>%
  group_by(track_id, system) %>%
  summarise(
    n_positions = n(),
    
    n_finite_steps = sum(
      is.finite(step)
    ),
    
    n_finite_angles = sum(
      is.finite(angle)
    ),
    
    n_zero_steps = sum(
      step == 0,
      na.rm = TRUE
    ),
    
    n_usable_movement_rows = sum(
      usable_movement
    ),
    
    .groups = "drop"
  )

print(
  movement_check,
  n = Inf,
  width = Inf
)