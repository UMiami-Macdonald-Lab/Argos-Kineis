library(dplyr)
library(lubridate)

# Ensure timestamps are interpreted as UTC.
parse_utc <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(lubridate::with_tz(x, "UTC"))
  }
  
  lubridate::ymd_hms(
    x,
    tz = "UTC",
    quiet = TRUE
  )
}

# Combine the two prediction datasets.
daily_tracks <- bind_rows(
  argos_daily %>%
    transmute(
      track_id = as.character(track_id),
      system = "Argos",
      date = parse_utc(date),
      lon,
      lat
    ),
  
  kineis_daily %>%
    transmute(
      track_id = as.character(track_id),
      system = "Kineis",
      date = parse_utc(date),
      lon,
      lat
    )
)

stopifnot(!anyNA(daily_tracks$date))

# Identify sharks with predictions from both systems.
paired_ids <- daily_tracks %>%
  distinct(track_id, system) %>%
  count(track_id) %>%
  filter(n == 2) %>%
  select(track_id)

# Find the intersection of their tracking periods.
shared_windows <- daily_tracks %>%
  semi_join(paired_ids, by = "track_id") %>%
  group_by(track_id, system) %>%
  summarise(
    first_prediction = min(date),
    last_prediction = max(date),
    .groups = "drop"
  ) %>%
  group_by(track_id) %>%
  summarise(
    shared_start = max(first_prediction),
    shared_end = min(last_prediction),
    .groups = "drop"
  ) %>%
  mutate(
    overlap_days = as.numeric(
      difftime(
        shared_end,
        shared_start,
        units = "days"
      )
    ),
    has_overlap = shared_end > shared_start
  )

print(
  shared_windows,
  n = Inf,
  width = Inf
)

# Count predictions within each shared window.
overlap_counts <- daily_tracks %>%
  inner_join(
    shared_windows %>% filter(has_overlap),
    by = "track_id"
  ) %>%
  filter(
    date >= shared_start,
    date <= shared_end
  ) %>%
  group_by(track_id, system) %>%
  summarise(
    n_daily_positions = n(),
    first_retained = min(date),
    last_retained = max(date),
    .groups = "drop"
  )

print(
  overlap_counts,
  n = Inf,
  width = Inf
)