library(dplyr)
library(lubridate)
library(tibble)

# Exclude 45272 from the paired segmented analysis because
# its shared window contains only three positions per system.
analysis_windows <- shared_windows %>%
  filter(
    has_overlap,
    track_id != "45272"
  ) %>%
  mutate(
    track_id = as.character(track_id)
  )

# Build a shared daily prediction grid at midnight UTC.
# Every requested time stays inside the shared window.
common_prediction_times <- bind_rows(
  lapply(
    seq_len(nrow(analysis_windows)),
    function(i) {
      
      first_time <- ceiling_date(
        analysis_windows$shared_start[i],
        unit = "day",
        change_on_boundary = FALSE
      )
      
      last_time <- floor_date(
        analysis_windows$shared_end[i],
        unit = "day"
      )
      
      if (first_time > last_time) {
        return(
          tibble(
            id = character(),
            date = as.POSIXct(character(), tz = "UTC")
          )
        )
      }
      
      tibble(
        id = analysis_windows$track_id[i],
        date = seq(
          from = first_time,
          to = last_time,
          by = "24 hours"
        )
      )
    }
  )
)

# Summarize the requested grid.
common_grid_summary <- common_prediction_times %>%
  group_by(id) %>%
  summarise(
    n_common_positions = n(),
    first_prediction = min(date),
    last_prediction = max(date),
    
    # Upper bounds before accounting for zero-length steps
    # or other unusable movement values.
    n_possible_steps = n() - 1L,
    n_possible_angles = pmax(n() - 2L, 0L),
    
    .groups = "drop"
  )

print(
  common_grid_summary,
  n = Inf,
  width = Inf
)