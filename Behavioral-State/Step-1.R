library(readr)
library(dplyr)

# 1. Check which SSM fits returned results.
ssm_status <- read_csv(
  "SSM_fit_status.csv",
  show_col_types = FALSE
)

print(
  ssm_status,
  n = Inf,
  width = Inf
)

# 2. Check the prediction intervals.
interval_check <- read_csv(
  "SSM_prediction_interval_check.csv",
  show_col_types = FALSE
)

print(
  interval_check,
  n = Inf,
  width = Inf
)

# 3. Check how many sharks have predictions from both systems.
argos_daily <- read_csv(
  "Argos-SSM_daily.csv",
  col_types = cols(track_id = col_character()),
  show_col_types = FALSE
)

kineis_daily <- read_csv(
  "Doppler-SSM_daily.csv",
  col_types = cols(track_id = col_character()),
  show_col_types = FALSE
)

track_availability <- full_join(
  argos_daily %>%
    count(track_id, name = "n_Argos"),
  
  kineis_daily %>%
    count(track_id, name = "n_Kineis"),
  
  by = "track_id"
) %>%
  mutate(
    both_systems =
      !is.na(n_Argos) &
      !is.na(n_Kineis)
  ) %>%
  arrange(track_id)

print(
  track_availability,
  n = Inf,
  width = Inf
)