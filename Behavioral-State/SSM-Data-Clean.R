# ============================================================
# STEP 1: CREATE DAILY SSM TRACKS FOR BAYESMOVE
#
# Required files in the working directory:
#   Argos-CLEAN_postdeploy_water.csv
#   Doppler-CLEAN_postdeploy_water.csv
#
# Main outputs:
#   Argos-SSM_daily.csv
#   Doppler-SSM_daily.csv
#
# Additional outputs:
#   SSM_track_eligibility.csv
#   SSM_fit_status.csv
#   SSM_invalid_rows.csv
#   SSM_duplicate_times.csv
#   SSM_input_locations.csv
#   SSM_diagnostics/
#
# Model:
#   aniMotum correlated random walk: model = "crw"
#   Prediction interval: 24 hours
#
# Eligibility rules retained from your earlier SSM script:
#   - At least 10 valid, unique positions per track/system
#   - No gap greater than 10 days
#
# IMPORTANT:
#   Kinéis location classes are treated as Argos-equivalent
#   error categories, matching the earlier script.
#   This assumption requires validation.
#
#   Kinéis error radii are NOT used as coordinate standard
#   errors without a justified conversion.
# ============================================================


# ------------------------------------------------------------
# 0. INSTALL PACKAGES ONCE IF NEEDED
# ------------------------------------------------------------

# install.packages(c(
#   "aniMotum",
#   "dplyr",
#   "readr",
#   "lubridate",
#   "tibble"
# ))


# ------------------------------------------------------------
# 1. SETTINGS
# ------------------------------------------------------------

# Change this if your input files are in another folder.
input_dir <- "."

# Save the daily CSVs in this folder.
# Keep this the same as input_dir in the bayesmove script.
output_dir <- input_dir

min_positions <- 10L
max_gap_days <- 10

# aniMotum uses HOURS here.
prediction_interval_hours <- 24

# Retain aniMotum's default speed/angle prefilter.
use_prefilter <- TRUE

# Maximum travel speed used by that prefilter, in m/s.
# Review this biological assumption for your study.
vmax_ms <- 5


# ------------------------------------------------------------
# 2. LOAD PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "aniMotum",
  "dplyr",
  "readr",
  "lubridate",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(tibble)
})

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Create a separate folder for this run's model objects.
run_label <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

diagnostic_dir <- file.path(
  output_dir,
  "SSM_diagnostics",
  run_label
)

dir.create(
  diagnostic_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_output <- function(x, filename) {
  readr::write_csv(
    x,
    file.path(output_dir, filename)
  )
}

writeLines(
  capture.output(sessionInfo()),
  file.path(
    diagnostic_dir,
    "Session_info.txt"
  )
)

warning(
  "Kinéis classes will be treated as Argos-equivalent ",
  "location-error categories. This reproduces the earlier ",
  "class-based workflow; it does not validate that assumption."
)


# ------------------------------------------------------------
# 3. READ THE CLEANED LOCATION FILES
# ------------------------------------------------------------

read_cleaned <- function(
    filename,
    system_name,
    quality_column
) {
  
  path <- file.path(
    input_dir,
    filename
  )
  
  if (!file.exists(path)) {
    stop("Missing input file: ", path)
  }
  
  raw <- readr::read_csv(
    path,
    col_types = readr::cols(
      .default = readr::col_character()
    )
  )
  
  required_columns <- c(
    "track_id",
    "timestamp",
    "lon",
    "lat",
    quality_column
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(raw)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      filename,
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  raw %>%
    transmute(
      input_row = row_number(),
      
      track_id = trimws(
        as.character(track_id)
      ),
      
      date = lubridate::ymd_hms(
        timestamp,
        tz = "UTC",
        quiet = TRUE
      ),
      
      lc = toupper(
        trimws(
          as.character(
            .data[[quality_column]]
          )
        )
      ),
      
      lon = suppressWarnings(
        as.numeric(lon)
      ),
      
      lat = suppressWarnings(
        as.numeric(lat)
      )
    ) %>%
    mutate(
      system = system_name
    )
}

argos <- read_cleaned(
  filename = "Argos-CLEAN_postdeploy_water.csv",
  system_name = "Argos",
  quality_column = "loc_quality"
)

kineis <- read_cleaned(
  filename = "Doppler-CLEAN_postdeploy_water.csv",
  system_name = "Kineis",
  quality_column = "location_class"
)

locations <- bind_rows(
  argos,
  kineis
)


# ------------------------------------------------------------
# 4. CHECK REQUIRED VALUES AND LOCATION CLASSES
# ------------------------------------------------------------

# These are the Argos-type classes used in this workflow.
allowed_classes <- c(
  "3", "2", "1", "0",
  "A", "B", "Z"
)

locations <- locations %>%
  mutate(
    valid_id =
      !is.na(track_id) &
      nzchar(track_id),
    
    valid_time = !is.na(date),
    
    valid_coordinates =
      is.finite(lon) &
      is.finite(lat) &
      abs(lon) <= 180 &
      abs(lat) <= 90,
    
    valid_class =
      !is.na(lc) &
      lc %in% allowed_classes,
    
    valid =
      valid_id &
      valid_time &
      valid_coordinates &
      valid_class
  )

write_output(
  locations %>%
    filter(!valid),
  "SSM_invalid_rows.csv"
)

# Stop for unrecognized classes rather than silently
# assigning them a different measurement-error category.
unknown_classes <- locations %>%
  filter(
    !is.na(lc),
    nzchar(lc),
    !lc %in% allowed_classes
  ) %>%
  distinct(system, lc)

if (nrow(unknown_classes) > 0) {
  
  print(
    unknown_classes,
    n = Inf
  )
  
  stop(
    "Unrecognized location classes. ",
    "Review SSM_invalid_rows.csv before fitting."
  )
}

locations <- locations %>%
  filter(valid) %>%
  select(
    input_row,
    track_id,
    system,
    date,
    lc,
    lon,
    lat
  )

if (nrow(locations) == 0) {
  stop("No valid input locations remain.")
}


# ------------------------------------------------------------
# 5. RESOLVE DUPLICATE TIMESTAMPS
# ------------------------------------------------------------

# This uses the same class ordering for both systems,
# consistent with the class-based assumption in this script.
class_order <- c(
  "3", "2", "1", "0", "A", "B", "Z"
)

ranked_locations <- locations %>%
  mutate(
    class_rank = match(lc, class_order)
  ) %>%
  group_by(
    system,
    track_id,
    date
  ) %>%
  arrange(
    class_rank,
    input_row,
    .by_group = TRUE
  ) %>%
  mutate(
    n_at_time = n(),
    retained = row_number() == 1L
  ) %>%
  ungroup()

# Save all duplicated timestamps, showing which fix was kept.
write_output(
  ranked_locations %>%
    filter(n_at_time > 1),
  "SSM_duplicate_times.csv"
)

# Save discarded fixes separately.
write_output(
  ranked_locations %>%
    filter(!retained),
  "SSM_discarded_duplicate_fixes.csv"
)

message(
  "Removed ",
  sum(!ranked_locations$retained),
  " duplicate-timestamp fixes."
)

# Keep one fix per system, shark, and timestamp.
locations <- ranked_locations %>%
  filter(retained) %>%
  select(
    -class_rank,
    -n_at_time,
    -retained
  ) %>%
  arrange(
    system,
    track_id,
    date
  )

stopifnot(
  anyDuplicated(
    locations[c("system", "track_id", "date")]
  ) == 0L
)
# ------------------------------------------------------------
# 6. APPLY TRACK ELIGIBILITY RULES
# ------------------------------------------------------------

eligibility <- locations %>%
  group_by(
    system,
    track_id
  ) %>%
  summarise(
    n_positions = n(),
    
    start_date = min(date),
    
    end_date = max(date),
    
    largest_gap_days = if (n() >= 2) {
      max(
        as.numeric(
          diff(date),
          units = "days"
        )
      )
    } else {
      0
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    fails_min_positions =
      n_positions < min_positions,
    
    fails_gap =
      largest_gap_days > max_gap_days,
    
    keep_track =
      !fails_min_positions &
      !fails_gap
  )

write_output(
  eligibility,
  "SSM_track_eligibility.csv"
)

print(
  eligibility,
  n = Inf
)

eligible_locations <- locations %>%
  semi_join(
    eligibility %>%
      filter(keep_track),
    
    by = c(
      "system",
      "track_id"
    )
  )

if (nrow(eligible_locations) == 0) {
  stop(
    "No tracks pass the eligibility rules. ",
    "Review SSM_track_eligibility.csv."
  )
}

write_output(
  eligible_locations,
  "SSM_input_locations.csv"
)


# ------------------------------------------------------------
# 7. FUNCTION TO FIT ONE TRACK
# ------------------------------------------------------------

fit_one_track <- function(track_data) {
  
  this_system <- unique(
    track_data$system
  )
  
  this_id <- unique(
    track_data$track_id
  )
  
  message(
    "\nFitting ",
    this_system,
    " track ",
    this_id,
    " with ",
    nrow(track_data),
    " locations..."
  )
  
  model_input <- track_data %>%
    transmute(
      id = as.character(track_id),
      date = date,
      lc = lc,
      lon = lon,
      lat = lat
    ) %>%
    arrange(date)
  
  captured_warnings <- character()
  
  result <- tryCatch(
    
    withCallingHandlers({
      
      fit <- aniMotum::fit_ssm(
        x = model_input,
        
        model = "crw",
        
        time.step = prediction_interval_hours,
        
        spdf = use_prefilter,
        
        vmax = vmax_ms
      )
      
      # Save the full fit before extracting predictions.
      safe_id <- gsub(
        "[^A-Za-z0-9_-]",
        "_",
        this_id
      )
      
      file_prefix <- paste(
        this_system,
        safe_id,
        sep = "_"
      )
      
      saveRDS(
        fit,
        file.path(
          diagnostic_dir,
          paste0(
            file_prefix,
            "_fit.rds"
          )
        )
      )
      
      predicted <- aniMotum::grab(
        fit,
        what = "predicted",
        as_sf = FALSE
      ) %>%
        as_tibble()
      
      required_prediction_columns <- c(
        "date",
        "lon",
        "lat"
      )
      
      if (
        !all(
          required_prediction_columns %in%
          names(predicted)
        )
      ) {
        stop(
          "Predicted output does not contain date, lon, lat."
        )
      }
      
      if (nrow(predicted) == 0) {
        stop("Model returned no predicted locations.")
      }
      
      if (
        anyNA(predicted$date) ||
        any(!is.finite(predicted$lon)) ||
        any(!is.finite(predicted$lat))
      ) {
        stop(
          "Model returned missing or nonfinite predictions. ",
          "Inspect the saved fit."
        )
      }
      
      predicted <- predicted %>%
        mutate(
          track_id = this_id,
          system = this_system
        ) %>%
        select(
          track_id,
          system,
          date,
          lon,
          lat,
          everything()
        ) %>%
        arrange(date)
      
      list(
        predictions = predicted,
        
        status = tibble(
          system = this_system,
          track_id = this_id,
          status = "fit_returned_review_diagnostics",
          n_input = nrow(model_input),
          n_predictions = nrow(predicted),
          error_message = NA_character_
        )
      )
      
    }, warning = function(w) {
      
      captured_warnings <<- c(
        captured_warnings,
        conditionMessage(w)
      )
      
      invokeRestart("muffleWarning")
      
    }),
    
    error = function(e) {
      
      list(
        predictions = NULL,
        
        status = tibble(
          system = this_system,
          track_id = this_id,
          status = "failed",
          n_input = nrow(model_input),
          n_predictions = 0L,
          error_message = conditionMessage(e)
        )
      )
    }
  )
  
  result$status$warnings <- if (
    length(captured_warnings) > 0
  ) {
    paste(
      unique(captured_warnings),
      collapse = " | "
    )
  } else {
    NA_character_
  }
  
  result
}


# ------------------------------------------------------------
# 8. FIT ALL ELIGIBLE TRACKS
# ------------------------------------------------------------

track_keys <- eligible_locations %>%
  distinct(
    system,
    track_id
  ) %>%
  arrange(
    system,
    track_id
  )

results <- vector(
  mode = "list",
  length = nrow(track_keys)
)

for (i in seq_len(nrow(track_keys))) {
  
  current_system <- track_keys$system[i]
  current_id <- track_keys$track_id[i]
  
  current_track <- eligible_locations %>%
    filter(
      system == current_system,
      track_id == current_id
    )
  
  results[[i]] <- fit_one_track(
    current_track
  )
  
  # Save the status table after every track.
  # This preserves progress if a later fit is interrupted.
  current_status <- bind_rows(
    lapply(
      results[seq_len(i)],
      function(z) z$status
    )
  )
  
  write_output(
    current_status,
    "SSM_fit_status.csv"
  )
}


# ------------------------------------------------------------
# 9. COMBINE AND EXPORT DAILY PREDICTIONS
# ------------------------------------------------------------

fit_status <- bind_rows(
  lapply(
    results,
    function(z) z$status
  )
)

daily_predictions <- bind_rows(
  lapply(
    results,
    function(z) z$predictions
  )
)

write_output(
  fit_status,
  "SSM_fit_status.csv"
)

# Define empty outputs too, so a failed system does not
# leave an older CSV looking like this run's result.
empty_predictions <- tibble(
  track_id = character(),
  system = character(),
  date = as.POSIXct(
    character(),
    tz = "UTC"
  ),
  lon = double(),
  lat = double()
)

if (nrow(daily_predictions) == 0) {
  
  argos_daily <- empty_predictions
  kineis_daily <- empty_predictions
  
} else {
  
  argos_daily <- daily_predictions %>%
    filter(system == "Argos") %>%
    arrange(
      track_id,
      date
    )
  
  kineis_daily <- daily_predictions %>%
    filter(system == "Kineis") %>%
    arrange(
      track_id,
      date
    )
}

write_output(
  argos_daily,
  "Argos-SSM_daily.csv"
)

write_output(
  kineis_daily,
  "Doppler-SSM_daily.csv"
)


# ------------------------------------------------------------
# 10. CHECK DAILY PREDICTION INTERVALS
# ------------------------------------------------------------

if (nrow(daily_predictions) > 0) {
  
  interval_check <- daily_predictions %>%
    arrange(
      system,
      track_id,
      date
    ) %>%
    group_by(
      system,
      track_id
    ) %>%
    summarise(
      n_predictions = n(),
      
      minimum_interval_hours = if (n() >= 2) {
        min(
          as.numeric(
            diff(date),
            units = "hours"
          )
        )
      } else {
        NA_real_
      },
      
      maximum_interval_hours = if (n() >= 2) {
        max(
          as.numeric(
            diff(date),
            units = "hours"
          )
        )
      } else {
        NA_real_
      },
      
      .groups = "drop"
    )
  
  write_output(
    interval_check,
    "SSM_prediction_interval_check.csv"
  )
  
  print(
    interval_check,
    n = Inf
  )
}


# ------------------------------------------------------------
# 11. REPORT RESULTS
# ------------------------------------------------------------

print(
  fit_status,
  n = Inf
)

message(
  "\nSaved Argos-SSM_daily.csv: ",
  nrow(argos_daily),
  " predicted locations."
)

message(
  "Saved Doppler-SSM_daily.csv: ",
  nrow(kineis_daily),
  " predicted locations."
)

message(
  "Full model objects saved in: ",
  normalizePath(diagnostic_dir)
)

if (
  nrow(argos_daily) == 0 ||
  nrow(kineis_daily) == 0
) {
  warning(
    "At least one system has no predictions. ",
    "Review SSM_fit_status.csv before running bayesmove."
  )
}

message(
  "\nReview model warnings and saved fits before ",
  "using the predictions for behavioral inference."
)