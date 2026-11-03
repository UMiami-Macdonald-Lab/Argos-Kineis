# ============================================================
# NEXT STEP: SSM PREDICTIONS ON THE COMMON DAILY GRID
#
# Requires:
#   common_prediction_times from the preceding chunk
#   SSM_input_locations.csv from the earlier SSM workflow
#
# Retains the earlier class-based measurement-error approach.
# ============================================================

library(aniMotum)
library(dplyr)
library(readr)
library(tibble)

stopifnot(
  exists("common_prediction_times"),
  inherits(common_prediction_times$date, "POSIXct")
)

# Use the folder containing SSM_input_locations.csv.
ssm_folder <- "."

# Match the settings used for the original SSM fits.
aligned_prefilter <- TRUE
aligned_vmax <- 5

# Read the original eligible observations, NOT daily predictions.
ssm_inputs <- read_csv(
  file.path(ssm_folder, "SSM_input_locations.csv"),
  col_types = cols(
    track_id = col_character(),
    system = col_character(),
    lc = col_character(),
    date = col_datetime()
  ),
  show_col_types = FALSE
)

# Separate folder for this run's model objects.
aligned_run_dir <- file.path(
  ssm_folder,
  paste0(
    "SSM_aligned_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

dir.create(
  aligned_run_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Every retained shark needs one fit from each system.
fit_keys <- expand.grid(
  track_id = unique(common_prediction_times$id),
  system = c("Argos", "Kineis"),
  stringsAsFactors = FALSE
)

aligned_results <- vector(
  "list",
  nrow(fit_keys)
)

for (i in seq_len(nrow(fit_keys))) {
  
  this_id <- fit_keys$track_id[i]
  this_system <- fit_keys$system[i]
  
  # Retain the full eligible observation series for fitting.
  model_input <- ssm_inputs %>%
    filter(
      track_id == this_id,
      system == this_system
    ) %>%
    transmute(
      id = as.character(track_id),
      date,
      lc,
      lon,
      lat
    ) %>%
    arrange(date)
  
  prediction_schedule <- common_prediction_times %>%
    filter(id == this_id) %>%
    select(id, date) %>%
    arrange(date) %>%
    as.data.frame()
  
  message(
    "Fitting ", this_system, " shark ", this_id,
    "; requesting ", nrow(prediction_schedule),
    " aligned predictions."
  )
  
  fit_warnings <- character()
  
  aligned_results[[i]] <- tryCatch(
    
    withCallingHandlers({
      
      if (nrow(model_input) == 0) {
        stop("No input observations found.")
      }
      
      fit <- aniMotum::fit_ssm(
        x = model_input,
        model = "crw",
        time.step = prediction_schedule,
        spdf = aligned_prefilter,
        vmax = aligned_vmax
      )
      
      saveRDS(
        fit,
        file.path(
          aligned_run_dir,
          paste0(this_system, "_", this_id, "_fit.rds")
        )
      )
      
      predicted <- aniMotum::grab(
        fit,
        what = "predicted",
        as_sf = FALSE
      ) %>%
        as_tibble() %>%
        mutate(
          track_id = this_id,
          system = this_system
        ) %>%
        arrange(date)
      
      # Verify every requested time was returned exactly once.
      if (
        nrow(predicted) != nrow(prediction_schedule) ||
        anyNA(predicted$date) ||
        anyDuplicated(predicted$date) > 0
      ) {
        stop("Returned prediction count/timestamps need review.")
      }
      
      time_difference <- abs(
        as.numeric(predicted$date) -
          as.numeric(prediction_schedule$date)
      )
      
      if (any(time_difference > 0.001)) {
        stop("Returned timestamps do not match the requested grid.")
      }
      
      if (
        any(!is.finite(predicted$lon)) ||
        any(!is.finite(predicted$lat))
      ) {
        stop("Nonfinite predicted coordinates.")
      }
      
      list(
        predictions = predicted,
        error = NA_character_
      )
      
    }, warning = function(w) {
      
      fit_warnings <<- c(
        fit_warnings,
        conditionMessage(w)
      )
      
      invokeRestart("muffleWarning")
      
    }),
    
    error = function(e) {
      list(
        predictions = NULL,
        error = conditionMessage(e)
      )
    }
  )
  
  aligned_results[[i]]$status <- tibble(
    track_id = this_id,
    system = this_system,
    n_input = nrow(model_input),
    n_requested = nrow(prediction_schedule),
    n_returned = nrow(
      aligned_results[[i]]$predictions
    ) %||% 0L,
    error_message = aligned_results[[i]]$error,
    warnings = if (length(fit_warnings)) {
      paste(unique(fit_warnings), collapse = " | ")
    } else {
      NA_character_
    }
  )
}


aligned_status <- bind_rows(
  lapply(aligned_results, function(z) z$status)
)

print(
  aligned_status,
  n = Inf,
  width = Inf
)

write_csv(
  aligned_status,
  file.path(aligned_run_dir, "Aligned_SSM_status.csv")
)

# Export prediction files only if all eight fits passed
# the prediction-count, timestamp, and coordinate checks.
if (all(is.na(aligned_status$error_message))) {
  
  aligned_daily <- bind_rows(
    lapply(aligned_results, function(z) z$predictions)
  ) %>%
    select(
      track_id,
      system,
      date,
      lon,
      lat,
      everything()
    ) %>%
    arrange(track_id, system, date)
  
  write_csv(
    filter(aligned_daily, system == "Argos"),
    file.path(ssm_folder, "Argos-SSM_daily_aligned.csv")
  )
  
  write_csv(
    filter(aligned_daily, system == "Kineis"),
    file.path(ssm_folder, "Doppler-SSM_daily_aligned.csv")
  )
  
  message("Both aligned prediction files were saved.")
  
} else {
  
  message(
    "Some fits failed the checks. ",
    "Review aligned_status before continuing."
  )
}