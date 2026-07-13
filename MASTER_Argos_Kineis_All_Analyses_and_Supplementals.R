# =============================================================================
# MASTER ARGOS–KINEIS ANALYSIS PIPELINE
# =============================================================================
# This script consolidates the analyses used in the manuscript and supplements:
#   1. Data cleaning and deployment filtering
#   2. 24-hour aniMotum correlated-random-walk state-space models
#   3. Individual-level 50% and 90% KDE polygons and areas
#   4. Individual KDE mixed models, post-hoc tests, diagnostics, and supplements
#   5. Individual Jaccard overlap analyses
#   6. Transmission spacing, clustering, and transmission-yield models
#   7. Location-class distributions and continuous-error comparisons
#   8. Species-level pooled KDE areas and Jaccard overlap
#   9. Behavioral-state matched-day comparison, when BayesMove outputs exist
#  10. Supplemental output inventory and session information
#
# Place this script in the same working directory as:
#   Argos-Sharks.xls
#   Doppler-Sharks.xlsx
#   Detailed_Florida_State_Boundary.shp plus .dbf/.shx/.prj files
#
# Optional behavioral input files:
#   Argos_BayesMove_states.csv
#   Kineis_BayesMove_states.csv
#
# IMPORTANT:
# - Update deploy_tbl in Section 1 if additional tags are added.
# - The script writes analysis products to the current working directory.
# =============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260713)

RUN_BEHAVIOR_IF_AVAILABLE <- TRUE

required_packages <- c(
  "dplyr","tidyr","readr","readxl","writexl","lubridate","tibble",
  "sf","sp","adehabitatHR","purrr","ggplot2","patchwork","scales",
  "lme4","lmerTest","emmeans","multcomp","multcompView","effectsize",
  "performance","MASS","rstatix","aniMotum"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the master script:\n  ",
    paste(missing_packages, collapse = ", "),
    "\n\nExample:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(readxl)
  library(writexl); library(lubridate); library(tibble); library(sf)
  library(sp); library(adehabitatHR); library(purrr); library(ggplot2)
  library(patchwork); library(scales); library(lme4); library(lmerTest)
  library(emmeans); library(multcomp); library(multcompView)
  library(effectsize); library(performance); library(MASS); library(rstatix)
  library(aniMotum)
})

dir.create("Main_Figures", showWarnings = FALSE)
dir.create("Supplemental_Figures", showWarnings = FALSE)
dir.create("Supplemental_Tables", showWarnings = FALSE)
dir.create("Model_Objects", showWarnings = FALSE)

safe_write_csv <- function(x, filename) {
  readr::write_csv(x, file.path("Supplemental_Tables", filename))
}

save_model <- function(x, filename) {
  saveRDS(x, file.path("Model_Objects", filename))
}

cat("\nStarting master Argos–Kineis analysis pipeline...\n")


# ============================================================================
# 01_DATA_CLEANING
# ============================================================================
# ==========================================================
# MAKE TWO CLEAN EXPORT FILES (WITH QUALITY COLUMNS):
#   1) Argos-CLEAN_postdeploy_water.xlsx (and .csv)
#      - includes loc_quality
#      - filters OUT loc_quality == "Z"
#   2) Doppler-CLEAN_postdeploy_water.xlsx (and .csv)
#      - includes doppler_error_m + location_class
#
# RULES (BOTH):
# - keep only deployed track_ids in deploy_tbl
# - drop NA id/time/lat/lon
# - keep only post-deployment
#   (timestamp >= deployment_date for that animal)
# - remove points that fall on Florida land polygon
#   using Detailed_Florida_State_Boundary.shp
#   after transforming to EPSG:3086
#
# NOTES:
# - nickname is NOT expected in the raw files; it is joined
#   from deploy_tbl
# - all dplyr verbs are explicitly namespaced to avoid
#   masking issues
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(sf)
  library(lubridate)
  library(tibble)
})

# optional, for xlsx output
if (!requireNamespace("writexl", quietly = TRUE)) {
}

# only load aniMotum if you truly need it elsewhere
if (!requireNamespace("aniMotum", quietly = TRUE)) {
}
library(aniMotum)


library(dplyr)
library(readxl)
library(sf)
library(lubridate)
library(tibble)
library(writexl) 
# ----------------------------------------------------------
# KEY USER PATH FILES
# ----------------------------------------------------------
argos_path   <- "Argos-Sharks.xls"
doppler_path <- "Doppler-Sharks.xlsx"

# Florida shapefile base map / land mask
florida_shp <- "Detailed_Florida_State_Boundary.shp"

if (!file.exists(argos_path)) {
  stop("Missing Argos file: ", argos_path)
}
if (!file.exists(doppler_path)) {
  stop("Missing Doppler file: ", doppler_path)
}
if (!file.exists(florida_shp)) {
  stop("Missing Florida shapefile: ", florida_shp,
       "\nMake sure the .shp file is in your working directory, along with the matching .shx, .dbf, and .prj files.")
}

# ----------------------------------------------------------
# DEPLOYMENT DATA
# ----------------------------------------------------------
deploy_tbl <- tibble::tribble(
  ~track_id, ~species_code, ~nickname,             ~deployment_date,      ~TL_cm,
  34675,     "CLIM",        "Blacktip-B",          as.Date("2025-11-13"), 169,
  34893,     "CLIM",        "Blacktip-A",          as.Date("2025-11-12"), 182,
  34884,     "GCUV",        "Tiger-B",             as.Date("2025-10-13"), 254,
  34885,     "GCUV",        "Tiger-A",             as.Date("2025-10-02"), 236,
  34890,     "GCUV",        "Tiger-C",             as.Date("2025-10-28"), 204,
  34892,     "GCUV",        "Tiger-D",             as.Date("2025-11-06"), 166,
  34894,     "SMOK",        "GreatHammerhead-A",   as.Date("2025-11-22"), 286,
  34895,     "SMOK",        "GreatHammerhead-B",   as.Date("2025-12-04"), 282
)

# ----------------------------------------------------------
# READ IN RAW DATA
# ----------------------------------------------------------
read_argos <- function(path) {
  x <- readxl::read_excel(path)
  
  needed <- c("Platform ID No.", "Loc. date", "Latitude", "Longitude", "Loc. quality")
  miss <- setdiff(needed, names(x))
  if (length(miss) > 0) {
    stop("Argos file missing columns: ", paste(miss, collapse = ", "))
  }
  
  x %>%
    dplyr::transmute(
      source      = "Argos",
      track_id    = suppressWarnings(as.numeric(`Platform ID No.`)),
      timestamp   = lubridate::as_datetime(`Loc. date`, tz = "UTC"),
      lat         = suppressWarnings(as.numeric(Latitude)),
      lon         = suppressWarnings(as.numeric(Longitude)),
      loc_quality = as.character(`Loc. quality`)
    )
}

read_doppler <- function(path) {
  x <- readxl::read_excel(path)
  
  needed <- c(
    "Device ID", "Location date (UTC)", "Latitude", "Longitude",
    "Doppler Error radius (m)", "Location class"
  )
  miss <- setdiff(needed, names(x))
  if (length(miss) > 0) {
    stop("Doppler file missing columns: ", paste(miss, collapse = ", "))
  }
  
  x %>%
    dplyr::transmute(
      source           = "Doppler",
      track_id         = suppressWarnings(as.numeric(`Device ID`)),
      timestamp        = lubridate::as_datetime(`Location date (UTC)`, tz = "UTC"),
      lat              = suppressWarnings(as.numeric(Latitude)),
      lon              = suppressWarnings(as.numeric(Longitude)),
      doppler_error_m  = suppressWarnings(as.numeric(`Doppler Error radius (m)`)),
      location_class   = as.character(`Location class`)
    )
}

argos_raw <- read_argos(argos_path) %>%
  dplyr::filter(!is.na(loc_quality), loc_quality != "Z")

doppler_raw <- read_doppler(doppler_path)

# ----------------------------------------------------------
# FLORIDA LAND POLYGON FOR LAND FILTER
# Uses user shapefile as the default basemap / land mask
# ----------------------------------------------------------
florida <- sf::st_read(florida_shp, quiet = TRUE)

# validate geometry and combine to one land mask
florida_land <- florida %>%
  sf::st_make_valid() %>%
  sf::st_union() %>%
  sf::st_as_sf()

# transform to projected CRS for land intersection tests
florida_land_3086 <- florida_land %>%
  sf::st_transform(3086)

# ----------------------------------------------------------
# CLEANING FUNCTION
# ----------------------------------------------------------
clean_one_source <- function(df_src, deploy_tbl, florida_land_3086) {
  
  # keep only IDs in deployment table and attach metadata
  df0 <- df_src %>%
    dplyr::filter(
      !is.na(track_id),
      !is.na(timestamp),
      !is.na(lat),
      !is.na(lon)
    ) %>%
    dplyr::semi_join(deploy_tbl, by = "track_id") %>%
    dplyr::left_join(
      deploy_tbl %>%
        dplyr::select(track_id, species_code, nickname, deployment_date, TL_cm),
      by = "track_id"
    ) %>%
    dplyr::filter(as.Date(timestamp) >= deployment_date) %>%
    dplyr::arrange(track_id, timestamp)
  
  if (nrow(df0) == 0) {
    return(df0 %>% dplyr::mutate(on_land = logical(0)))
  }
  
  # convert to sf points in WGS84, keep lon/lat columns
  pts_wgs <- sf::st_as_sf(
    df0,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  # transform points to projected CRS for land filtering
  pts_3086 <- pts_wgs %>%
    sf::st_transform(3086)
  
  # identify points intersecting Florida land polygon
  on_land_mat <- sf::st_intersects(pts_3086, florida_land_3086, sparse = FALSE)
  on_land <- apply(on_land_mat, 1, any)
  
  pts_water <- pts_3086 %>%
    dplyr::mutate(on_land = on_land) %>%
    dplyr::filter(!on_land)
  
  if (nrow(pts_water) == 0) {
    out <- df0[0, ] %>% dplyr::mutate(on_land = logical(0))
    return(out)
  }
  
  pts_water %>%
    sf::st_drop_geometry() %>%
    dplyr::select(
      source,
      track_id,
      nickname,
      species_code,
      TL_cm,
      deployment_date,
      timestamp,
      lat,
      lon,
      dplyr::any_of(c("loc_quality", "doppler_error_m", "location_class")),
      on_land
    )
}

# ----------------------------------------------------------
# CLEAN EACH SOURCE
# ----------------------------------------------------------
argos_clean <- clean_one_source(
  df_src = argos_raw,
  deploy_tbl = deploy_tbl,
  florida_land_3086 = florida_land_3086
)

doppler_clean <- clean_one_source(
  df_src = doppler_raw,
  deploy_tbl = deploy_tbl,
  florida_land_3086 = florida_land_3086
)

# ----------------------------------------------------------
# EXPORT
# ----------------------------------------------------------
utils::write.csv(
  argos_clean,
  "Argos-CLEAN_postdeploy_water.csv",
  row.names = FALSE
)

utils::write.csv(
  doppler_clean,
  "Doppler-CLEAN_postdeploy_water.csv",
  row.names = FALSE
)

writexl::write_xlsx(
  argos_clean,
  "Argos-CLEAN_postdeploy_water.xlsx"
)

writexl::write_xlsx(
  doppler_clean,
  "Doppler-CLEAN_postdeploy_water.xlsx"
)

message("Wrote: Argos-CLEAN_postdeploy_water.(csv/xlsx) with ", nrow(argos_clean), " rows")
message("Wrote: Doppler-CLEAN_postdeploy_water.(csv/xlsx) with ", nrow(doppler_clean), " rows")

# ----------------------------------------------------------
# QUICK CHECKS ON DATA
# ----------------------------------------------------------

# A) Raw counts (nickname not present yet)
cat("\nRAW COUNTS (no nickname in raw files):\n")
argos_raw %>% dplyr::count(track_id, sort = TRUE) %>% print(n = Inf)
doppler_raw %>% dplyr::count(track_id, sort = TRUE) %>% print(n = Inf)

# B) Raw counts WITH nickname (join from deploy_tbl)
cat("\nRAW COUNTS (after joining nickname from deploy_tbl):\n")

argos_raw %>%
  dplyr::semi_join(deploy_tbl, by = "track_id") %>%
  dplyr::left_join(
    deploy_tbl %>% dplyr::select(track_id, nickname),
    by = "track_id"
  ) %>%
  dplyr::count(track_id, nickname, sort = TRUE) %>%
  print(n = Inf)

doppler_raw %>%
  dplyr::semi_join(deploy_tbl, by = "track_id") %>%
  dplyr::left_join(
    deploy_tbl %>% dplyr::select(track_id, nickname),
    by = "track_id"
  ) %>%
  dplyr::count(track_id, nickname, sort = TRUE) %>%
  print(n = Inf)

# C) Final cleaned counts
cat("\nCLEANED COUNTS (post-deploy + water-only):\n")
argos_clean %>% dplyr::count(track_id, nickname, sort = TRUE) %>% print(n = Inf)
doppler_clean %>% dplyr::count(track_id, nickname, sort = TRUE) %>% print(n = Inf)

# D) Confirm Argos has no Z left
if ("loc_quality" %in% names(argos_clean)) {
  cat("\nARGOS loc_quality distribution (cleaned):\n")
  argos_clean %>% dplyr::count(loc_quality, sort = TRUE) %>% print(n = Inf)
}

# E) Doppler quality summaries
if ("location_class" %in% names(doppler_clean)) {
  cat("\nDOPPLER location_class distribution (cleaned):\n")
  doppler_clean %>% dplyr::count(location_class, sort = TRUE) %>% print(n = Inf)
}

if ("doppler_error_m" %in% names(doppler_clean)) {
  cat("\nDOPPLER error radius summary (m):\n")
  print(summary(doppler_clean$doppler_error_m))
}

# F) Optional: confirm shapefile CRS and geometry status
cat("\nFLORIDA SHAPEFILE CRS:\n")
print(sf::st_crs(florida))

cat("\nFLORIDA LAND MASK SUMMARY:\n")
print(florida_land_3086)

# ============================================================================
# 02_STATE_SPACE_MODELS
# ============================================================================
# ==========================================================
# aniMotum SSM REGULARIZATION TO DAILY (24-HOUR) TIMESTEP
#
# Exclusion rules (PER track_id, PER source):
#   - Exclude if < 10 positions
#   - Exclude if any gap > 10 days
#
# Outputs:
#   Argos-BEFORE_SSM.csv
#   Doppler-BEFORE_SSM.csv
#   Argos-SSM_daily.csv
#   Doppler-SSM_daily.csv
#   SSM_excluded_tracks.csv
#
# NOTE:
#   In aniMotum, time.step is in HOURS.
#   So daily = time.step = 24
# ==========================================================


suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
  library(tibble)
})

library(aniMotum)


# -------------------------
# 0) INPUT FILES
# -------------------------
argos_in   <- "Argos-CLEAN_postdeploy_water.csv"
doppler_in <- "Doppler-CLEAN_postdeploy_water.csv"

if (!file.exists(argos_in))   stop("Missing: ", argos_in)
if (!file.exists(doppler_in)) stop("Missing: ", doppler_in)

argos_df   <- readr::read_csv(argos_in, show_col_types = FALSE)
doppler_df <- readr::read_csv(doppler_in, show_col_types = FALSE)

# -------------------------
# 1) HELPER: safely parse timestamps
# -------------------------
parse_time_utc <- function(x) {
  if (inherits(x, "POSIXct")) {
    return(x)
  }
  suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
}

# -------------------------
# 2) TRACK FILTERING
# -------------------------
summarize_track_filters <- function(df,
                                    id_col = "track_id",
                                    time_col = "timestamp",
                                    min_n = 10,
                                    max_gap_days_limit = 10) {
  
  df %>%
    dplyr::mutate(.time = parse_time_utc(.data[[time_col]])) %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(.time)) %>%
    dplyr::arrange(.data[[id_col]], .time) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      n_pos = dplyr::n(),
      max_gap_days = dplyr::if_else(
        dplyr::n() < 2,
        0,
        max(as.numeric(difftime(.time[-1], .time[-dplyr::n()], units = "days")), na.rm = TRUE)
      ),
      fails_min_n = n_pos < min_n,
      fails_gap   = max_gap_days > max_gap_days_limit,
      keep_track  = !(fails_min_n | fails_gap),
      .groups = "drop"
    ) %>%
    dplyr::rename(track_id = 1)
}

argos_filter_tbl   <- summarize_track_filters(argos_df,   min_n = 10, max_gap_days_limit = 10)
doppler_filter_tbl <- summarize_track_filters(doppler_df, min_n = 10, max_gap_days_limit = 10)

excluded <- dplyr::bind_rows(
  argos_filter_tbl   %>% dplyr::mutate(source = "Argos"),
  doppler_filter_tbl %>% dplyr::mutate(source = "Doppler")
) %>%
  dplyr::filter(!keep_track) %>%
  dplyr::arrange(source, track_id)

readr::write_csv(excluded, "SSM_excluded_tracks.csv")

# -------------------------
# 3) "BEFORE" DATA = exact subset used for SSM
# -------------------------
argos_before <- argos_df %>%
  dplyr::semi_join(
    argos_filter_tbl %>% dplyr::filter(keep_track),
    by = "track_id"
  ) %>%
  dplyr::mutate(timestamp = parse_time_utc(timestamp)) %>%
  dplyr::filter(!is.na(timestamp), !is.na(lon), !is.na(lat), !is.na(loc_quality)) %>%
  dplyr::arrange(track_id, timestamp)

doppler_before <- doppler_df %>%
  dplyr::semi_join(
    doppler_filter_tbl %>% dplyr::filter(keep_track),
    by = "track_id"
  ) %>%
  dplyr::mutate(timestamp = parse_time_utc(timestamp)) %>%
  dplyr::filter(!is.na(timestamp), !is.na(lon), !is.na(lat), !is.na(location_class)) %>%
  dplyr::arrange(track_id, timestamp)

readr::write_csv(argos_before,   "Argos-BEFORE_SSM.csv")
readr::write_csv(doppler_before, "Doppler-BEFORE_SSM.csv")

# -------------------------
# 4) FORMAT FOR aniMotum
# required columns: id, date, lc, lon, lat
# -------------------------
argos_am <- argos_before %>%
  dplyr::transmute(
    id   = as.character(track_id),
    date = timestamp,
    lc   = as.character(loc_quality),
    lon  = as.numeric(lon),
    lat  = as.numeric(lat)
  ) %>%
  dplyr::filter(!is.na(id), !is.na(date), !is.na(lc), !is.na(lon), !is.na(lat)) %>%
  dplyr::arrange(id, date)

doppler_am <- doppler_before %>%
  dplyr::transmute(
    id   = as.character(track_id),
    date = timestamp,
    lc   = as.character(location_class),
    lon  = as.numeric(lon),
    lat  = as.numeric(lat)
  ) %>%
  dplyr::filter(!is.na(id), !is.na(date), !is.na(lc), !is.na(lon), !is.na(lat)) %>%
  dplyr::arrange(id, date)

# -------------------------
# 5) FIT SSM
# -------------------------
argos_fit <- aniMotum::fit_ssm(
  argos_am,
  model = "crw",
  time.step = 24
)

doppler_fit <- aniMotum::fit_ssm(
  doppler_am,
  model = "crw",
  time.step = 24
)

# -------------------------
# 6) HELPER: extract predicted locations safely
# -------------------------
extract_predicted <- function(fit_obj) {
  
  out <- aniMotum::grab(fit_obj, what = "predicted") %>%
    tibble::as_tibble()
  
  # create numeric track_id from id if present
  if ("id" %in% names(out)) {
    out <- out %>%
      dplyr::mutate(track_id = suppressWarnings(as.numeric(id)))
  }
  
  # expected key columns first, then everything else
  first_cols <- c("track_id", "id", "date", "lon", "lat")
  first_cols <- first_cols[first_cols %in% names(out)]
  
  out %>%
    dplyr::select(
      dplyr::all_of(first_cols),
      dplyr::everything()
    ) %>%
    dplyr::distinct()
}

argos_daily   <- extract_predicted(argos_fit)
doppler_daily <- extract_predicted(doppler_fit)

# -------------------------
# 7) EXPORT
# -------------------------
readr::write_csv(argos_daily,   "Argos-SSM_daily.csv")
readr::write_csv(doppler_daily, "Doppler-SSM_daily.csv")

message("Saved:")
message("  Argos-BEFORE_SSM.csv")
message("  Doppler-BEFORE_SSM.csv")
message("  Argos-SSM_daily.csv (24-hour step)")
message("  Doppler-SSM_daily.csv (24-hour step)")
message("  SSM_excluded_tracks.csv")

# -------------------------
# 8) SANITY CHECKS
# -------------------------
cat("\nARGOS tracks kept for SSM:\n")
argos_filter_tbl %>%
  dplyr::count(keep_track) %>%
  print(n = Inf)

cat("\nDOPPLER tracks kept for SSM:\n")
doppler_filter_tbl %>%
  dplyr::count(keep_track) %>%
  print(n = Inf)

cat("\nARGOS predicted output columns:\n")
print(names(argos_daily))

cat("\nDOPPLER predicted output columns:\n")
print(names(doppler_daily))

# Optional span check
if (all(c("track_id", "date") %in% names(argos_daily))) {
  cat("\nARGOS daily predicted row counts:\n")
  argos_daily %>%
    dplyr::mutate(date = parse_time_utc(date)) %>%
    dplyr::group_by(track_id) %>%
    dplyr::summarise(
      min_dt = min(date, na.rm = TRUE),
      max_dt = max(date, na.rm = TRUE),
      n_rows = dplyr::n(),
      span_days = as.integer(as.Date(max_dt) - as.Date(min_dt)) + 1,
      .groups = "drop"
    ) %>%
    print(n = Inf)
}

if (all(c("track_id", "date") %in% names(doppler_daily))) {
  cat("\nDOPPLER daily predicted row counts:\n")
  doppler_daily %>%
    dplyr::mutate(date = parse_time_utc(date)) %>%
    dplyr::group_by(track_id) %>%
    dplyr::summarise(
      min_dt = min(date, na.rm = TRUE),
      max_dt = max(date, na.rm = TRUE),
      n_rows = dplyr::n(),
      span_days = as.integer(as.Date(max_dt) - as.Date(min_dt)) + 1,
      .groups = "drop"
    ) %>%
    print(n = Inf)
}

# ============================================================================
# 03_INDIVIDUAL_KDE
# ============================================================================
# ==========================================================
# KDE POLYGONS + AREAS FOR 4 METHODS
#   Argos-Raw   = Argos-CLEAN_postdeploy_water.xlsx
#   Argos-SSM   = argos_daily
#   Kineis-Raw  = Doppler-CLEAN_postdeploy_water.xlsx
#   Kineis-SSM  = doppler_daily
#
# Outputs:
#   - KDE_polygons_50_90_all_methods.gpkg
#   - KDE_areas_50_90_all_methods.csv
#   - KDE_area_table_by_shark_wide.csv
#   - Fig_KDE_boxplot_grouped_colors.png/.pdf
#   - Fig_KDE_boxplot_grouped_colors_LOG.png/.pdf
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(sp)
  library(adehabitatHR)
  library(purrr)
  library(readr)
  library(readxl)
  library(tidyr)
  library(ggplot2)
})

# -------------------------
# USER SETTINGS
# -------------------------
crs_proj   <- 3086
grid       <- 500
min_pts    <- 10
min_unique <- 5
levels_ud  <- c(50, 90)

argos_file   <- "Argos-CLEAN_postdeploy_water.xlsx"
doppler_file <- "Doppler-CLEAN_postdeploy_water.xlsx"

out_gpkg <- "KDE_polygons_50_90_all_methods.gpkg"
out_csv  <- "KDE_areas_50_90_all_methods.csv"

# ==========================================================
# 1) HELPERS
# ==========================================================

standardize_track_df <- function(df, label) {
  nm <- names(df)
  
  if (!("track_id" %in% nm)) {
    stop(label, ": df must have a 'track_id' column.")
  }
  
  if (!all(c("lon", "lat") %in% nm)) {
    stop(label, ": df must have 'lon' and 'lat' columns.")
  }
  
  df %>%
    dplyr::mutate(
      track_id = as.character(track_id),
      lon = as.numeric(lon),
      lat = as.numeric(lat)
    ) %>%
    dplyr::filter(
      !is.na(track_id),
      !is.na(lon),
      !is.na(lat)
    )
}

make_spdf_one_id <- function(sf_pts, tid) {
  sp_pts <- as(sf_pts, "Spatial")
  
  sp::SpatialPointsDataFrame(
    coords = sp_pts,
    data = data.frame(
      id = factor(rep(tid, length(sp_pts)))
    )
  )
}

make_kde_poly <- function(df, tid, level,
                          crs_proj = 3086,
                          grid = 500,
                          min_pts = 10,
                          min_unique = 5) {
  
  one <- df %>% dplyr::filter(track_id == tid)
  
  if (nrow(one) < min_pts) return(NULL)
  
  n_unique <- dplyr::n_distinct(paste(round(one$lon, 6), round(one$lat, 6)))
  if (n_unique < min_unique) return(NULL)
  
  pts_sf <- sf::st_as_sf(
    one,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  ) %>%
    sf::st_transform(crs_proj)
  spdf <- make_spdf_one_id(pts_sf, tid)
  
  kud <- suppressWarnings(
    tryCatch(
      adehabitatHR::kernelUD(
        spdf[, "id", drop = FALSE],
        h = "href",
        grid = grid
      ),
      error = function(e) NULL
    )
  )
  
  if (is.null(kud)) return(NULL)
  
  poly <- suppressWarnings(
    tryCatch(
      adehabitatHR::getverticeshr(kud, percent = level),
      error = function(e) NULL
    )
  )
  
  if (is.null(poly)) return(NULL)
  
  sf::st_as_sf(poly) %>%
    sf::st_set_crs(crs_proj) %>%
    sf::st_make_valid() %>%
    dplyr::mutate(
      track_id = tid,
      level = level,
      n_pts = nrow(one),
      n_unique = n_unique
    )
}

add_labels <- function(x, source, processing) {
  if (is.null(x)) return(NULL)
  x %>% dplyr::mutate(source = source, processing = processing)
}

# ==========================================================
# 2) READ RAW FILES
# ==========================================================

argos_raw_in <- readxl::read_excel(argos_file, sheet = 1)
kineis_raw_in <- readxl::read_excel(doppler_file, sheet = 1)

# standardize raw files
argos_raw   <- standardize_track_df(argos_raw_in,   "Argos raw file")
kineis_raw  <- standardize_track_df(kineis_raw_in,  "Doppler raw file")

# standardize daily SSM objects already in memory
argos_ssm   <- standardize_track_df(argos_daily,    "argos_daily (SSM)")
kineis_ssm  <- standardize_track_df(doppler_daily,  "doppler_daily (SSM)")

all_ids <- sort(unique(c(
  argos_raw$track_id,
  argos_ssm$track_id,
  kineis_raw$track_id,
  kineis_ssm$track_id
)))

# ==========================================================
# 3) OPTIONAL SIMPLE CHECK:
# which sharks pass KDE thresholds for raw Argos?
# ==========================================================

argos_raw_check <- argos_raw %>%
  dplyr::group_by(track_id) %>%
  dplyr::summarise(
    n_pts = dplyr::n(),
    n_unique = dplyr::n_distinct(paste(round(lon, 6), round(lat, 6))),
    passes_kde = n_pts >= min_pts & n_unique >= min_unique,
    .groups = "drop"
  ) %>%
  dplyr::arrange(track_id)

print(argos_raw_check)
readr::write_csv(argos_raw_check, "Argos_raw_KDE_point_check.csv")

# ==========================================================
# 4) BUILD KDE POLYGONS
# ==========================================================

kde_polys <- purrr::map(all_ids, function(tid) {
  
  polys <- list(
    purrr::map(
      levels_ud,
      ~ add_labels(
        make_kde_poly(argos_raw, tid, .x, crs_proj, grid, min_pts, min_unique),
        "Argos", "Raw"
      )
    ),
    purrr::map(
      levels_ud,
      ~ add_labels(
        make_kde_poly(argos_ssm, tid, .x, crs_proj, grid, min_pts, min_unique),
        "Argos", "SSM"
      )
    ),
    purrr::map(
      levels_ud,
      ~ add_labels(
        make_kde_poly(kineis_raw, tid, .x, crs_proj, grid, min_pts, min_unique),
        "Kineis", "Raw"
      )
    ),
    purrr::map(
      levels_ud,
      ~ add_labels(
        make_kde_poly(kineis_ssm, tid, .x, crs_proj, grid, min_pts, min_unique),
        "Kineis", "SSM"
      )
    )
  ) %>%
    purrr::flatten() %>%
    purrr::compact()
  
  if (length(polys) == 0) return(NULL)
  
  dplyr::bind_rows(polys)
}) %>%
  purrr::compact() %>%
  dplyr::bind_rows()

if (nrow(kde_polys) == 0) {
  stop("No KDE polygons produced. Try lowering min_pts/min_unique or check lon/lat.")
}

# output polygons in WGS84
kde_polys_wgs <- sf::st_transform(kde_polys, 4326)

# ==========================================================
# 5) SAVE POLYGONS
# ==========================================================

if (file.exists(out_gpkg)) file.remove(out_gpkg)

sf::st_write(
  kde_polys_wgs,
  out_gpkg,
  layer = "kde_50_90_all_methods",
  quiet = TRUE
)

# ==========================================================
# 6) BUILD AREA TABLE
# ==========================================================

kde_areas <- kde_polys %>%
  dplyr::mutate(
    area_km2 = as.numeric(sf::st_area(geometry)) / 1e6
  ) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(track_id, source, processing, level, area_km2, n_pts, n_unique) %>%
  dplyr::arrange(track_id, source, processing, level)

readr::write_csv(kde_areas, out_csv)

cat("\nKDE area table summary:\n")
print(kde_areas %>% dplyr::count(source, processing, level))
print(kde_areas)

# ==========================================================
# 7) MAKE WIDE TABLE BY SHARK
# keeps all sharks, fills missing KDEs with NA
# ==========================================================

all_combos <- tidyr::expand_grid(
  track_id = sort(unique(all_ids)),
  source = c("Argos", "Kineis"),
  processing = c("Raw", "SSM"),
  level = c(50, 90)
)

kde_table_wide <- all_combos %>%
  dplyr::left_join(
    kde_areas,
    by = c("track_id", "source", "processing", "level")
  ) %>%
  dplyr::mutate(
    method = paste(source, processing, paste0(level, "%"), sep = "_")
  ) %>%
  dplyr::select(track_id, method, area_km2) %>%
  tidyr::pivot_wider(
    names_from = method,
    values_from = area_km2
  ) %>%
  dplyr::arrange(track_id)

print(kde_table_wide)
readr::write_csv(kde_table_wide, "KDE_area_table_by_shark_wide.csv")

# ==========================================================
# 8) PREP DATA FOR PLOTTING
# ==========================================================

plot_df <- kde_areas %>%
  dplyr::mutate(
    method = paste(source, processing, sep = "-"),
    method = factor(
      method,
      levels = c("Argos-Raw", "Argos-SSM", "Kineis-Raw", "Kineis-SSM")
    ),
    level = factor(level, levels = c(50, 90), labels = c("50%", "90%")),
    method_level = paste(method, level, sep = "_")
  )

fill_vals <- c(
  "Argos-Raw_50%"   = "#9ecae1",
  "Argos-Raw_90%"   = "#3182bd",
  "Argos-SSM_50%"   = "#a1d99b",
  "Argos-SSM_90%"   = "#31a354",
  "Kineis-Raw_50%"  = "#fdae6b",
  "Kineis-Raw_90%"  = "#e6550d",
  "Kineis-SSM_50%"  = "#bcbddc",
  "Kineis-SSM_90%"  = "#756bb1"
)

# ==========================================================
# 9) RAW AREA FIGURE
# ==========================================================

p_kde <- ggplot(plot_df, aes(x = method, y = area_km2, fill = method_level)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  geom_point(
    aes(group = level),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    alpha = 0.7,
    size = 2
  ) +
  scale_fill_manual(
    values = fill_vals,
    breaks = names(fill_vals),
    labels = c(
      "Argos-Raw 50%", "Argos-Raw 90%",
      "Argos-SSM 50%", "Argos-SSM 90%",
      "Kineis-Raw 50%", "Kineis-Raw 90%",
      "Kineis-SSM 50%", "Kineis-SSM 90%"
    ),
    name = "Method and KDE level"
  ) +
  labs(
    x = "Method",
    y = expression("KDE area (km"^2*")"),
    title = "KDE area by method",
    subtitle = "Within each method, 50% and 90% KDEs are shown side by side"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

print(p_kde)

ggplot2::ggsave(
  "Fig_KDE_boxplot_grouped_colors.png",
  p_kde,
  width = 10,
  height = 5,
  dpi = 400
)

ggplot2::ggsave(
  "Fig_KDE_boxplot_grouped_colors.pdf",
  p_kde,
  width = 10,
  height = 5
)

# ==========================================================
# 10) LOG AREA FIGURE
# ==========================================================

plot_df_log <- plot_df %>%
  dplyr::mutate(log_area = log(area_km2))

p_kde_log <- ggplot(plot_df_log, aes(x = method, y = log_area, fill = method_level)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  geom_point(
    aes(group = level),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    alpha = 0.7,
    size = 2
  ) +
  scale_fill_manual(
    values = fill_vals,
    breaks = names(fill_vals),
    labels = c(
      "Argos-Raw 50%", "Argos-Raw 90%",
      "Argos-SSM 50%", "Argos-SSM 90%",
      "Kineis-Raw 50%", "Kineis-Raw 90%",
      "Kineis-SSM 50%", "Kineis-SSM 90%"
    ),
    name = "Method and KDE level"
  ) +
  labs(
    x = "Method",
    y = "log(KDE area, km²)"
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )

print(p_kde_log)

ggplot2::ggsave(
  "Fig_KDE_boxplot_grouped_colors_LOG.png",
  p_kde_log,
  width = 10,
  height = 5,
  dpi = 400
)

ggplot2::ggsave(
  "Fig_KDE_boxplot_grouped_colors_LOG.pdf",
  p_kde_log,
  width = 10,
  height = 5
)

# ==========================================================
# KDE SUMMARY STATS BY SHARK × METHOD × UD LEVEL
# ==========================================================

kde_summary_by_shark <- kde_areas %>%
  dplyr::mutate(
    method = paste(source, processing, sep = "-"),
    level = paste0(level, "%")
  ) %>%
  dplyr::group_by(track_id, method, level) %>%
  dplyr::summarise(
    n_locations = first(n_pts),
    n_unique_locations = first(n_unique),
    mean_area_km2 = mean(area_km2, na.rm = TRUE),
    sd_area_km2 = sd(area_km2, na.rm = TRUE),
    median_area_km2 = median(area_km2, na.rm = TRUE),
    min_area_km2 = min(area_km2, na.rm = TRUE),
    max_area_km2 = max(area_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(track_id, method, level)

print(kde_summary_by_shark)

readr::write_csv(
  kde_summary_by_shark,
  "KDE_summary_stats_by_shark_method_level.csv"
)

kde_table_clean <- kde_areas %>%
  dplyr::mutate(
    method = paste(source, processing, sep = "-"),
    level = paste0(level, "%")
  ) %>%
  dplyr::select(
    track_id,
    method,
    level,
    area_km2,
    n_pts,
    n_unique
  ) %>%
  dplyr::arrange(track_id, method, level)

print(kde_table_clean)

readr::write_csv(
  kde_table_clean,
  "KDE_area_by_shark_method_level.csv"
)

# ==========================================================
# ONE ROW PER SHARK: KDE AREA TABLE
# ==========================================================

kde_table_by_shark <- kde_areas %>%
  dplyr::mutate(
    method_level = paste(source, processing, paste0("UD", level), sep = "_")
  ) %>%
  dplyr::select(track_id, method_level, area_km2) %>%
  tidyr::pivot_wider(
    names_from = method_level,
    values_from = area_km2
  ) %>%
  dplyr::arrange(track_id)

print(kde_table_by_shark)

readr::write_csv(
  kde_table_by_shark,
  "KDE_area_table_one_row_per_shark.csv"
)

# ==========================================================
# FINAL TABLE: ONE ROW PER SHARK
# ==========================================================

all_combos <- tidyr::expand_grid(
  track_id = sort(unique(all_ids)),
  source = c("Argos", "Kineis"),
  processing = c("Raw", "SSM"),
  level = c(50, 90)
)

kde_table_by_shark <- all_combos %>%
  dplyr::left_join(
    kde_areas,
    by = c("track_id", "source", "processing", "level")
  ) %>%
  dplyr::mutate(
    method_level = paste(source, processing, paste0("UD", level), sep = "_")
  ) %>%
  dplyr::select(track_id, method_level, area_km2) %>%
  tidyr::pivot_wider(
    names_from = method_level,
    values_from = area_km2
  ) %>%
  dplyr::arrange(track_id)

print(kde_table_by_shark)

readr::write_csv(
  kde_table_by_shark,
  "KDE_area_table_one_row_per_shark.csv"
)

# ============================================================================
# 04_INDIVIDUAL_KDE_STATISTICS
# ============================================================================
# ==========================================================
# FULL HOME-RANGE STATISTICS + LETTERED FIGURE + SUPPLEMENTALS
# Uses existing object:
#   kde_areas
#
# Expected columns in kde_areas:
#   track_id, source, processing, level, area_km2, n_pts, n_unique
#
# Model:
#   log(area_km2) ~ workflow * ud_level + (1 | shark)
#
# Outputs:
#   KDE_LMM_TypeIII.csv
#   KDE_EMMEANS_table.csv
#   KDE_EMMEANS_workflow_within_level.csv
#   KDE_CLD_letters.csv
#   Fig_KDE_boxplot_letters_raw.png/.pdf
#   Fig_KDE_boxplot_letters_log.png/.pdf
#   Fig_KDE_LMM_diagnostics.png/.pdf
#   Supplementary_Table1_sample_sizes.csv
#   Supplementary_Table2_LMM_ANOVA.csv
#   Supplementary_Table3_Tukey_contrasts.csv
#   Supplementary_Table4_emmeans.csv
#   Supplementary_Table5_KDE_summary_stats.csv
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(multcomp)
  library(multcompView)
  library(readr)
  library(tibble)
})

# -------------------------
# 0) OPTIONAL PACKAGE INSTALLS
# -------------------------
pkgs_needed <- c("lme4", "lmerTest", "emmeans", "multcomp", "multcompView")
for (p in pkgs_needed) {
}

# -------------------------
# 1) PREP DATA
# -------------------------
plot_df <- kde_areas %>%
  dplyr::mutate(
    shark = factor(track_id),
    workflow = paste(source, processing, sep = "-"),
    workflow = factor(
      workflow,
      levels = c("Argos-Raw", "Argos-SSM", "Kineis-Raw", "Kineis-SSM")
    ),
    ud_level = factor(level, levels = c(50, 90), labels = c("50%", "90%")),
    method_level = paste(workflow, ud_level, sep = "_"),
    log_area = log(area_km2)
  ) %>%
  dplyr::filter(
    !is.na(shark),
    !is.na(workflow),
    !is.na(ud_level),
    !is.na(area_km2),
    area_km2 > 0,
    !is.na(log_area)
  )

cat("\nSample sizes by workflow and UD level:\n")
print(plot_df %>% dplyr::count(workflow, ud_level))

# -------------------------
# 2) SET CONTRASTS FOR TYPE III TESTS
# -------------------------
old_contrasts <- options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))

# -------------------------
# 3) FIT LMM
# -------------------------
m_kde <- lmerTest::lmer(
  log_area ~ workflow * ud_level + (1 | shark),
  data = plot_df,
  REML = TRUE
)

cat("\nMODEL SUMMARY:\n")
print(summary(m_kde))

cat("\nTYPE III ANOVA:\n")
anova_kde <- stats::anova(m_kde, type = 3)
print(anova_kde)

anova_out <- as.data.frame(anova_kde) %>%
  tibble::rownames_to_column(var = "term")
readr::write_csv(anova_out, "KDE_LMM_TypeIII.csv")

# -------------------------
# 4) POST-HOC COMPARISONS
#    workflow contrasts within each UD level
# -------------------------
emm_workflow_by_level <- emmeans::emmeans(
  m_kde,
  ~ workflow | ud_level
)

cat("\nEMMEANS (workflow within UD level):\n")
print(emm_workflow_by_level)

pairs_workflow_by_level <- pairs(
  emm_workflow_by_level,
  adjust = "tukey"
)

cat("\nPAIRWISE CONTRASTS (Tukey-adjusted):\n")
print(pairs_workflow_by_level)

pairs_out <- as.data.frame(pairs_workflow_by_level)
readr::write_csv(pairs_out, "KDE_EMMEANS_workflow_within_level.csv")

# -------------------------
# 5) COMPACT LETTER DISPLAYS
# -------------------------
cld_tbl <- multcomp::cld(
  emm_workflow_by_level,
  adjust = "tukey",
  Letters = letters,
  alpha = 0.05
) %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    .group = gsub(" ", "", .group)
  )

cat("\nCOMPACT LETTER DISPLAY:\n")
print(cld_tbl)

readr::write_csv(cld_tbl, "KDE_CLD_letters.csv")

emm_out <- as.data.frame(emm_workflow_by_level)
readr::write_csv(emm_out, "KDE_EMMEANS_table.csv")

# -------------------------
# 6) MODEL DIAGNOSTICS
# -------------------------
png("Fig_KDE_LMM_diagnostics.png", width = 1400, height = 700, res = 200)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

plot(
  fitted(m_kde), resid(m_kde),
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "A. Residuals vs Fitted",
  pch = 16
)
abline(h = 0, lty = 2)
lines(stats::lowess(fitted(m_kde), resid(m_kde)), col = "blue", lwd = 2)

qqnorm(resid(m_kde), main = "B. Normal Q-Q Plot", pch = 16)
qqline(resid(m_kde), lty = 2, lwd = 2)

dev.off()

pdf("Fig_KDE_LMM_diagnostics.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

plot(
  fitted(m_kde), resid(m_kde),
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "A. Residuals vs Fitted",
  pch = 16
)
abline(h = 0, lty = 2)
lines(stats::lowess(fitted(m_kde), resid(m_kde)), col = "blue", lwd = 2)

qqnorm(resid(m_kde), main = "B. Normal Q-Q Plot", pch = 16)
qqline(resid(m_kde), lty = 2, lwd = 2)

dev.off()

# -------------------------
# 7) COLORS FOR FIGURES
#    50% lighter, 90% darker within workflow
# -------------------------
fill_vals <- c(
  "Argos-Raw_50%"   = "#9ecae1",
  "Argos-Raw_90%"   = "#3182bd",
  "Argos-SSM_50%"   = "#a1d99b",
  "Argos-SSM_90%"   = "#31a354",
  "Kineis-Raw_50%"  = "#fdae6b",
  "Kineis-Raw_90%"  = "#e6550d",
  "Kineis-SSM_50%"  = "#bcbddc",
  "Kineis-SSM_90%"  = "#756bb1"
)

# -------------------------
# 8) LETTER POSITIONS
# -------------------------
letter_pos_raw <- plot_df %>%
  dplyr::group_by(workflow, ud_level) %>%
  dplyr::summarise(
    y = max(area_km2, na.rm = TRUE) * 1.08,
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    cld_tbl %>% dplyr::select(workflow, ud_level, .group),
    by = c("workflow", "ud_level")
  )

letter_pos_log <- plot_df %>%
  dplyr::group_by(workflow, ud_level) %>%
  dplyr::summarise(
    y = max(log_area, na.rm = TRUE) + 0.08 * diff(range(log_area, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    cld_tbl %>% dplyr::select(workflow, ud_level, .group),
    by = c("workflow", "ud_level")
  )

# -------------------------
# 9) RAW-SCALE FIGURE WITH LETTERS
# -------------------------
p_raw <- ggplot(
  plot_df,
  aes(x = workflow, y = area_km2, fill = method_level)
) +
  geom_boxplot(
    aes(group = interaction(workflow, ud_level)),
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  geom_point(
    aes(group = ud_level),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    alpha = 0.7,
    size = 2
  ) +
  geom_text(
    data = letter_pos_raw,
    aes(x = workflow, y = y, label = .group, group = ud_level),
    position = position_dodge(width = 0.75),
    vjust = 0,
    size = 5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = fill_vals,
    breaks = names(fill_vals),
    labels = c(
      "Argos-Raw 50%", "Argos-Raw 90%",
      "Argos-SSM 50%", "Argos-SSM 90%",
      "Kineis-Raw 50%", "Kineis-Raw 90%",
      "Kineis-SSM 50%", "Kineis-SSM 90%"
    ),
    name = "Workflow and UD level"
  ) +
  labs(
    x = "Workflow",
    y = expression("KDE area (km"^2*")"),
    title = "KDE area by tracking workflow",
    subtitle = "Letters indicate Tukey-adjusted pairwise groupings within each UD level"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  ) +
  expand_limits(y = max(letter_pos_raw$y, na.rm = TRUE) * 1.08)

print(p_raw)

ggplot2::ggsave(
  "Fig_KDE_boxplot_letters_raw.png",
  p_raw,
  width = 10,
  height = 5,
  dpi = 400
)
ggplot2::ggsave(
  "Fig_KDE_boxplot_letters_raw.pdf",
  p_raw,
  width = 10,
  height = 5
)

# -------------------------
# 10) LOG-SCALE FIGURE WITH LETTERS
# -------------------------
p_log <- ggplot(
  plot_df,
  aes(x = workflow, y = log_area, fill = method_level)
) +
  geom_boxplot(
    aes(group = interaction(workflow, ud_level)),
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  geom_point(
    aes(group = ud_level),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    alpha = 0.7,
    size = 2
  ) +
  geom_text(
    data = letter_pos_log,
    aes(x = workflow, y = y, label = .group, group = ud_level),
    position = position_dodge(width = 0.75),
    vjust = 0,
    size = 5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = fill_vals,
    breaks = names(fill_vals),
    labels = c(
      "Argos-Raw 50%", "Argos-Raw 90%",
      "Argos-SSM 50%", "Argos-SSM 90%",
      "Kineis-Raw 50%", "Kineis-Raw 90%",
      "Kineis-SSM 50%", "Kineis-SSM 90%"
    ),
    name = "Workflow and UD level"
  ) +
  labs(
    x = "Workflow",
    y = "log(KDE area, km²)"
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  ) +
  expand_limits(y = max(letter_pos_log$y, na.rm = TRUE) + 0.05)

print(p_log)

ggplot2::ggsave(
  "Fig_KDE_boxplot_letters_log.png",
  p_log,
  width = 10,
  height = 5,
  dpi = 400
)
ggplot2::ggsave(
  "Fig_KDE_boxplot_letters_log.pdf",
  p_log,
  width = 10,
  height = 5
)

# -------------------------
# 11) SUPPLEMENTARY TABLES
# -------------------------
supp_table1 <- plot_df %>%
  dplyr::group_by(workflow, ud_level) %>%
  dplyr::summarise(
    sharks = dplyr::n_distinct(shark),
    observations = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(ud_level, workflow)

supp_table2 <- as.data.frame(anova_kde) %>%
  tibble::rownames_to_column("term")

supp_table3 <- as.data.frame(pairs_workflow_by_level)

supp_table4 <- as.data.frame(emm_workflow_by_level)

supp_table5 <- plot_df %>%
  dplyr::group_by(workflow, ud_level) %>%
  dplyr::summarise(
    mean_km2 = mean(area_km2, na.rm = TRUE),
    sd_km2 = sd(area_km2, na.rm = TRUE),
    median_km2 = median(area_km2, na.rm = TRUE),
    min_km2 = min(area_km2, na.rm = TRUE),
    max_km2 = max(area_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(ud_level, workflow)

readr::write_csv(supp_table1, "Supplementary_Table1_sample_sizes.csv")
readr::write_csv(supp_table2, "Supplementary_Table2_LMM_ANOVA.csv")
readr::write_csv(supp_table3, "Supplementary_Table3_Tukey_contrasts.csv")
readr::write_csv(supp_table4, "Supplementary_Table4_emmeans.csv")
readr::write_csv(supp_table5, "Supplementary_Table5_KDE_summary_stats.csv")

# -------------------------
# 12) WORD-FRIENDLY CONSOLE TABLES
# -------------------------
cat("\n==================================================\n")
cat("Supplementary Table 1. Sample sizes by workflow and UD level\n")
cat("==================================================\n")
print(supp_table1, n = Inf)

cat("\n==================================================\n")
cat("Supplementary Table 2. Type III ANOVA for LMM\n")
cat("==================================================\n")
print(supp_table2, n = Inf)

cat("\n==================================================\n")
cat("Supplementary Table 3. Tukey-adjusted pairwise contrasts within UD level\n")
cat("==================================================\n")
print(supp_table3, n = Inf)

cat("\n==================================================\n")
cat("Supplementary Table 4. Estimated marginal means\n")
cat("==================================================\n")
print(supp_table4, n = Inf)

cat("\n==================================================\n")
cat("Supplementary Table 5. Observed KDE area summary statistics\n")
cat("==================================================\n")
print(supp_table5, n = Inf)

# -------------------------
# 13) RESTORE CONTRAST OPTIONS
# -------------------------
options(contrasts = old_contrasts$contrasts)

cat("\nAll analyses complete.\n")

# ============================================================================
# 05_INDIVIDUAL_JACCARD
# ============================================================================
# ==========================================================
# JACCARD OVERLAP AMONG KDE POLYGONS
# Argos/Kineis × Raw/SSM
# ==========================================================

# -------------------------
# 0) PACKAGES
# -------------------------
library(sf)
library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(effectsize)

# -------------------------
# 1) SETTINGS
# -------------------------
out_overlap_csv   <- "Jaccard_overlap_by_track_level.csv"
out_anova_csv     <- "Jaccard_LMM_TypeIII.csv"
out_emm_comp_csv  <- "Jaccard_EMMEANS_comparison_within_level.csv"
out_emm_level_csv <- "Jaccard_EMMEANS_level_within_comparison.csv"
out_eta_csv       <- "Jaccard_effect_sizes_partial_eta2.csv"

eps <- 0.001

# -------------------------
# 2) CHECK INPUT
# -------------------------
if (!exists("kde_polys")) {
  stop("Object 'kde_polys' not found.")
}

if (!inherits(kde_polys, "sf")) {
  stop("kde_polys must be an sf object in projected coordinates.")
}

required_cols <- c("track_id", "source", "processing", "level", "geometry")
missing_cols <- setdiff(required_cols, names(kde_polys))
if (length(missing_cols) > 0) {
  stop("kde_polys is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

# -------------------------
# 3) PREP KDE POLYGONS
# -------------------------
kde_union <- kde_polys %>%
  mutate(
    track_id   = as.character(track_id),
    source     = as.character(source),
    processing = as.character(processing),
    workflow   = paste(source, processing, sep = "-"),
    level      = as.numeric(level)
  ) %>%
  filter(
    workflow %in% c("Argos-Raw", "Argos-SSM", "Kineis-Raw", "Kineis-SSM"),
    level %in% c(50, 90)
  ) %>%
  group_by(track_id, workflow, level) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_as_sf() %>%
  st_make_valid()

print(kde_union %>% st_drop_geometry() %>% count(workflow, level))

if (nrow(kde_union) == 0) {
  stop("No polygons left after filtering. Check source/processing/level values.")
}

# -------------------------
# 4) HELPER FUNCTION
# -------------------------
calc_jaccard_one <- function(poly1, poly2) {
  poly1 <- st_make_valid(poly1)
  poly2 <- st_make_valid(poly2)
  
  inter <- suppressWarnings(st_intersection(poly1, poly2))
  union <- suppressWarnings(st_union(poly1, poly2))
  
  inter_area <- if (nrow(inter) == 0) 0 else as.numeric(st_area(inter))
  union_area <- as.numeric(st_area(union))
  
  if (is.na(union_area) || union_area == 0) return(NA_real_)
  
  inter_area / union_area
}

# -------------------------
# 5) DEFINE COMPARISONS
# -------------------------
comparison_tbl <- tibble::tribble(
  ~comparison,                 ~workflow1,    ~workflow2,
  "Argos Raw vs Argos SSM",    "Argos-Raw",   "Argos-SSM",
  "Kineis Raw vs Kineis SSM",  "Kineis-Raw",  "Kineis-SSM",
  "Argos Raw vs Kineis Raw",   "Argos-Raw",   "Kineis-Raw",
  "Argos SSM vs Kineis SSM",   "Argos-SSM",   "Kineis-SSM"
)

# -------------------------
# 6) CALCULATE JACCARD
# -------------------------
all_tracks <- sort(unique(kde_union$track_id))
all_levels <- c(50, 90)

jaccard_df <- map_dfr(all_tracks, function(tid) {
  map_dfr(all_levels, function(lev) {
    
    one_set <- kde_union %>%
      filter(track_id == tid, level == lev)
    
    pmap_dfr(comparison_tbl, function(comparison, workflow1, workflow2) {
      
      poly1 <- one_set %>% filter(workflow == workflow1)
      poly2 <- one_set %>% filter(workflow == workflow2)
      
      if (nrow(poly1) == 0 || nrow(poly2) == 0) {
        return(tibble(
          track_id   = tid,
          level      = lev,
          comparison = comparison,
          workflow1  = workflow1,
          workflow2  = workflow2,
          jaccard    = NA_real_
        ))
      }
      
      tibble(
        track_id   = tid,
        level      = lev,
        comparison = comparison,
        workflow1  = workflow1,
        workflow2  = workflow2,
        jaccard    = calc_jaccard_one(poly1, poly2)
      )
    })
  })
}) %>%
  filter(!is.na(jaccard)) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c(
        "Argos Raw vs Argos SSM",
        "Kineis Raw vs Kineis SSM",
        "Argos Raw vs Kineis Raw",
        "Argos SSM vs Kineis SSM"
      )
    ),
    level = factor(level, levels = c(50, 90), labels = c("50%", "90%"))
  )

if (nrow(jaccard_df) == 0) {
  stop("jaccard_df is empty. Check whether all required workflow pairs exist for each track and level.")
}

write_csv(jaccard_df, out_overlap_csv)

cat("\nJaccard overlap summary:\n")
print(jaccard_df %>% count(comparison, level))
print(summary(jaccard_df$jaccard))

# -------------------------
# 7) FIT MODEL
# -------------------------
m_jaccard <- lmer(
  jaccard ~ comparison * level + (1 | track_id),
  data = jaccard_df,
  REML = TRUE
)

cat("\nModel summary:\n")
print(summary(m_jaccard))

# -------------------------
# 8) TYPE III ANOVA
# -------------------------
anova_jaccard <- anova(m_jaccard, type = 3, ddf = "Satterthwaite")
print(anova_jaccard)

anova_jaccard_df <- as.data.frame(anova_jaccard) %>%
  rownames_to_column("term")

write_csv(anova_jaccard_df, out_anova_csv)

# -------------------------
# 9) EFFECT SIZES
# -------------------------
eta_jaccard <- effectsize::eta_squared(
  anova_jaccard,
  partial = TRUE
)

print(eta_jaccard)

eta_jaccard_df <- as.data.frame(eta_jaccard)
write_csv(eta_jaccard_df, out_eta_csv)

# -------------------------
# 10) FIGURE
# -------------------------
fill_vals <- c(
  "50%" = "#9ecae1",
  "90%" = "#3182bd"
)

p_jaccard <- ggplot(
  jaccard_df,
  aes(x = comparison, y = jaccard, fill = level)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    alpha = 0.7,
    size = 2
  ) +
  scale_fill_manual(values = fill_vals, name = "UD level") +
  labs(
    x = "Tracking comparison",
    y = "Jaccard similarity",
    title = "Spatial overlap among KDE polygons",
    subtitle = "Overlap shown for 50% and 90% utilization distributions"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

print(p_jaccard)

ggsave("Fig_Jaccard_boxplot.png", p_jaccard, width = 10, height = 5, dpi = 400)
ggsave("Fig_Jaccard_boxplot.pdf", p_jaccard, width = 10, height = 5)

# -------------------------
# 11) MODEL DIAGNOSTICS
# -------------------------    
png(
  "Supplementary_Figure_Jaccard_model_diagnostics.png",
  width = 1800, height = 900, res = 300
)

par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

plot(
  fitted(m_jaccard),
  resid(m_jaccard),
  pch = 16,
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "A. Residuals vs Fitted"
)
abline(h = 0, lty = 2)
lines(lowess(fitted(m_jaccard), resid(m_jaccard)), col = "blue", lwd = 2)

qqnorm(resid(m_jaccard), pch = 16, main = "B. Normal Q-Q Plot")
qqline(resid(m_jaccard), lwd = 2, lty = 2)

dev.off()

pdf(
  "Supplementary_Figure_Jaccard_model_diagnostics.pdf",
  width = 10, height = 5
)

par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

plot(
  fitted(m_jaccard),
  resid(m_jaccard),
  pch = 16,
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "A. Residuals vs Fitted"
)
abline(h = 0, lty = 2)
lines(lowess(fitted(m_jaccard), resid(m_jaccard)), col = "blue", lwd = 2)

qqnorm(resid(m_jaccard), pch = 16, main = "B. Normal Q-Q Plot")
qqline(resid(m_jaccard), lwd = 2, lty = 2)

dev.off()

cat("\nDone. Output files written to working directory.\n") 

# ============================================================================
# 06_TRANSMISSION_SPACING_AND_YIELD
# ============================================================================
# ============================================================
# Argos vs Kineis transmission spacing and transmission yield
# Cleaned script with explicit namespaces
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(MASS)
  library(patchwork)
})

# ------------------------------------------------------------
# 1) Files
# ------------------------------------------------------------
argos_file  <- "Argos-CLEAN_postdeploy_water.csv"
kineis_file <- "Doppler-CLEAN_postdeploy_water.csv"

if (!file.exists(argos_file)) stop("Missing file: ", argos_file)
if (!file.exists(kineis_file)) stop("Missing file: ", kineis_file)

# ------------------------------------------------------------
# 2) Read data
# ------------------------------------------------------------
argos <- readr::read_csv(argos_file, show_col_types = FALSE) %>%
  dplyr::mutate(
    system = "Argos",
    timestamp = lubridate::ymd_hms(timestamp, tz = "UTC")
  )

kineis <- readr::read_csv(kineis_file, show_col_types = FALSE) %>%
  dplyr::mutate(
    system = "Kineis",
    timestamp = lubridate::ymd_hms(timestamp, tz = "UTC")
  )

# Keep only columns that exist
keep_cols <- c(
  "system", "track_id", "nickname", "species_code",
  "deployment_date", "timestamp", "lat", "lon"
)

argos_keep  <- intersect(names(argos), keep_cols)
kineis_keep <- intersect(names(kineis), keep_cols)

dat <- dplyr::bind_rows(
  argos[, argos_keep, drop = FALSE],
  kineis[, kineis_keep, drop = FALSE]
) %>%
  dplyr::filter(!is.na(track_id), !is.na(timestamp)) %>%
  dplyr::mutate(
    track_id = factor(track_id),
    system   = factor(system, levels = c("Argos", "Kineis"))
  )

# ------------------------------------------------------------
# 3) Remove duplicate timestamps within shark x system
# ------------------------------------------------------------
dat_unique <- dat %>%
  dplyr::arrange(track_id, system, timestamp) %>%
  dplyr::group_by(track_id, system) %>%
  dplyr::distinct(timestamp, .keep_all = TRUE) %>%
  dplyr::ungroup()

# ------------------------------------------------------------
# 4) Calculate time intervals between transmissions
# ------------------------------------------------------------
dt_df <- dat_unique %>%
  dplyr::arrange(track_id, system, timestamp) %>%
  dplyr::group_by(track_id, system) %>%
  dplyr::mutate(
    dt_hours = as.numeric(difftime(timestamp, dplyr::lag(timestamp), units = "hours"))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(dt_hours), dt_hours > 0)

# Per-shark summary
dt_animal <- dt_df %>%
  dplyr::group_by(track_id, system) %>%
  dplyr::summarise(
    n_intervals     = dplyr::n(),
    median_dt_hours = median(dt_hours, na.rm = TRUE),
    mean_dt_hours   = mean(dt_hours, na.rm = TRUE),
    sd_dt_hours     = sd(dt_hours, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    system = factor(system, levels = c("Argos", "Kineis"))
  )

# Wide version for paired Wilcoxon
dt_wide <- dt_animal %>%
  dplyr::select(track_id, system, median_dt_hours) %>%
  tidyr::pivot_wider(names_from = system, values_from = median_dt_hours)

paired <- dt_wide %>%
  dplyr::filter(!is.na(Argos), !is.na(Kineis))

# ------------------------------------------------------------
# 5) Statistics for interval comparisons
# ------------------------------------------------------------
wilcox_out <- wilcox.test(
  paired$Argos,
  paired$Kineis,
  paired = TRUE,
  exact = FALSE
)

dt_model <- dt_animal %>%
  dplyr::filter(median_dt_hours > 0) %>%
  dplyr::mutate(
    log_median_dt = log(median_dt_hours)
  )

m_dt <- lmerTest::lmer(
  log_median_dt ~ system + (1 | track_id),
  data = dt_model,
  REML = TRUE
)

aov_dt <- anova(m_dt)
p_lmm <- aov_dt["system", "Pr(>F)"]
p_wilcox <- wilcox_out$p.value

# ------------------------------------------------------------
# 6) Shared days at liberty for panel C
# ------------------------------------------------------------
# This makes x-values the same for Argos and Kineis within shark
shark_days <- dat_unique %>%
  dplyr::group_by(track_id) %>%
  dplyr::summarise(
    first_time_all = min(timestamp, na.rm = TRUE),
    last_time_all  = max(timestamp, na.rm = TRUE),
    days_at_liberty = as.numeric(difftime(last_time_all, first_time_all, units = "days")),
    .groups = "drop"
  ) %>%
  dplyr::filter(days_at_liberty > 0)

# Number of unique transmissions per shark x system
tx_summary <- dat_unique %>%
  dplyr::group_by(track_id, system) %>%
  dplyr::summarise(
    n_transmissions = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::left_join(shark_days, by = "track_id") %>%
  dplyr::mutate(
    transmissions_per_day = n_transmissions / days_at_liberty,
    system = factor(system, levels = c("Argos", "Kineis"))
  )

# ------------------------------------------------------------
# 7) Count model for transmission totals
# ------------------------------------------------------------
m_pois <- lme4::glmer(
  n_transmissions ~ system + offset(log(days_at_liberty)) + (1 | track_id),
  family = poisson(link = "log"),
  data = tx_summary
)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp  <- residuals(model, type = "pearson")
  chisq <- sum(rp^2)
  ratio <- chisq / rdf
  c(chisq = chisq, ratio = ratio, rdf = rdf)
}

od <- overdisp_fun(m_pois)

if (is.finite(od["ratio"]) && od["ratio"] > 1.5) {
  final_model <- lme4::glmer.nb(
    n_transmissions ~ system + offset(log(days_at_liberty)) + (1 | track_id),
    data = tx_summary
  )
  model_name <- "NB-GLMM"
} else {
  final_model <- m_pois
  model_name <- "Poisson GLMM"
}

coef_tab <- summary(final_model)$coefficients
coef_name <- grep("^system", rownames(coef_tab), value = TRUE)

p_glmm <- coef_tab[coef_name, "Pr(>|z|)"]
beta_tx <- lme4::fixef(final_model)[coef_name]
rate_ratio <- exp(beta_tx)

# ------------------------------------------------------------
# 8) Prediction data for panel C
# ------------------------------------------------------------
pred_grid <- expand.grid(
  days_at_liberty = seq(
    min(tx_summary$days_at_liberty, na.rm = TRUE),
    max(tx_summary$days_at_liberty, na.rm = TRUE),
    length.out = 200
  ),
  system = factor(c("Argos", "Kineis"), levels = c("Argos", "Kineis"))
)

X_grid <- model.matrix(~ system, data = pred_grid)
beta_hat <- lme4::fixef(final_model)
V_beta <- as.matrix(vcov(final_model))

eta_grid <- as.vector(X_grid %*% beta_hat) + log(pred_grid$days_at_liberty)
se_grid <- sqrt(diag(X_grid %*% V_beta %*% t(X_grid)))

pred_grid$fit   <- exp(eta_grid)
pred_grid$lower <- exp(eta_grid - 1.96 * se_grid)
pred_grid$upper <- exp(eta_grid + 1.96 * se_grid)

# ------------------------------------------------------------
# 9) Formatting helpers
# ------------------------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, format = "f", digits = 3))
}

lab_A <- paste0("LMM: ", fmt_p(p_lmm))
lab_B <- paste0("Wilcoxon: ", fmt_p(p_wilcox))
lab_C <- paste0(model_name, ": ", fmt_p(p_glmm), "\nRR = ", round(rate_ratio, 2))

pal <- c("Argos" = "#D55E00", "Kineis" = "#009E73")

theme_pub <- ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 18),
    axis.text = ggplot2::element_text(size = 16, color = "black"),
    legend.position = "none",
    plot.tag = ggplot2::element_text(face = "bold", size = 16)
  )

# ------------------------------------------------------------
# 10) Panel A
# ------------------------------------------------------------
pA <- ggplot2::ggplot(dt_animal, ggplot2::aes(x = system, y = median_dt_hours)) +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = system),
    outlier.shape = NA,
    width = 0.55,
    alpha = 0.85,
    color = "black"
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = system),
    width = 0.08,
    size = 3,
    alpha = 0.9
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = max(dt_animal$median_dt_hours, na.rm = TRUE) * 1.12,
    label = lab_A,
    size = 5
  ) +
  ggplot2::scale_fill_manual(values = pal) +
  ggplot2::scale_color_manual(values = pal) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = NULL,
    y = "Median time between transmissions (hours)"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_pub

# ------------------------------------------------------------
# 11) Panel B
# ------------------------------------------------------------
pB <- ggplot2::ggplot(
  dt_animal,
  ggplot2::aes(x = system, y = median_dt_hours, group = track_id)
) +
  ggplot2::geom_line(color = "grey65", linewidth = 0.7) +
  ggplot2::geom_point(
    ggplot2::aes(fill = system),
    shape = 21,
    size = 3.2,
    color = "black"
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = max(dt_animal$median_dt_hours, na.rm = TRUE) * 1.12,
    label = lab_B,
    size = 5
  ) +
  ggplot2::scale_fill_manual(values = pal) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = NULL,
    y = "Median time between transmissions (hours)"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_pub

# ------------------------------------------------------------
# 12) Panel C
# ------------------------------------------------------------
pC <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = pred_grid,
    ggplot2::aes(
      x = days_at_liberty,
      ymin = lower,
      ymax = upper,
      fill = system
    ),
    alpha = 0.18
  ) +
  ggplot2::geom_line(
    data = pred_grid,
    ggplot2::aes(
      x = days_at_liberty,
      y = fit,
      color = system
    ),
    linewidth = 1.2
  ) +
  ggplot2::geom_point(
    data = tx_summary,
    ggplot2::aes(
      x = days_at_liberty,
      y = n_transmissions,
      color = system
    ),
    size = 3,
    alpha = 0.65
  ) +
  ggplot2::annotate(
    "text",
    x = min(tx_summary$days_at_liberty, na.rm = TRUE) + 2,
    y = max(tx_summary$n_transmissions, na.rm = TRUE) * 0.95,
    hjust = 0,
    label = lab_C,
    size = 5
  ) +
  ggplot2::scale_color_manual(values = pal) +
  ggplot2::scale_fill_manual(values = pal) +
  ggplot2::labs(
    x = "Days at liberty",
    y = "Number of transmissions"
  ) +
  theme_pub

# ------------------------------------------------------------
# 13) Combine figure
# ------------------------------------------------------------
p_multi <- (pA | pB | pC) +
  patchwork::plot_annotation(tag_levels = "A")

print(p_multi)

# Save
ggplot2::ggsave(
  "Fig_transmission_comparison_cleaned.png",
  p_multi,
  width = 12,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  "Fig_transmission_comparison_cleaned.pdf",
  p_multi,
  width = 12,
  height = 4.5,
  bg = "white"
)

# ------------------------------------------------------------
# 14) Export summaries
# ------------------------------------------------------------
readr::write_csv(dt_df, "Transmission_intervals_by_fix.csv")
readr::write_csv(dt_animal, "Transmission_intervals_by_shark.csv")
readr::write_csv(tx_summary, "Transmission_counts_by_shark_system.csv")

# ------------------------------------------------------------
# 15) Print results
# ------------------------------------------------------------
cat("\n==============================\n")
cat("Panel A: LMM on log median interval\n")
cat("==============================\n")
print(summary(m_dt))
print(anova(m_dt))

cat("\n==============================\n")
cat("Panel B: Paired Wilcoxon\n")
cat("==============================\n")
print(wilcox_out)

cat("\n==============================\n")
cat("Panel C: Count model\n")
cat("==============================\n")
print(summary(final_model))
cat("\nRate ratio (Kineis vs Argos): ", round(rate_ratio, 3), "\n")
cat("Model p-value: ", p_glmm, "\n")
cat("Overdispersion ratio from Poisson check: ", round(od["ratio"], 3), "\n")
# =============================================================================
# 07_LOCATION_CLASS_DISTRIBUTION_AND_CONTINUOUS_ERROR
# =============================================================================

lc_levels <- c("3", "2", "1", "0", "A", "B")

argos_lc <- argos_clean %>%
  dplyr::transmute(
    track_id = as.character(track_id),
    system = "Argos",
    quality_class = as.character(loc_quality)
  )

kineis_lc <- doppler_clean %>%
  dplyr::transmute(
    track_id = as.character(track_id),
    system = "Kineis",
    quality_class = as.character(location_class)
  )

shared_lc_ids <- intersect(unique(argos_lc$track_id), unique(kineis_lc$track_id))

lc_df <- dplyr::bind_rows(argos_lc, kineis_lc) %>%
  dplyr::filter(
    track_id %in% shared_lc_ids,
    quality_class %in% lc_levels
  ) %>%
  dplyr::mutate(
    system = factor(system, levels = c("Argos", "Kineis")),
    quality_class = factor(quality_class, levels = lc_levels)
  )

lc_counts <- lc_df %>%
  dplyr::count(system, quality_class, name = "n") %>%
  dplyr::group_by(system) %>%
  dplyr::mutate(proportion = n / sum(n)) %>%
  dplyr::ungroup()

lc_table <- xtabs(n ~ system + quality_class, data = lc_counts)
lc_chisq <- stats::chisq.test(lc_table)
lc_cramers_v <- sqrt(
  unname(lc_chisq$statistic) /
    (sum(lc_table) * min(nrow(lc_table) - 1, ncol(lc_table) - 1))
)

lc_residuals <- as.data.frame(as.table(lc_chisq$stdres)) %>%
  dplyr::rename(
    system = Var1,
    quality_class = Var2,
    standardized_residual = Freq
  )

safe_write_csv(lc_counts, "Supplemental_LC_counts_and_proportions.csv")
safe_write_csv(lc_residuals, "Supplemental_LC_standardized_residuals.csv")
safe_write_csv(
  tibble::tibble(
    chi_square = unname(lc_chisq$statistic),
    df = unname(lc_chisq$parameter),
    p_value = lc_chisq$p.value,
    cramers_v = lc_cramers_v
  ),
  "Supplemental_LC_chisquare_and_effect_size.csv"
)

p_lc_prop <- ggplot2::ggplot(
  lc_counts,
  ggplot2::aes(x = quality_class, y = proportion, fill = system)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.78),
    width = 0.70,
    color = "black",
    linewidth = 0.25
  ) +
  ggplot2::scale_fill_manual(values = c("Argos" = "#D55E00", "Kineis" = "#009E73")) +
  ggplot2::scale_y_continuous(labels = scales::label_percent()) +
  ggplot2::labs(
    x = "Location class",
    y = "Proportion of locations",
    fill = "System"
  ) +
  ggplot2::theme_classic(base_size = 16)

p_lc_resid <- ggplot2::ggplot(
  lc_residuals,
  ggplot2::aes(x = quality_class, y = standardized_residual, fill = system)
) +
  ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dashed") +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.78),
    width = 0.70,
    color = "black",
    linewidth = 0.25
  ) +
  ggplot2::scale_fill_manual(values = c("Argos" = "#D55E00", "Kineis" = "#009E73")) +
  ggplot2::labs(
    x = "Location class",
    y = "Standardized residual",
    fill = "System"
  ) +
  ggplot2::theme_classic(base_size = 16)

ggplot2::ggsave(
  file.path("Main_Figures", "Figure_Location_Class_Distribution.png"),
  p_lc_prop / p_lc_resid,
  width = 9, height = 9, dpi = 500, bg = "white"
)
ggplot2::ggsave(
  file.path("Main_Figures", "Figure_Location_Class_Distribution.pdf"),
  p_lc_prop / p_lc_resid,
  width = 9, height = 9, bg = "white"
)

argos_lookup <- tibble::tibble(
  quality_class = factor(lc_levels, levels = lc_levels),
  argos_expected_m = c(125, 375, 1000, 2500, 1000, 5000)
)

compare_df <- doppler_clean %>%
  dplyr::filter(track_id %in% shared_lc_ids) %>%
  dplyr::transmute(
    track_id = factor(track_id),
    quality_class = factor(as.character(location_class), levels = lc_levels),
    kineis_error_m = as.numeric(doppler_error_m)
  ) %>%
  dplyr::filter(
    !is.na(quality_class),
    !is.na(kineis_error_m),
    kineis_error_m > 0
  ) %>%
  dplyr::left_join(argos_lookup, by = "quality_class") %>%
  dplyr::mutate(
    ratio = kineis_error_m / argos_expected_m,
    log10_ratio = log10(ratio),
    raw_bias_m = argos_expected_m - kineis_error_m,
    observed_direction = dplyr::case_when(
      ratio < 1 ~ "Kineis lower than Argos expected",
      ratio > 1 ~ "Kineis greater than Argos expected",
      TRUE ~ "Equal"
    )
  )

error_summary <- compare_df %>%
  dplyr::group_by(quality_class) %>%
  dplyr::summarise(
    n = dplyr::n(),
    argos_expected_m = dplyr::first(argos_expected_m),
    mean_kineis_m = mean(kineis_error_m),
    median_kineis_m = median(kineis_error_m),
    sd_kineis_m = sd(kineis_error_m),
    mean_ratio = mean(ratio),
    median_ratio = median(ratio),
    proportion_kineis_lower = mean(ratio < 1),
    proportion_kineis_greater = mean(ratio > 1),
    .groups = "drop"
  )

m_ratio <- lmerTest::lmer(
  log10_ratio ~ 0 + quality_class + (1 | track_id),
  data = compare_df,
  REML = TRUE
)

ratio_anova_model <- lmerTest::lmer(
  log10_ratio ~ quality_class + (1 | track_id),
  data = compare_df,
  REML = TRUE
)

ratio_anova <- stats::anova(ratio_anova_model)
ratio_emm_log <- emmeans::emmeans(m_ratio, ~ quality_class)

ratio_model_df <- as.data.frame(ratio_emm_log) %>%
  dplyr::mutate(
    ratio = 10^emmean,
    lower.CL = 10^lower.CL,
    upper.CL = 10^upper.CL,
    model_status = dplyr::case_when(
      lower.CL > 1 ~ "Kineis greater than Argos expected",
      upper.CL < 1 ~ "Kineis lower than Argos expected",
      TRUE ~ "Not different from Argos expected"
    ),
    sig_label = dplyr::if_else(lower.CL > 1 | upper.CL < 1, "*", "")
  )

ratio_tests <- as.data.frame(
  emmeans::test(ratio_emm_log, null = 0, adjust = "none")
)

safe_write_csv(error_summary, "Supplemental_Error_Ratio_Observed_Summary.csv")
safe_write_csv(
  as.data.frame(ratio_anova) %>% tibble::rownames_to_column("term"),
  "Supplemental_Error_Ratio_TypeIII.csv"
)
safe_write_csv(ratio_model_df, "Supplemental_Error_Ratio_Model_Estimates.csv")
safe_write_csv(ratio_tests, "Supplemental_Error_Ratio_Tests_Against_One.csv")
save_model(m_ratio, "Error_Ratio_LMM.rds")

prop_df <- compare_df %>%
  dplyr::count(quality_class, observed_direction, name = "n") %>%
  dplyr::group_by(quality_class) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

error_cols <- c(
  "Kineis lower than Argos expected" = "#0072B2",
  "Not different from Argos expected" = "grey65",
  "Kineis greater than Argos expected" = "#D55E00",
  "Equal" = "grey65"
)

pA_error <- ggplot2::ggplot(
  ratio_model_df,
  ggplot2::aes(x = quality_class, y = ratio, color = model_status)
) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.12, linewidth = 0.8
  ) +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = sig_label, y = upper.CL * 1.13),
    color = "black", size = 6
  ) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_color_manual(values = error_cols, name = NULL) +
  ggplot2::labs(
    x = "Location class",
    y = "Model-estimated Kineis-to-Argos\nexpected-error ratio"
  ) +
  ggplot2::theme_classic(base_size = 16)

pB_error <- ggplot2::ggplot(
  compare_df,
  ggplot2::aes(x = quality_class, y = ratio, fill = observed_direction)
) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  ggplot2::geom_jitter(width = 0.12, alpha = 0.08, size = 0.6) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_fill_manual(values = error_cols, name = NULL) +
  ggplot2::labs(
    x = "Location class",
    y = "Observed Kineis-to-Argos\nexpected-error ratio"
  ) +
  ggplot2::theme_classic(base_size = 16)

pC_error <- ggplot2::ggplot(
  prop_df,
  ggplot2::aes(x = quality_class, y = prop, fill = observed_direction)
) +
  ggplot2::geom_col(color = "black", linewidth = 0.25) +
  ggplot2::scale_y_continuous(labels = scales::label_percent()) +
  ggplot2::scale_fill_manual(values = error_cols, name = NULL) +
  ggplot2::labs(
    x = "Location class",
    y = "Proportion of locations"
  ) +
  ggplot2::theme_classic(base_size = 16)

fig_error <- (pA_error / pB_error / pC_error) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 20, face = "bold"),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path("Main_Figures", "Figure_Error_Ratio_Comparison.png"),
  fig_error, width = 9, height = 14, dpi = 500, bg = "white"
)
ggplot2::ggsave(
  file.path("Main_Figures", "Figure_Error_Ratio_Comparison.pdf"),
  fig_error, width = 9, height = 14, bg = "white"
)

# =============================================================================
# 08_SPECIES_LEVEL_POOLED_KDE_AND_OVERLAP
# =============================================================================

species_lookup <- dplyr::bind_rows(
  argos_clean %>% dplyr::select(track_id, species_code),
  doppler_clean %>% dplyr::select(track_id, species_code)
) %>%
  dplyr::distinct(track_id, species_code) %>%
  dplyr::mutate(track_id = as.character(track_id))

add_species <- function(df) {
  df %>%
    dplyr::mutate(track_id = as.character(track_id)) %>%
    dplyr::left_join(species_lookup, by = "track_id")
}

species_inputs <- list(
  `Argos-Raw` = add_species(argos_raw),
  `Argos-SSM` = add_species(argos_ssm),
  `Kineis-Raw` = add_species(kineis_raw),
  `Kineis-SSM` = add_species(kineis_ssm)
)

make_species_kde <- function(df, species_value, workflow_value, level_value) {
  one <- df %>%
    dplyr::filter(species_code == species_value) %>%
    dplyr::filter(!is.na(lon), !is.na(lat))

  n_pts <- nrow(one)
  n_unique <- dplyr::n_distinct(paste(round(one$lon, 6), round(one$lat, 6)))
  if (n_pts < min_pts || n_unique < min_unique) return(NULL)

  pts <- sf::st_as_sf(one, coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    sf::st_transform(crs_proj)

  sp_pts <- as(pts, "Spatial")
  spdf <- sp::SpatialPointsDataFrame(
    coords = sp_pts,
    data = data.frame(id = factor(rep(species_value, length(sp_pts))))
  )

  kud <- tryCatch(
    adehabitatHR::kernelUD(spdf[, "id", drop = FALSE], h = "href", grid = grid),
    error = function(e) NULL
  )
  if (is.null(kud)) return(NULL)

  poly <- tryCatch(
    adehabitatHR::getverticeshr(kud, percent = level_value),
    error = function(e) NULL
  )
  if (is.null(poly)) return(NULL)

  sf::st_as_sf(poly) %>%
    sf::st_set_crs(crs_proj) %>%
    sf::st_make_valid() %>%
    dplyr::mutate(
      species_code = species_value,
      workflow = workflow_value,
      level = level_value,
      n_pts = n_pts,
      n_unique = n_unique
    )
}

species_codes <- sort(unique(stats::na.omit(species_lookup$species_code)))

species_kde_polys <- purrr::imap_dfr(species_inputs, function(df, workflow_name) {
  purrr::map_dfr(species_codes, function(sp_code) {
    purrr::map_dfr(levels_ud, function(lev) {
      make_species_kde(df, sp_code, workflow_name, lev)
    })
  })
})

if (nrow(species_kde_polys) > 0) {
  species_kde_areas <- species_kde_polys %>%
    dplyr::mutate(area_km2 = as.numeric(sf::st_area(geometry)) / 1e6) %>%
    sf::st_drop_geometry()

  safe_write_csv(species_kde_areas, "Supplemental_Species_KDE_Areas.csv")

  if (file.exists("Species_KDE_polygons.gpkg")) file.remove("Species_KDE_polygons.gpkg")
  sf::st_write(
    sf::st_transform(species_kde_polys, 4326),
    "Species_KDE_polygons.gpkg",
    layer = "species_kde",
    quiet = TRUE
  )

  species_model_df <- species_kde_areas %>%
    dplyr::mutate(
      species_code = factor(species_code),
      workflow = factor(
        workflow,
        levels = c("Argos-Raw", "Argos-SSM", "Kineis-Raw", "Kineis-SSM")
      ),
      ud_level = factor(level, levels = c(50, 90), labels = c("50%", "90%")),
      log_area = log(area_km2)
    )

  m_species_kde <- stats::lm(
    log_area ~ species_code * workflow + ud_level,
    data = species_model_df
  )

  species_drop1 <- stats::drop1(m_species_kde, test = "F") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term")

  safe_write_csv(species_drop1, "Supplemental_Species_KDE_DropOne_F_Tests.csv")
  save_model(m_species_kde, "Species_KDE_LM.rds")

  p_species_area <- ggplot2::ggplot(
    species_model_df,
    ggplot2::aes(x = workflow, y = area_km2, fill = ud_level)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      color = "black", linewidth = 0.25
    ) +
    ggplot2::scale_y_log10(labels = scales::label_number()) +
    ggplot2::facet_wrap(~ species_code, scales = "free_y") +
    ggplot2::labs(
      x = "Tracking workflow",
      y = expression("Species-level KDE area (km"^2*")"),
      fill = "UD level"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  ggplot2::ggsave(
    file.path("Main_Figures", "Figure_Species_KDE_Areas.png"),
    p_species_area, width = 11, height = 8, dpi = 500, bg = "white"
  )

  species_union <- species_kde_polys %>%
    dplyr::group_by(species_code, workflow, level) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    sf::st_make_valid()

  species_jaccard <- purrr::map_dfr(species_codes, function(sp_code) {
    purrr::map_dfr(c(50, 90), function(lev) {
      one <- species_union %>% dplyr::filter(species_code == sp_code, level == lev)
      purrr::pmap_dfr(comparison_tbl, function(comparison, workflow1, workflow2) {
        p1 <- one %>% dplyr::filter(workflow == workflow1)
        p2 <- one %>% dplyr::filter(workflow == workflow2)
        if (nrow(p1) == 0 || nrow(p2) == 0) return(tibble::tibble())
        tibble::tibble(
          species_code = sp_code,
          level = lev,
          comparison = comparison,
          jaccard = calc_jaccard_one(p1, p2)
        )
      })
    })
  })

  safe_write_csv(species_jaccard, "Supplemental_Species_Jaccard_Overlap.csv")

  if (nrow(species_jaccard) >= 8) {
    species_jaccard_model_df <- species_jaccard %>%
      dplyr::mutate(
        species_code = factor(species_code),
        comparison = factor(comparison),
        ud_level = factor(level, levels = c(50, 90), labels = c("50%", "90%"))
      )

    m_species_jaccard <- stats::lm(
      jaccard ~ species_code * comparison + ud_level,
      data = species_jaccard_model_df
    )

    species_jaccard_drop1 <- stats::drop1(m_species_jaccard, test = "F") %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term")

    safe_write_csv(
      species_jaccard_drop1,
      "Supplemental_Species_Jaccard_DropOne_F_Tests.csv"
    )
    save_model(m_species_jaccard, "Species_Jaccard_LM.rds")
  }
}

# =============================================================================
# 09_BEHAVIORAL_STATE_COMPARISON
# =============================================================================
# This section runs when two state-assignment files are available. Each file
# should contain track ID, timestamp/date, and a behavioral state. Optional
# projected x/y coordinates allow step-length and turning-angle summaries.

find_first_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}

standardize_behavior <- function(path, system_name) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  id_col <- find_first_col(names(x), c("track_id", "id", "shark", "tag_id"))
  time_col <- find_first_col(names(x), c("timestamp", "date", "datetime", "time"))
  state_col <- find_first_col(names(x), c("state", "behavior", "behaviour", "state_label"))
  x_col <- find_first_col(names(x), c("x", "x_coord", "utm_x", "easting"))
  y_col <- find_first_col(names(x), c("y", "y_coord", "utm_y", "northing"))

  if (any(is.na(c(id_col, time_col, state_col)))) {
    stop(path, " must include ID, timestamp/date, and state columns.")
  }

  out <- x %>%
    dplyr::transmute(
      track_id = as.character(.data[[id_col]]),
      timestamp = lubridate::as_datetime(.data[[time_col]], tz = "UTC"),
      state = as.character(.data[[state_col]]),
      x = if (!is.na(x_col)) as.numeric(.data[[x_col]]) else NA_real_,
      y = if (!is.na(y_col)) as.numeric(.data[[y_col]]) else NA_real_,
      system = system_name
    ) %>%
    dplyr::filter(!is.na(track_id), !is.na(timestamp), !is.na(state))

  out
}

argos_behavior_file <- "Argos_BayesMove_states.csv"
kineis_behavior_file <- "Kineis_BayesMove_states.csv"

if (
  RUN_BEHAVIOR_IF_AVAILABLE &&
  file.exists(argos_behavior_file) &&
  file.exists(kineis_behavior_file)
) {
  behavior_df <- dplyr::bind_rows(
    standardize_behavior(argos_behavior_file, "Argos"),
    standardize_behavior(kineis_behavior_file, "Kineis")
  ) %>%
    dplyr::mutate(day = as.Date(timestamp))

  daily_states <- behavior_df %>%
    dplyr::count(track_id, system, day, state, name = "n_state") %>%
    dplyr::group_by(track_id, system, day) %>%
    dplyr::slice_max(n_state, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(track_id, system, day, dominant_state = state) %>%
    tidyr::pivot_wider(names_from = system, values_from = dominant_state) %>%
    dplyr::filter(!is.na(Argos), !is.na(Kineis)) %>%
    dplyr::mutate(agreement = Argos == Kineis)

  behavior_agreement_by_shark <- daily_states %>%
    dplyr::group_by(track_id) %>%
    dplyr::summarise(
      n_matched_days = dplyr::n(),
      n_agree = sum(agreement),
      proportion_agreement = mean(agreement),
      .groups = "drop"
    )

  behavior_state_composition <- behavior_df %>%
    dplyr::count(track_id, system, state, name = "n") %>%
    dplyr::group_by(track_id, system) %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup()

  movement_metrics <- behavior_df %>%
    dplyr::filter(!is.na(x), !is.na(y)) %>%
    dplyr::arrange(track_id, system, timestamp) %>%
    dplyr::group_by(track_id, system) %>%
    dplyr::mutate(
      dx = x - dplyr::lag(x),
      dy = y - dplyr::lag(y),
      step_length = sqrt(dx^2 + dy^2),
      heading = atan2(dy, dx),
      turning_angle = atan2(
        sin(heading - dplyr::lag(heading)),
        cos(heading - dplyr::lag(heading))
      )
    ) %>%
    dplyr::summarise(
      median_log_step_length = median(log1p(step_length), na.rm = TRUE),
      median_absolute_turning_angle = median(abs(turning_angle), na.rm = TRUE),
      .groups = "drop"
    )

  safe_write_csv(daily_states, "Supplemental_Behavior_Matched_Days.csv")
  safe_write_csv(
    behavior_agreement_by_shark,
    "Supplemental_Behavior_Agreement_By_Shark.csv"
  )
  safe_write_csv(
    behavior_state_composition,
    "Supplemental_Behavior_State_Composition.csv"
  )
  safe_write_csv(
    movement_metrics,
    "Supplemental_Behavior_Movement_Metrics.csv"
  )
} else {
  message(
    "\nBehavioral comparison skipped. Add ",
    argos_behavior_file, " and ", kineis_behavior_file,
    " to run Section 09.\n"
  )
}

# =============================================================================
# 10_SUPPLEMENTAL_OUTPUT_INVENTORY_AND_SESSION_INFO
# =============================================================================

all_outputs <- list.files(
  path = ".",
  recursive = TRUE,
  full.names = FALSE
)

output_inventory <- tibble::tibble(file = all_outputs) %>%
  dplyr::filter(
    grepl(
      "Supplement|Figure|KDE|Jaccard|Transmission|Error|Location|SSM|Behavior|Model_Objects",
      file,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::mutate(
    category = dplyr::case_when(
      grepl("Supplemental_Figures|diagnostic", file, ignore.case = TRUE) ~
        "Supplemental figure/diagnostic",
      grepl("Supplemental_Tables|\\.csv$", file, ignore.case = TRUE) ~
        "Table/data output",
      grepl("Main_Figures|Figure", file, ignore.case = TRUE) ~
        "Main figure",
      grepl("Model_Objects|\\.rds$", file, ignore.case = TRUE) ~
        "Saved model",
      TRUE ~ "Other analysis output"
    )
  ) %>%
  dplyr::arrange(category, file)

readr::write_csv(output_inventory, "MASTER_Output_Inventory.csv")
writeLines(capture.output(sessionInfo()), "MASTER_sessionInfo.txt")

cat("\nMaster analysis pipeline complete.\n")
cat("See MASTER_Output_Inventory.csv for the generated files.\n")
