library(CDECRetrieve)
library(tidyverse)
library(tidyr)
library(lubridate)
library(dplyr)
library(purrr)
library(here)

# This is a function to dowload multiple CDEC stations at once and save as dataframes, using
# the cdecquery function 

# Also uses a matching function to simultaneously add lat/lon for stations to data

#########################################################################
# Function to downlaod multiple FLOW stations from CDEC
# A few caveats: 
#1. All stations must have the same sensor code. Event duration can be grouped together
download_cdec_flow <- function(sta_list,
  start_date,
  end_date,
  sensor_num = 20,
  dur_code,
  coords_file = here("data_raw", "env_map_coords_2.csv"), # path to coord csv with location_id, lat, lon
  assign_global = TRUE) {

# Load and validate coordinates lookup
if (!is.null(coords_file)) {
coords <- read.csv(coords_file, stringsAsFactors = FALSE)

# Flexible column detection: tolerate varied capitalisation / naming
names(coords) <- tolower(trimws(names(coords)))

# Expect columns: location_id (or station_id / station), latitude (or lat), longitude (or lon / long)
id_col  <- grep("^(location_id|station.id|station)$", names(coords), value = TRUE)[1]
lat_col <- grep("^(latitude|lat)$",                   names(coords), value = TRUE)[1]
lon_col <- grep("^(longitude|lon|long)$",             names(coords), value = TRUE)[1]

if (any(is.na(c(id_col, lat_col, lon_col)))) {
warning("coords_file must contain columns for station ID, latitude, and longitude. ",
"Recognised names: location_id/station.id/station, latitude/lat, longitude/lon/long. ",
"Coordinates will NOT be joined.")
coords <- NULL
} else {
coords <- coords[, c(id_col, lat_col, lon_col)]
names(coords) <- c("location_id", "latitude", "longitude")
coords$location_id <- toupper(trimws(coords$location_id))
}
} else {
coords <- NULL
}

# download loop
results <- map(
sta_list,
function(sta) {

message(paste("Downloading", sta))

tryCatch({

dat <- cdec_query(
station    = sta,
sensor_num = sensor_num,
dur_code   = dur_code,
start_date = start_date,
end_date   = end_date
)

if (nrow(dat) == 0) {
warning(paste("No data returned for", sta))
}

dat <- dat %>%
  rename(flow.cfs = parameter_value)

# --- Join coordinates -------------------------------------------------
if (!is.null(coords)) {
dat <- dat %>%
mutate(location_id = toupper(trimws(location_id))) %>%
left_join(coords, by = "location_id")

if (all(is.na(dat$latitude))) {
warning(paste("No coordinate match found for station", sta,
"— check that the location_id in your CSV matches CDEC station codes."))
}
}

dat

}, error = function(e) {
warning(paste("Failed for", sta, ":", e$message))
NULL
})
}
)

names(results) <- paste0(sta_list, "_flow_data")

if (assign_global) {
list2env(results, envir = .GlobalEnv)
}

return(results)
}
####################################################################

#########################################################################
# Function to downlaod multiple TEMP stations from CDEC
# A few caveats: 
#1. All stations must have the same sensor code AND event duration

download_cdec_wtemp <- function(sta_list,
  start_date,
  end_date,
  sensor_num  = 25,       # wtemp in °F
  dur_code,
  coords_file = here("data_raw", "env_map_coords_2.csv"), # path to coord csv with location_id, lat, lon
  assign_global = TRUE) {

# Load and validate coordinates lookup
if (!is.null(coords_file)) {
coords <- read.csv(coords_file, stringsAsFactors = FALSE)

# Flexible column detection: tolerate varied capitalisation / naming
names(coords) <- tolower(trimws(names(coords)))

# Expect columns: location_id (or station_id / station), latitude (or lat), longitude (or lon / long)
id_col  <- grep("^(location_id|station.id|station)$", names(coords), value = TRUE)[1]
lat_col <- grep("^(latitude|lat)$",                   names(coords), value = TRUE)[1]
lon_col <- grep("^(longitude|lon|long)$",             names(coords), value = TRUE)[1]

if (any(is.na(c(id_col, lat_col, lon_col)))) {
warning("coords_file must contain columns for station ID, latitude, and longitude. ",
"Recognised names: location_id/station.id/station, latitude/lat, longitude/lon/long. ",
"Coordinates will NOT be joined.")
coords <- NULL
} else {
coords <- coords[, c(id_col, lat_col, lon_col)]
names(coords) <- c("location_id", "latitude", "longitude")
coords$location_id <- toupper(trimws(coords$location_id))
}
} else {
coords <- NULL
}

# download loop
results <- map(
sta_list,
function(sta) {

message(paste("Downloading", sta))

tryCatch({

dat <- cdec_query(
station    = sta,
sensor_num = sensor_num,
dur_code   = dur_code,
start_date = start_date,
end_date   = end_date
)

if (nrow(dat) == 0) {
warning(paste("No data returned for", sta))
}

dat <- dat %>%
rename(wtemp.f = parameter_value) %>%
mutate(wtemp.c = (wtemp.f - 32) * (5 / 9))

# --- Join coordinates -------------------------------------------------
if (!is.null(coords)) {
dat <- dat %>%
mutate(location_id = toupper(trimws(location_id))) %>%
left_join(coords, by = "location_id")

if (all(is.na(dat$latitude))) {
warning(paste("No coordinate match found for station", sta,
"— check that the location_id in your CSV matches CDEC station codes."))
}
}

dat

}, error = function(e) {
warning(paste("Failed for", sta, ":", e$message))
NULL
})
}
)

names(results) <- paste0(sta_list, "_wtemp_data")

if (assign_global) {
list2env(results, envir = .GlobalEnv)
}

return(results)
}
####################################################################


