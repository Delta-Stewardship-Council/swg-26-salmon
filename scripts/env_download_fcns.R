library(CDECRetrieve)
library(tidyverse)
library(lubridate)
library(dplyr)
library(purrr)

# This is a function to dowload multiple CDEC stations at once as save as dataframes, using
# the 

#########################################################################
# Function to downlaod multiple FLOW stations from CDEC
download_cdec_flow <- function(sta_list,
  start_date,
  end_date,
  sensor_num = 20,
  dur_code = "E",
  assign_global = TRUE) {

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

dat

}, error = function(e) {

warning(paste("Failed for", sta, ":", e$message))
NULL
})
}
)

names(results) <- paste0(sta_list, "_data")

# Create separate objects in global environment
if (assign_global) {
list2env(results, envir = .GlobalEnv)
}

return(results)
}
####################################################################


# Example test using Feather river stations

#sets the dates to be pulled from cdec, specific to Feather River

start.date.f <- "2013-01-01"
end.date.f <- "2025-06-30"

# set vector of stations
stations.flow.f <- c("ORF", "GRL", "FSB", "VON")

# use the fcn
flow_feather <- download_cdec_flow(
  sta_list= stations.flow.f,
  start_date = start.date.f,
  end_date = end.date.f
)
