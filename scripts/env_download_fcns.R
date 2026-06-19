library(CDECRetrieve)
library(tidyverse)
library(tidyr)
library(lubridate)
library(dplyr)
library(purrr)

# This is a function to dowload multiple CDEC stations at once as save as dataframes, using
# the 

#########################################################################
# Function to downlaod multiple FLOW stations from CDEC
# A few caveats: 
#1. All stations must have the same sensor code AND event duration
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

dat <- dat %>%
  rename(flow.cfs = parameter_value)
dat

}, error = function(e) {

warning(paste("Failed for", sta, ":", e$message))
NULL
})
}
)

names(results) <- paste0(sta_list, "_flow_data")

# Create separate objects in global environment
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
  sensor_num = 25, #temp w *F
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

dat <- dat %>% 
  rename(wtemp.f=parameter_value) %>% 
  mutate(wtemp.c = ((wtemp.f - 32)*(5/9)))
dat

}, error = function(e) {

warning(paste("Failed for", sta, ":", e$message))
NULL
})
}
)

names(results) <- paste0(sta_list, "_wtemp_data")

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
stations.wtemp.f <- c("ORF", "GRL", "VON")

# use the fcn for flow data
flow_feather <- download_cdec_flow(
  sta_list= stations.flow.f,
  start_date = start.date.f,
  end_date = end.date.f
)

# use the fcn for temp data
wtemp_feather <- download_cdec_wtemp(
  sta_list= stations.wtemp.f,
  start_date = start.date.f,
  end_date = end.date.f
)

test.flow <- ORF_flow_data %>% 
  drop_na(datetime) %>% 
  select(-parameter_cd)
test.temp <- ORF_wtemp_data %>% 
  drop_na(datetime) %>% 
  select(datetime, wtemp.f, wtemp.c)

ORF_env <- full_join(test.flow, test.temp, by = join_by(datetime))

test.flow <- GRL_flow_data %>% 
  drop_na(datetime) %>% 
  select(-parameter_cd)
test.temp <- GRL_wtemp_data %>% 
  drop_na(datetime) %>% 
  select(datetime, wtemp.f, wtemp.c)

GRL_env <- full_join(test.flow, test.temp, by = join_by(datetime))

test.flow <- VON_flow_data %>% 
  drop_na(datetime) %>% 
  select(-parameter_cd)
test.temp <- VON_wtemp_data %>% 
  drop_na(datetime) %>% 
  select(datetime, wtemp.f, wtemp.c)

VON_env <- full_join(test.flow, test.temp, by = join_by(datetime))



# Write files to data folder
write_csv(ORF_env, "data_raw/feat_env_ORF.csv")
write_csv(FSB_flow_data, "data_raw/feat_env_FSB.csv")
write_csv(GRL_env, "data_raw/feat_env_GRL.csv")
write_csv(VON_env, "data_raw/feat_env_VON.csv")

