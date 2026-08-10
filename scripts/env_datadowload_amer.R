# Environmental Data download for the American River

library(CDECRetrieve)
library(tidyverse)
library(tidyr)
library(lubridate)
library(dplyr)
library(purrr)
library(here)

# NEED TO HAVE RUN FCNS FROM env_download_fcns.R

project <- here()
source(here(project, '/scripts/env_download_fcns.R'), echo = FALSE)

# Example test using american river stations

#sets the dates to be pulled from cdec, specific to American River

start.date <- "2013-01-01"
end.date <- "2025-06-30"

# set vector of stations
stations.flow.a <- c("AFO", "IST")
stations.wtemp.a <- c("AFO", "AWP", "AWB", "IST")

# Use for single staitons that didn't run correctly
#sta.flow.f <- c("GRL")
#sta.temp.f <- c("FOW")

# use the fcn for flow data
flow_amer <- download_cdec_flow(
  sta_list= stations.flow.a,
  start_date = start.date,
  end_date = end.date,
  dur_code = "E"
)

# use the fcn for temp data
wtemp_amer <- download_cdec_wtemp(
  sta_list= stations.wtemp.a,
  start_date = start.date,
  end_date = end.date,
  dur_code = "E"
)

test.flow <- AFO_flow_data %>% 
  drop_na(datetime) %>% 
  select(-parameter_cd)
test.temp <- AFO_wtemp_data %>% 
  drop_na(datetime) %>% 
  select(datetime, wtemp.f, wtemp.c)

AFO_env <- full_join(test.flow, test.temp, by = join_by(datetime))

test.flow <- IST_flow_data %>% 
  drop_na(datetime) %>% 
  select(-parameter_cd)
test.temp <- IST_wtemp_data %>% 
  drop_na(datetime) %>% 
  select(datetime, wtemp.f, wtemp.c)

IST_env <- full_join(test.flow, test.temp, by = join_by(datetime))


# Write files to data folder
write_csv(AFO_env, "data_raw/amer_env_ORF.csv")
write_csv(IST_env, "data_raw/amer_env_GRL.csv")

write_csv(AWB_wtemp_data, "data_raw/amer_env_AWB_wtemp.csv")
write_csv(AWP_wtemp_data, "data_raw/amer_env_AWP_wtemp.csv")
