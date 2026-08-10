# Environmental Data download for the Feather River

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

# Example test using Feather river stations

#sets the dates to be pulled from cdec, specific to Feather River

start.date.f <- "2013-01-01"
end.date.f <- "2025-06-30"

# set vector of stations
stations.flow.f <- c("ORF", "GRL", "FSB", "VON")
stations.wtemp.f <- c("ORF", "FRA", "FTA", "FOW", "GRL", "VON")

# Use for single staitons that didn't run correctly
#sta.flow.f <- c("GRL")
#sta.temp.f <- c("FOW")

# use the fcn for flow data
flow_feather <- download_cdec_flow(
  sta_list= stations.flow.f,
  start_date = start.date.f,
  end_date = end.date.f,
  dur_code = "E"
)

# use the fcn for temp data
wtemp_feather <- download_cdec_wtemp(
  sta_list= stations.wtemp.f,
  start_date = start.date.f,
  end_date = end.date.f,
  dur_code = "E"
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
write_csv(GRL_env, "data_raw/feat_env_GRL.csv")
write_csv(VON_env, "data_raw/feat_env_VON.csv")

write_csv(FOW_wtemp_data, "data_raw/feat_env_FOW_wtemp.csv")
write_csv(FRA_wtemp_data, "data_raw/feat_env_FRA_wtemp.csv")
write_csv(FTA_wtemp_data, "data_raw/feat_env_FTA_wtemp.csv")
write_csv(FSB_flow_data, "data_raw/feat_env_FSB_flow.csv")
