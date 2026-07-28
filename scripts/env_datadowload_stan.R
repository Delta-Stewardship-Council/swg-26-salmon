# Environmental Data download for the Stanislaus River

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


#sets the dates to be pulled from cdec

start.date <- "2013-01-01"
end.date <- "2025-06-30"

# set vector of stations
# hourly
stations.wtemp.sh <- c("GDC", "KFS", "OBS", "ORA", "GMB", "JMP", "SRJ", "RPN")

# daily
stations.flow.sd <- c("GDW")

# event (preferred)
stations.flow.se <- c("RIP", "KOT")
stations.wtemp.se <- c("SOK", "RIP")

# Use for single staitons that didn't run correctly
#sta.flow.f <- c("GRL")
#sta.temp.f <- c("FOW")

# use the fcn for flow data- EVENT
flow_stan <- download_cdec_flow(
  sta_list= stations.flow.se,
  start_date = start.date,
  end_date = end.date,
  dur_code = "E"
)

# use the fcn for temp data -EVENT
wtemp_stan <- download_cdec_wtemp(
  sta_list= stations.wtemp.se,
  start_date = start.date,
  end_date = end.date,
  dur_code = "E"
)

test.flow <- RIP_flow_data %>% 
  drop_na(datetime) %>% 
  select(-parameter_cd)
test.temp <- RIP_wtemp_data %>% 
  drop_na(datetime) %>% 
  select(datetime, wtemp.f, wtemp.c)

RIP_env <- full_join(test.flow, test.temp, by = join_by(datetime))

# use the fcn for flow data- DAILY
flow_stan <- download_cdec_flow(
  sta_list= stations.flow.sd,
  start_date = start.date,
  end_date = end.date,
  dur_code = "D"
)

# use the fcn for temp data -HOURLY
wtemp_stan <- download_cdec_wtemp(
  sta_list= stations.wtemp.sh,
  start_date = start.date,
  end_date = end.date,
  dur_code = "E"
)



# Write files to data folder
write_csv(RIP_env, "data_raw/stan_env_RIP.csv")
#write_csv(IST_env, "data_raw/amer_env_GRL.csv")

write_csv(KOT_flow_data, "data_raw/stan_env_KOT_flow.csv")
write_csv(SOK_wtemp_data, "data_raw/stan_env_SOK_wtemp.csv")
write_csv(GDW_flow_data, "data_raw/stan_env_GDW_flow.csv")
write_csv(GDC_wtemp_data, "data_raw/stan_env_GDC_wtemp.csv")
write_csv(KFS_wtemp_data, "data_raw/stan_env_KFS_wtemp.csv")
write_csv(OBS_wtemp_data, "data_raw/stan_env_OBS_wtemp.csv")
write_csv(ORA_wtemp_data, "data_raw/stan_env_ORA_wtemp.csv")
write_csv(GMB_wtemp_data, "data_raw/stan_env_GMB_wtemp.csv")
write_csv(JMP_wtemp_data, "data_raw/stan_env_JMP_wtemp.csv")
write_csv(SRJ_wtemp_data, "data_raw/stan_env_SRJ_wtemp.csv")
write_csv(RPN_wtemp_data, "data_raw/stan_env_RPN_wtemp.csv")
