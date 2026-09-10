# Environmental Data download and QA/QC for the Feather River
# L McCormick updated 2026-09-10

library(CDECRetrieve)
library(tidyverse)
library(tidyr)
library(lubridate)
library(dplyr)
library(purrr)
library(here)
library(ggridges)
library(plotly)

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


# Write RAW data files to data folder
write_csv(ORF_env, "data_raw/feat_env_ORF.csv")
write_csv(GRL_env, "data_raw/feat_env_GRL.csv")
write_csv(VON_env, "data_raw/feat_env_VON.csv")

write_csv(FOW_wtemp_data, "data_raw/feat_env_FOW_wtemp.csv")
write_csv(FRA_wtemp_data, "data_raw/feat_env_FRA_wtemp.csv")
write_csv(FTA_wtemp_data, "data_raw/feat_env_FTA_wtemp.csv")
write_csv(FSB_flow_data, "data_raw/feat_env_FSB_flow.csv")




################################################################
## Data QAQC
################################################################

#############  ORF
ORF_env <- read.csv(here("data_raw/feat_env_ORF.csv")) %>% 
  mutate(agency_cd = "CDEC",
          location_id = "ORF",
        latitude = 39.52,
        longitude = -121.55) %>% 
  mutate(datetime= as.POSIXct(datetime, format= "%Y-%m-%dT%H:%M:%OSZ"))

ggplot(ORF_env, aes(x= datetime, y=wtemp.f))+
  geom_line()

#which(ORF_env$wtemp.f > 100)
#ORF_env$wtemp.f[which(ORF_env$wtemp.f > 100)]

# remove value of >10,000 *F
ORF_env$wtemp.c[which(ORF_env$wtemp.f > 100)] <- NA
ORF_env$wtemp.f[which(ORF_env$wtemp.f > 100)] <- NA

## weird spike of temp data in 2019, where temp goes up to the 90s. I'm not sure whether that's possible,
## so leaving here for now. Removed obvious outliers

ggplot(ORF_env, aes(x= datetime, y=flow.cfs))+
  geom_line()

#ORF_env$flow.cfs[which(ORF_env$flow.cfs < 80)]
# lots of values of <80 cfs, but I don't feel comfortable removing them yet

# remove value of 0 cfs
ORF_env$flow.cfs[which(ORF_env$flow.cfs == 0)] <- NA

# save QAQC file
write_csv(ORF_env, "data_processed/feat_env_ORF_clean.csv")



############## GRL
GRL_env <- read.csv(here("data_raw/feat_env_GRL.csv")) %>% 
  mutate(datetime= as.POSIXct(datetime, format= "%Y-%m-%dT%H:%M:%OSZ"))

ggplot(GRL_env, aes(x= datetime, y=wtemp.f))+
  geom_line()

# GRL_env$wtemp.f[which(GRL_env$wtemp.f<42)]

# remove value of < 42 *F
GRL_env$wtemp.c[which(GRL_env$wtemp.f < 42)] <- NA
GRL_env$wtemp.f[which(GRL_env$wtemp.f < 42)] <- NA

ggplot(GRL_env, aes(x= datetime, y=flow.cfs))+
  geom_line()

# GRL_env$flow.cfs[which(GRL_env$flow.cfs <0)]
# GRL_env$flow.cfs[which(GRL_env$flow.cfs > 150000)] # can't confirm this is an error

GRL_env$flow.cfs[which(GRL_env$flow.cfs < 0)] <- NA

# save QAQC file
write_csv(GRL_env, "data_processed/feat_env_GRL_clean.csv")



################ VON

VON_env <- read.csv(here("data_raw/feat_env_VON.csv")) %>% 
  mutate(agency_cd = "CDEC",
          location_id = "VON",
          latitude = 38.77,
          longitude = -121.60)%>% 
  mutate(datetime= as.POSIXct(datetime, format= "%Y-%m-%dT%H:%M:%OSZ"))

t <- ggplot(VON_env, aes(x= datetime, y=wtemp.f))+
  geom_line()

set.seed(100)
ggplotly(t)

#VON_env$wtemp.f[which(VON_env$wtemp.f < 42)]

# remove value of > 150 *F
VON_env$wtemp.c[which(VON_env$wtemp.f > 150)] <- NA
VON_env$wtemp.f[which(VON_env$wtemp.f > 150)] <- NA
# remove values of < 43
VON_env$wtemp.c[which(VON_env$wtemp.f < 42)] <- NA
VON_env$wtemp.f[which(VON_env$wtemp.f < 42)] <- NA

t <- ggplot(VON_env, aes(x= datetime, y=flow.cfs))+
  geom_line()

VON_env$flow.cfs[which(VON_env$flow.cfs <=0)] <- NA

# save QAQC file
write_csv(VON_env, "data_processed/feat_env_VON_clean.csv")



###########  FOW

FOW_wtemp <- read.csv(here("data_raw/feat_env_FOW_wtemp.csv")) %>% 
  mutate(flow.cfs = NA) %>% 
  mutate(datetime= as.POSIXct(datetime, format= "%Y-%m-%dT%H:%M:%OSZ")) %>% 
  select(agency_cd, location_id, datetime,flow.cfs, latitude, longitude, wtemp.f, wtemp.c)

t <- ggplot(FOW_wtemp, aes(x= datetime, y=wtemp.f))+
  geom_line()

set.seed(100)
ggplotly(t)

# remove value of < 45 *F
FOW_wtemp$wtemp.c[which(FOW_wtemp$wtemp.f <= 45)] <- NA
FOW_wtemp$wtemp.f[which(FOW_wtemp$wtemp.f <= 45)] <- NA

# save QAQC file
write_csv(FOW_wtemp, "data_processed/feat_env_FOW_wtemp_clean.csv")



###########   FRA
FRA_wtemp <- read.csv(here("data_raw/feat_env_FRA_wtemp.csv"))%>% 
  mutate(flow.cfs = NA) %>% 
  mutate(datetime= as.POSIXct(datetime, format= "%Y-%m-%dT%H:%M:%OSZ")) %>% 
  select(agency_cd, location_id, datetime,flow.cfs, latitude, longitude, wtemp.f, wtemp.c)

t <- ggplot(FRA_wtemp, aes(x= datetime, y=wtemp.f))+
  geom_line()

# save QAQC file
write_csv(FRA_wtemp, "data_processed/feat_env_FRA_wtemp_clean.csv")




################  FTA

FTA_wtemp <- read.csv(here("data_raw/feat_env_FTA_wtemp.csv"))%>% 
  mutate(flow.cfs = NA) %>% 
  mutate(datetime= as.POSIXct(datetime, format= "%Y-%m-%dT%H:%M:%OSZ")) %>% 
  select(agency_cd, location_id, datetime,flow.cfs, latitude, longitude, wtemp.f, wtemp.c)

t <- ggplot(FTA_wtemp, aes(x= datetime, y=wtemp.f))+
  geom_line()

set.seed(100)
ggplotly(t)

# remove value of < 45 *F
FTA_wtemp$wtemp.c[which(FTA_wtemp$wtemp.f <= 42)] <- NA
FTA_wtemp$wtemp.f[which(FTA_wtemp$wtemp.f <= 42)] <- NA

# save QAQC file
write_csv(FTA_wtemp, "data_processed/feat_env_FTA_wtemp_clean.csv")



##############  FSB

FSB_flow <- read.csv(here("data_raw/feat_env_FSB_flow.csv"))%>% 
  mutate(wtemp.f = NA,
          wtemp.c= NA) %>% 
  mutate(datetime= as.POSIXct(datetime, format= "%Y-%m-%dT%H:%M:%OSZ")) %>% 
  select(agency_cd, location_id, datetime, flow.cfs, latitude, longitude, wtemp.f, wtemp.c)

t <- ggplot(FSB_flow, aes(x= datetime, y=flow.cfs))+
  geom_line()

set.seed(100)
ggplotly(t)

# remove value of < 0
FSB_flow$flow.cfs[which(FSB_flow$flow.cfs <= 0)] <- NA


# save QAQC file
write_csv(FSB_flow, "data_processed/feat_env_FSB_flow_clean.csv")




# Combine 
FR_env_all <- rbind.data.frame(ORF_env, GRL_env, VON_env, FOW_wtemp, FRA_wtemp,
        FTA_wtemp, FSB_flow) #%>% 
  #mutate(datetime= as.POSIXct(datetime, format= "%Y-%m-%dT%H:%M:%OSZ"))

# order locations
FR_env_all$location_id <- as.factor(FR_env_all$location_id)
FR_env_all$location_id <- ordered(FR_env_all$location_id, levels= c("ORF", "FRA", "FTA", "FOW",
"GRL", "FSB", "VON"))

# Plot example
# flow
fr_flow_all_sta <- ggplot(FR_env_all, aes(x= datetime, y= flow.cfs))+
  geom_line()+
  facet_wrap(FR_env_all$location_id, scales = "free")+
    ylab("Discharge (cfs)")+
    xlab("Date/Time")+
    ggtitle("Feather River Flow")

ggsave(here("figures", "feath_env_flow_data.png"),
       plot = fr_flow_all_sta,
       width = 8, height = 5, dpi = 300)

# temp
fr_wtemp_all_sta <-ggplot(FR_env_all, aes(x= datetime, y= wtemp.c))+
  geom_line()+
  facet_wrap(FR_env_all$location_id)+
  ylab("Water Temperature (C)")+
  xlab("Date/Time")+
  ggtitle("Feather River Water Temp")

ggsave(here("figures", "feath_env_wtemp_data.png"),
       plot = fr_wtemp_all_sta,
       width = 8, height = 5, dpi = 300)



### Save full finali FR data - QAQC version
write_csv(FR_env_all, "data_processed/feat_env_all_clean.csv")





# Test code- probably doesn't work!
# plot alternative - ridgelines (code from amer_exploratory_figures.R)
# nah
ggplot(FR_env_all, aes(x = flow.cfs, y = as.factor(location_id))) +
  geom_density_ridges(fill = "steelblue", alpha = 0.6, scale = 0.9) +
  facet_wrap(facets = FR_env_all$location_id)
  scale_x_continuous(
    breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
    labels = month.abb
  ) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5))) +
  labs(x = NULL, y = "Water year",
       title = "Fall run catch timing by water year")
