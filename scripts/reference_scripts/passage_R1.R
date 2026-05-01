### Code from Burford 2025 for estimating survival

library(tidyverse)
library(lubridate)
library(ggmap)
library(mgcv)
library(performance)
library(DHARMa)
library(mgcViz)

setwd("~/Downloads")


# fit, diagnose, and interpret passage probability models
# this code produces:
# Figure 3, Figure 4, Figure 5
# Figure S1, Figure S2, Figure S6, Figure S7

# read in acoustic telemetry data
# response = 0, 1 successful passage from release to freshwater estuary exit
# predictors = migrant condition, migration distance, environmental factors
fish_red <- read.csv("passage_R1_data.csv")

# make map of case study area
# genlocs are estuary exit and ocean entry
# usgs are USGS gauges
# rbdd is red bluff diversion dam
genlocs <- data.frame(latitude = c(38.04337, 38.04123, 37.82794, 37.82126),
                      longitude = c(-122.1234, -122.1238, -122.4617, -122.4692),
                      location = c("estuary exit", "estuary exit", "ocean entry",  "ocean entry"))
usgs <- data.frame(latitude = 38.45601954,
                   longitude = -121.5013437,
                   location = "USGS gauges")
rbdd <- data.frame(latitude = 40.153465,
                   longitude = -122.202709,
                   location = "RBDD/release")
cdec <- data.frame(latitude = 40.28848836,
                   longitude = -122.1866645,
                   location = "USGS gauges")

# bind all these together
genlocs <- bind_rows(genlocs,usgs,rbdd,cdec)

# gather release locations from fish_red
# bind genlocs to this
rel_map <- fish_red %>%
  group_by(release_latitude,release_longitude) %>%
  slice(1) %>%
  rename(latitude = release_latitude,
         longitude = release_longitude) %>%
  select(latitude,longitude) %>%
  mutate(location = "release") %>%
  bind_rows(genlocs)

# get ggmap - will need to register with stadiamaps
# insert access key in below quotes
register_stadiamaps("")
myLocation <- c(-123,38,-121,40.5)
myMap <- get_stadiamap(myLocation, zoom = 8,maptype = "stamen_toner",crop = FALSE)

# order locations for consistency
rel_map$location = factor(rel_map$location, levels=c("release",
                                                     "RBDD/release",
                                                     "USGS gauges",
                                                     "estuary exit",
                                                     "ocean entry"))

# Figure 3 - map of case study area
myMap %>%
  ggmap(base_layer = ggplot(aes(longitude,latitude,col=location),data=rel_map)) +
  geom_point(size=5,alpha=1) +
  scale_y_continuous(breaks=seq(38,41,0.5),expand = c(0,0)) +
  scale_x_continuous(breaks=seq(-123.5,-121,0.5),expand = c(0,0)) +
  scale_color_viridis_d() +
  xlab(c("longitude (°W)")) +
  ylab(c("latitude (°N)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.26,0.6),legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

#ggsave("Figure2.png",width = 3.5,height = 5)

# fit, diagnose, and interpret time and temperature GAMS
# both models have satisfactory diagnostics no matter how many knots used (5-9)
# more than 5 knots gives peculiar predictions
# smooth or parametric interactions between time or temperature and flow give peculiar predictions
# thus, using 5 knots and no interactions in final models below

# time GAM for juvenile Chinook salmon passage success probability
g1 <- bam(survival ~
            s(cond,bs="cr",k=5) +
            s(release_river_km,bs="cr",k=5) +
            s(water_year_doy,bs="cr",k=5) +
            s(flow_cms,bs="cr",k=5),
          method = "fREML",
          select = T,
          family = binomial(link = "logit"),
          data = fish_red)

# temperature GAM for juvenile Chinook salmon passage success probability
g1t <- bam(survival ~
             s(cond,bs="cr",k=5) +
             s(release_river_km,bs="cr",k=5) +
             s(temp_c,bs="cr",k=5) +
             s(flow_cms,bs="cr",k=5),
           method = "fREML",
           select = T,
           family = binomial(link = "logit"),
           data = fish_red)

# check smooth functions and concurvity for both models
plot.gam(g1,pages=1)
g_concurv <- data.frame(concurvity(g1,full = FALSE)$observed) %>%
  select(-c(para)) %>%
  slice_tail(n=4)
g_concurv
#                        s.cond. s.release_river_km. s.water_year_doy. s.flow_cms.
#s(cond)             1.00000000          0.01411787        0.05440623 0.005609477
#s(release_river_km) 0.04716970          1.00000000        0.48527252 0.028935363
#s(water_year_doy)   0.04722393          0.35150618        1.00000000 0.017340951
#s(flow_cms)         0.01966619          0.20996449        0.31544573 1.000000000

plot.gam(g1t,pages=1)
g_concurv <- data.frame(concurvity(g1t,full = FALSE)$observed) %>%
  select(-c(para)) %>%
  slice_tail(n=4)
g_concurv
#                        s.cond. s.release_river_km.  s.temp_c. s.flow_cms.
#s(cond)             1.00000000         0.003308956 0.01771554  0.00591177
#s(release_river_km) 0.04702094         1.000000000 0.30488929  0.02521799
#s(temp_c)           0.01691686         0.220502513 1.00000000  0.32674168
#s(flow_cms)         0.01959056         0.179271529 0.40773504  1.00000000

# check residual diagnostics of both models
# Figure S1 - residual diagnostics of time GAM
testUniformity(g1)
testOutliers(g1)
testDispersion(g1)
testQuantiles(g1)
testZeroInflation(g1)

# Figure S2 - residual diagnostics of temperature GAM
testUniformity(g1t)
testOutliers(g1t)
testDispersion(g1t)
testQuantiles(g1t)
testZeroInflation(g1t)

# summary statistics for both models
# this info is presented in the results
summary(g1)
sum(influence(g1))
df.residual(g1)

summary(g1t)
sum(influence(g1t))
df.residual(g1t)

# compare maximum effect size of day of water year and flow for time GAM
# starting with day of water year
# build data to predict over
# mean migrant condition, median migration distance
# full range of everything else
newdat <- expand.grid(cond = mean(fish_red$cond,na.rm=T),
                      release_river_km = median(fish_red$release_river_km),
                      water_year_doy = seq(min(fish_red$water_year_doy),max(fish_red$water_year_doy),1),
                      flow_cms = c(seq(min(fish_red$flow_cms),max(fish_red$flow_cms),length.out=100)))

# make predictions and add to newdata
prdat <-  predict.gam(g1,newdat,se.fit = T)
newdat$fit <- prdat$fit

# gather inverse link function to convert from logit to probability
ilink <- family(g1)$linkinv

# gather minimum effect size of day for each flow
temp <- newdat %>%
  mutate(fit = ilink(fit)) %>%
  group_by(flow_cms) %>%
  slice_min(fit,n=1) %>%
  ungroup() %>%
  mutate(what="min")

# bind to maximum effect size of day for each flow
doyeff <- newdat %>%
  mutate(fit = ilink(fit)) %>%
  group_by(flow_cms) %>%
  slice_max(fit,n=1) %>%
  ungroup() %>%
  mutate(what="max") %>%
  bind_rows(temp) %>%
  arrange(flow_cms) %>%
  select(flow_cms,what,fit) %>%
  pivot_wider(names_from = what,values_from = fit) %>%
  mutate(eff = max-min)

# repeat this process for flow
newdat <- expand.grid(cond = mean(fish_red$cond,na.rm=T),
                      release_river_km = median(fish_red$release_river_km),
                      water_year_doy = seq(min(fish_red$water_year_doy),max(fish_red$water_year_doy),1),
                      flow_cms = c(seq(min(fish_red$flow_cms),max(fish_red$flow_cms),length.out=100)))

prdat <-  predict.gam(g1,newdat,se.fit = T)
newdat$fit <- prdat$fit
ilink <- family(g1)$linkinv

# gather min effect size of flow for each day
temp <- newdat %>%
  mutate(fit = ilink(fit)) %>%
  group_by(water_year_doy) %>%
  slice_min(fit,n=1) %>%
  ungroup() %>%
  mutate(what="min")

# bind to max
floweff <- newdat %>%
  mutate(fit = ilink(fit)) %>%
  group_by(water_year_doy) %>%
  slice_max(fit,n=1) %>%
  ungroup() %>%
  mutate(what="max") %>%
  bind_rows(temp) %>%
  arrange(water_year_doy) %>%
  select(water_year_doy,what,fit) %>%
  pivot_wider(names_from = what,values_from = fit) %>%
  mutate(eff = max-min)

# calculate how many more times greater the max effect size of flow was than that of day
# this info is presented in discussion
max(floweff$eff)/max(doyeff$eff)
# 1.574324

# Figure S7 - day of water year effect size vs. flow
doyeff %>%
  ggplot(aes(flow_cms,eff)) +
  geom_line(linewidth=3) +
  geom_vline(xintercept = median(fish_red$flow_cms),lty=2) +
  geom_vline(xintercept = mean(fish_red$flow_cms)) +
  theme_classic(base_size = 14) +
  ylab("day of water year effect size (Pr)") +
  xlab(expression(paste("flow (",m^3,s^-1,")")))

#ggsave("FigureS7.png",width = 4,height = 4)

# plot effect of migrant condition and migration distance from time GAM
# build data to predict over
# median migration distance, median day of water year, mean flow
# full range of migrant condition
newdat <- expand.grid(cond = seq(min(fish_red$cond,na.rm=T),max(fish_red$cond,na.rm=T),length.out=1000),
                      release_river_km = median(fish_red$release_river_km),
                      water_year_doy = median(fish_red$water_year_doy),
                      flow_cms = mean(fish_red$flow_cms))

# make predictions and add to newdata
prdat <-  predict.gam(g1,newdat,se.fit = T)
newdat$fit <- prdat$fit
newdat$conf <- prdat$se.fit*2
newdat$lwr <- newdat$fit-newdat$conf
newdat$upr <- newdat$fit+newdat$conf

# gather inverse link function and transform from logit to probability
ilink <- family(g1)$linkinv
newdat <- newdat %>%
  mutate(fit = ilink(fit),
         upr = ilink(upr),
         lwr = ilink(lwr))

# Figure S6A - smoothed function for migrant condition
newdat %>%
  ggplot(aes(cond,fit)) +
  geom_rug(aes(cond),fish_red,alpha=0.1,inherit.aes = F) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr),alpha=0.5) +
  theme_classic(base_size = 14) +
  xlab(expression(paste("migrant condition (100 x ",g~cm^-3,")"))) +
  ylab("probability of sucessful passage ± 95% CI")

#ggsave("FigureS6A.png",width = 4,height = 4.5)

# repeat for migration distance
# use mean migrant condition and full range of migration distance
newdat <- expand.grid(cond = mean(fish_red$cond,na.rm=T),
                      release_river_km = seq(min(fish_red$release_river_km,na.rm=T),max(fish_red$release_river_km,na.rm=T),length.out=1000),
                      water_year_doy = median(fish_red$water_year_doy),
                      flow_cms = mean(fish_red$flow_cms))

prdat <-  predict.gam(g1,newdat,se.fit = T)
newdat$fit <- prdat$fit
newdat$conf <- prdat$se.fit*2
newdat$lwr <- newdat$fit-newdat$conf
newdat$upr <- newdat$fit+newdat$conf

ilink <- family(g1)$linkinv
newdat <- newdat %>%
  mutate(fit = ilink(fit),
         upr = ilink(upr),
         lwr = ilink(lwr))

# Figure S6B - smoothed function for migrant condition
newdat %>%
  ggplot(aes(release_river_km,fit)) +
  geom_rug(aes(release_river_km),fish_red,alpha=0.1,inherit.aes = F) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr),alpha=0.5) +
  theme_classic(base_size = 14) +
  xlab("migration distance (km)") +
  ylab("probability of sucessful passage ± 95% CI")

#ggsave("FigureS6B.png",width = 4,height = 4.5)

# examine temporal pattern of passage success with respect to natural migration phenology
# read in RBDD initiation data
mig_start_full <- read.csv("migration_R1_data.csv")

# associate water years in RBDD data with years that had 365 or 366 days
max_day_prev <- data.frame(water_year = seq(2007,2019,1),
                           max_doy_prev = c(365,365,366,365,365,365,366,365,365,365,366,365,365))

# calculate the day of cumulative passage closest to 10, 50, and 90% of annual passage
# these are referred to as "initiation thresholds" in the manuscript (e.g., Figure 4)
# thresholds determined for all runs of both lifestages each year

# start with 10% threshold
# smolts of all runs first
mig_start <- mig_start_full %>%
  group_by(run,water_year_n) %>%
  mutate(passage_tot = sum(smolt,na.rm=T),
         passage_cum = cumsum(smolt),
         passage_perc = (passage_cum/passage_tot)*100) %>%
  filter(abs(passage_perc-10)==min(abs(passage_perc-10))) %>%
  ungroup() %>%
  group_by(run,water_year_n,passage_perc) %>%
  arrange(passage_perc) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(lifestage = "smolt")

# bind to max day
mig_start <- bind_rows(mig_start,max_day_prev)

# repeat for fry of all runs
temp <- mig_start_full %>%
  group_by(run,water_year_n) %>%
  mutate(passage_tot = sum(fry,na.rm=T),
         passage_cum = cumsum(fry),
         passage_perc = (passage_cum/passage_tot)*100) %>%
  filter(abs(passage_perc-10)==min(abs(passage_perc-10))) %>%
  ungroup() %>%
  group_by(run,water_year_n,passage_perc) %>%
  arrange(passage_perc) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(lifestage = "fry")

# bind to smolt data
mig_start <- bind_rows(mig_start,temp)

# filter out missing data
mig_start <- filter(mig_start,!is.na(lifestage))

# repeat this process for 50% threshold
mig_mid <- mig_start_full %>%
  group_by(run,water_year_n) %>%
  mutate(passage_tot = sum(smolt,na.rm=T),
         passage_cum = cumsum(smolt),
         passage_perc = (passage_cum/passage_tot)*100) %>%
  filter(abs(passage_perc-50)==min(abs(passage_perc-50))) %>%
  ungroup() %>%
  group_by(run,water_year_n,passage_perc) %>%
  arrange(passage_perc) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(lifestage = "smolt")

mig_mid <- bind_rows(mig_mid,max_day_prev)

temp <- mig_start_full %>%
  group_by(run,water_year_n) %>%
  mutate(passage_tot = sum(fry,na.rm=T),
         passage_cum = cumsum(fry),
         passage_perc = (passage_cum/passage_tot)*100) %>%
  filter(abs(passage_perc-50)==min(abs(passage_perc-50))) %>%
  ungroup() %>%
  group_by(run,water_year_n,passage_perc) %>%
  arrange(passage_perc) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(lifestage = "fry")

mig_mid <- bind_rows(mig_mid,temp)
mig_mid <- filter(mig_mid,!is.na(lifestage))

# repeat for 90% threshold
mig_end <- mig_start_full %>%
  group_by(run,water_year_n) %>%
  mutate(passage_tot = sum(smolt,na.rm=T),
         passage_cum = cumsum(smolt),
         passage_perc = (passage_cum/passage_tot)*100) %>%
  filter(abs(passage_perc-90)==min(abs(passage_perc-90))) %>%
  ungroup() %>%
  group_by(run,water_year_n,passage_perc) %>%
  arrange(passage_perc) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(lifestage = "smolt")

mig_end <- bind_rows(mig_end,max_day_prev)

temp <- mig_start_full %>%
  group_by(run,water_year_n) %>%
  mutate(passage_tot = sum(fry,na.rm=T),
         passage_cum = cumsum(fry),
         passage_perc = (passage_cum/passage_tot)*100) %>%
  filter(abs(passage_perc-90)==min(abs(passage_perc-90))) %>%
  ungroup() %>%
  group_by(run,water_year_n,passage_perc) %>%
  arrange(passage_perc) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(lifestage = "fry")

mig_end <- bind_rows(mig_end,temp)
mig_end <- filter(mig_end,!is.na(lifestage))

# for subsequent plotting, gather median 10, 50, and 90% thresholds by run and lifestage
# these will be vertical lines in plot, so defining y-axis values for endpoints (upr, lwr)
# fry lines will be relative to upper x-axis, while smolt lines will relative to lower x-axis
perc_10 <- mig_start %>%
  group_by(run,lifestage) %>%
  summarise(water_year_doy = median(water_year_doy)) %>%
  ungroup() %>%
  mutate(upr = ifelse(lifestage=="fry",0.3,0.125),
         lwr = ifelse(lifestage=="fry",0.125,0),
         threshold = "10%")

perc_50 <- mig_mid %>%
  group_by(run,lifestage) %>%
  summarise(water_year_doy = median(water_year_doy)) %>%
  ungroup() %>%
  mutate(upr = ifelse(lifestage=="fry",0.3,0.125),
         lwr = ifelse(lifestage=="fry",0.125,0),
         threshold = "50%")

perc_90 <- mig_end %>%
  group_by(run,lifestage) %>%
  summarise(water_year_doy = median(water_year_doy)) %>%
  ungroup() %>%
  mutate(upr = ifelse(lifestage=="fry",0.3,0.125),
         lwr = ifelse(lifestage=="fry",0.125,0),
         threshold = "90%")

# bind median passage thresholds together
perc_ <- bind_rows(perc_10,perc_50,perc_90)

# plot day of water year smoothed function from time GAM against initiation data and thresholds
# build data to predict over
# median migration distance, mean migrant condition
# median and mean flow
# full range of day of water year
newdat <- expand.grid(cond = mean(fish_red$cond,na.rm=T),
                      release_river_km = median(fish_red$release_river_km),
                      water_year_doy = seq(min(fish_red$water_year_doy),max(fish_red$water_year_doy),1),
                      flow_cms = c(round(median(fish_red$flow_cms)),round(mean(fish_red$flow_cms))))

# make predictions and add predictions to newdat
prdat <-  predict.gam(g1,newdat,se.fit = T)
newdat$fit <- prdat$fit
newdat$conf <- prdat$se.fit*2
newdat$lwr <- newdat$fit-newdat$conf
newdat$upr <- newdat$fit+newdat$conf

# use inverse link function to convert from logit to probability
ilink <- family(g1)$linkinv
newdat <- newdat %>%
  mutate(fit = ilink(fit),
         upr = ilink(upr),
         lwr = ilink(lwr))

# sum initiation over full RBDD timeseries by day for each lifestage
# these are going to be plotted as ribbons above (fry) and below (smolt) passage probability
# plotting robbins on 0 (min), 0.1 (max) scale
# start with fry, which will be plotted relative to top x-axis
fry <- mig_start_full %>%
  select(-c(smolt,passage)) %>%
  group_by(run,water_year_doy) %>%
  summarise(fry = sum(fry)) %>%
  ungroup() %>%
  group_by(run) %>%
  mutate(tot_fry = ((fry/max(fry))*-0.10)+0.3) %>%
  ungroup()

# repeat for smolts, which will be plotted relative to bottom x-axis
smolt <- mig_start_full %>%
  select(-c(fry,passage)) %>%
  group_by(run,water_year_doy) %>%
  summarise(smolt = sum(smolt)) %>%
  ungroup() %>%
  group_by(run) %>%
  mutate(tot_smolt = (smolt/max(smolt))*0.10) %>%
  ungroup()

# make water year calendar date conversion data
# start with full sequence of dates in dataset
wy <- data.frame(Date = as.POSIXlt(seq(as.Date(min(fish_red$mid_mig_date)), as.Date(max(fish_red$mid_mig_date)), by="1 days")))

# water years straddle calendar years, beginning 10/1 of one year and ending 9/30 the next
# water years are named for the calendar year that they end in
# for each date, determine water year, week, and day
wy <- wy %>%
  mutate(year=year(Date),
         year_next = year+1) %>%
  group_by(year) %>%
  mutate(oct_1 = as.Date(str_c(year,"-",10,"-",1)),
         sept_30 = as.Date(str_c(year,"-",9,"-",30)),
         jan_1 = as.Date(str_c(year_next,"-",1,"-",1)),
         water_year = case_when(Date<oct_1~year,
                                Date>=oct_1~year_next)) %>%
  ungroup() %>%
  mutate(day_diff = jan_1-oct_1,
         temp_date = Date+day_diff,
         water_year_week = week(temp_date),
         water_year_doy = yday(temp_date)) %>%
  select(Date,water_year,water_year_doy)

# make labels for subsequent plots
# just gather a single year's data
# this label was for date range of RBDD initiation data
wy_lab <- wy %>%
  filter(water_year_doy%in%seq(min(mig_start_full$water_year_doy),max(mig_start_full$water_year_doy),50)) %>%
  mutate(d = str_sub(Date,-5,-1)) %>%
  filter(water_year==2019)

# order all data going into plot for consistency
smolt$run = factor(smolt$run, levels=c("winter",
                                       "fall",
                                       "latefall",
                                       "spring"))
fry$run = factor(fry$run, levels=c("winter",
                                   "fall",
                                   "latefall",
                                   "spring"))
perc_$run = factor(perc_$run, levels=c("winter",
                                       "fall",
                                       "latefall",
                                       "spring"))
perc_$threshold = factor(perc_$threshold, levels=c("10%",
                                                   "50%",
                                                   "90%"))

# Figure 4 - passage probability vs. natural migration phenology
newdat %>%
  ggplot(aes(water_year_doy,fit,by=factor(flow_cms))) +
  geom_line(aes(water_year_doy,tot_fry,by=run),fry) +
  geom_ribbon(aes(water_year_doy,tot_fry,ymin = tot_fry, ymax = 0.3,fill=run),
              fry,alpha=0.5,inherit.aes = F) +
  geom_line(aes(water_year_doy,tot_smolt,by=run),smolt) +
  geom_ribbon(aes(water_year_doy,tot_smolt,ymin = 0, ymax = tot_smolt,fill=run),
              smolt,alpha=0.5,inherit.aes = F) +
  geom_segment(aes(x=water_year_doy,xend=water_year_doy,y=upr,yend=lwr,col=run,lty=threshold),
               perc_,inherit.aes = F) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),alpha=0.7) +
  geom_line() +
  scale_x_continuous(limits = c(min(mig_start_full$water_year_doy),max(mig_start_full$water_year_doy)),
                     breaks = seq(min(mig_start_full$water_year_doy),max(mig_start_full$water_year_doy),50),
                     sec.axis = dup_axis(labels = c(wy_lab$d),name = "calendar date"),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),breaks=c(0,0.1,0.2,0.3)) +
  scale_fill_viridis_d(end = 0.9) +
  scale_color_viridis_d(end = 0.9) +
  ylab(c("passage (Pr) ± 95% CI")) +
  xlab(c("day of water year")) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_text(face = "bold",size = 11),
        legend.position = "bottom",
        legend.box.spacing = unit(0.2,"cm")) +
  guides(fill=guide_legend(title="population",direction = "vertical",nrow = 1,keywidth = 0.7),
         col=guide_legend(title="population",direction = "vertical",nrow = 1,keywidth = 0.7), 
         lty=guide_legend(title="initiation threshold",direction = "vertical",nrow = 1,keywidth = 1))

#ggsave("Figure4.png",height = 5,width = 5)

# examine how passage probability varies with flow and day of water year
# build data to predict over
# median migration distance, mean migrant condition
# full range of flow and day of water year
newdat <- expand.grid(cond = mean(fish_red$cond,na.rm=T),
                      release_river_km = median(fish_red$release_river_km),
                      water_year_doy = seq(min(fish_red$water_year_doy),max(fish_red$water_year_doy),7),
                      flow_cms = seq(min(fish_red$flow_cms)-10,max(fish_red$flow_cms)+10,100))

# make predictions and add to newdat
prdat <-  predict.gam(g1,newdat,se.fit = T)
newdat$fit <- prdat$fit
newdat$conf <- prdat$se.fit*2
newdat$lwr <- newdat$fit-newdat$conf
newdat$upr <- newdat$fit+newdat$conf

# use inverse link function to convert from logit to probability
# if the fit is less than ± 95% confidence interval, make NA
# not many other ways to show uncertainty with this plot type
newdat <- newdat %>%
  mutate(fit = ilink(fit),
         upr = ilink(upr),
         lwr = ilink(lwr),
         conf = (upr-lwr)/2,
         fit = ifelse(conf<fit,fit,NA))

# make water year labels for subsequent plots
# note this label is only for the date range of telemetry data
# last label was for date range of RBDD initiation data
wy_lab <- wy %>%
  filter(water_year_doy%in%seq(min(fish_red$water_year_doy),max(fish_red$water_year_doy),50)) %>%
  mutate(d = str_sub(Date,-5,-1)) %>%
  filter(water_year==2019)

# Figure 5A - passage probability with respect to flow and time
newdat %>%
  ggplot(aes(water_year_doy,flow_cms)) +
  geom_raster(aes(fill = fit), interpolate=F) +
  scale_fill_gradientn(colours=c("red","blanchedalmond","blue"),
                       name=expression(bold(atop("passage", paste("(Pr)")))),
                       limits = c(0,0.4),breaks = seq(0,0.4,0.1),na.value = "grey70") +
  geom_point(aes(water_year_doy,flow_cms,pch=origin),data=fish_red,alpha=0.2) +
  scale_x_continuous(limits = c(min(newdat$water_year_doy)-10,max(newdat$water_year_doy)+10),
                     expand = c(0,0),
                     breaks = seq(min(fish_red$water_year_doy),max(fish_red$water_year_doy),50),
                     sec.axis = dup_axis(labels = c(wy_lab$d),name = "calendar date")) +
  scale_y_continuous(limits = c(min(newdat$flow_cms)-100,max(newdat$flow_cms)+100),
                     expand = c(0,0),
                     breaks = seq(500,2500,500)) +
  geom_contour(aes(z=fit),col="black",breaks = c(0.1),lty=2) +
  theme_classic(base_size = 15) +
  theme(legend.title = element_text(face = "bold",size = 12),
        legend.position = "bottom",
        legend.box.spacing = unit(0.1,"cm")) +
  guides(shape=guide_legend(keywidth = 0,direction = "vertical",nrow = 1)) +
  xlab(c("day of water year")) +
  ylab(expression(paste("flow "~"(",m^3,s^-1,")")))

#ggsave("Figure5A.png",height = 5,width = 4.5)

# predict how passage probability changes with 10,000 cfs flow addition
# just predicting passage probability from 10,000 cfs above current values in newdat
newdat <- mutate(newdat,flow_cms = flow_cms+(10000/35.315))

# make predictions and add to newdat
prdat <-  predict.gam(g1,newdat,se.fit = T)
newdat$fit_inc <- prdat$fit
newdat$conf_inc <- prdat$se.fit*2
newdat$lwr_inc <- newdat$fit_inc-newdat$conf_inc
newdat$upr_inc <- newdat$fit_inc+newdat$conf_inc

# convert to probability
# if the ± 95% confidence interval increase was greater than the fit increase, make NA
# same way of showing uncertainty as last plot
newdat <- newdat %>%
  mutate(fit_inc = ilink(fit_inc),
         upr_inc = ilink(upr_inc),
         lwr_inc = ilink(lwr_inc),
         conf_inc = (upr_inc-lwr_inc)/2,
         fit_inc = ifelse(conf_inc<fit_inc,fit_inc,NA))

# determine passage probability difference
# since made predictions for 10,000 cfs addition, need to drop x-axis back down for plotting
newdat <- newdat %>%
  mutate(delta_surv = fit_inc-fit,
         flow_cms = flow_cms-(10000/35.315))
  
# Figure 5C - effect of 10,000 cfs flow addition
newdat %>%
  ggplot(aes(water_year_doy,flow_cms)) +
  geom_raster(aes(fill = delta_surv), interpolate=F) +
  scale_fill_gradientn(colours=c("red","blanchedalmond","blue"),
                       name=expression(bold(atop("passage"~Delta~"(Pr)",paste("+283"~m^3,s^-1)))),
                       limits = c(0,0.15),breaks = c(0,0.05,0.1,0.15),na.value = "grey70") +
  scale_x_continuous(limits = c(min(newdat$water_year_doy)-10,max(newdat$water_year_doy)+10),
                     expand = c(0,0),
                     breaks = seq(min(fish_red$water_year_doy),max(fish_red$water_year_doy),50),
                     sec.axis = dup_axis(labels = c(wy_lab$d),name = "calendar date")) +
  scale_y_continuous(limits = c(min(newdat$flow_cms)-100,max(newdat$flow_cms)+100),
                     expand = c(0,0),
                     breaks = seq(500,2500,500)) +
  geom_contour(aes(z=delta_surv),col="black",breaks = c(0.1),lty=2) +
  theme_classic(base_size = 15) +
  theme(legend.title = element_text(face = "bold",size = 12),
        legend.position = "bottom",
        legend.box.spacing = unit(0.1,"cm")) +
  xlab(c("day of water year")) +
  ylab(expression(paste("flow "~"(",m^3,s^-1,")")))

#ggsave("Figure5C.png",height = 5,width = 4.5)

# predict how passage probability changes with 30,000 cfs addition
# same steps as previous two plots: build newdat, predict, add flow, predict
newdat <- expand.grid(cond = mean(fish_red$cond,na.rm=T),
                      release_river_km = median(fish_red$release_river_km),
                      water_year_doy = seq(min(fish_red$water_year_doy),max(fish_red$water_year_doy),7),
                      flow_cms = seq(min(fish_red$flow_cms)-10,max(fish_red$flow_cms)+10,100))

prdat <-  predict.gam(g1,newdat,se.fit = T)
newdat$fit <- prdat$fit
newdat$conf <- prdat$se.fit*2
newdat$lwr <- newdat$fit-newdat$conf
newdat$upr <- newdat$fit+newdat$conf

newdat <- newdat %>%
  mutate(fit = ilink(fit),
         upr = ilink(upr),
         lwr = ilink(lwr),
         conf = (upr-lwr)/2,
         fit = ifelse(conf<fit,fit,NA))

newdat <- mutate(newdat,flow_cms = flow_cms+(30000/35.315))
prdat <-  predict.gam(g1,newdat,se.fit = T)

newdat$fit_inc <- prdat$fit
newdat$conf_inc <- prdat$se.fit*2
newdat$lwr_inc <- newdat$fit_inc-newdat$conf_inc
newdat$upr_inc <- newdat$fit_inc+newdat$conf_inc

newdat <- newdat %>%
  mutate(fit_inc = ilink(fit_inc),
         upr_inc = ilink(upr_inc),
         lwr_inc = ilink(lwr_inc),
         conf_inc = (upr_inc-lwr_inc)/2,
         fit_inc = ifelse(conf_inc<fit_inc,fit_inc,NA))

newdat <- newdat %>%
  mutate(delta_surv = fit_inc-fit,
         flow_cms = flow_cms-(30000/35.315))

# Figure 5E - effect of 30,000 cfs flow addition
newdat %>%
  ggplot(aes(water_year_doy,flow_cms)) +
  geom_raster(aes(fill = delta_surv), interpolate=F) +
  scale_fill_gradientn(colours=c("red","blanchedalmond","blue"),
                       name=expression(bold(atop("passage"~Delta~"(Pr)",paste("+849"~m^3,s^-1)))),
                       limits = c(0,0.4),breaks = seq(0,0.4,0.1),na.value = "grey70") +
  scale_x_continuous(limits = c(min(newdat$water_year_doy)-10,max(newdat$water_year_doy)+10),
                     expand = c(0,0),
                     breaks = seq(min(fish_red$water_year_doy),max(fish_red$water_year_doy),50),
                     sec.axis = dup_axis(labels = c(wy_lab$d),name = "calendar date")) +
  scale_y_continuous(limits = c(min(newdat$flow_cms)-100,max(newdat$flow_cms)+100),
                     expand = c(0,0),
                     breaks = seq(500,2500,500)) +
  geom_contour(aes(z=delta_surv),col="black",breaks = c(0.1),lty=2) +
  theme_classic(base_size = 15) +
  theme(legend.title = element_text(face = "bold",size = 12),
        legend.position = "bottom",
        legend.box.spacing = unit(0.1,"cm")) +
  xlab(c("day of water year")) +
  ylab(expression(paste("flow "~"(",m^3,s^-1,")")))

#ggsave("Figure5E.png",height = 5,width = 4.5)

# repeat the 3 previous visualizations but for temperature GAM
# only difference in building data to predict over is full range of temperature, not day of water year
newdat <- expand.grid(cond = mean(fish_red$cond,na.rm=T),
                      release_river_km = median(fish_red$release_river_km),
                      temp_c = seq(min(fish_red$temp_c),max(fish_red$temp_c)+0.5,0.5),
                      flow_cms = seq(min(fish_red$flow_cms)-10,max(fish_red$flow_cms)+10,100))

prdat <-  predict.gam(g1t,newdat,se.fit = T)
newdat$fit <- prdat$fit
newdat$conf <- prdat$se.fit*2
newdat$lwr <- newdat$fit-newdat$conf
newdat$upr <- newdat$fit+newdat$conf

ilink <- family(g1t)$linkinv
newdat <- newdat %>%
  mutate(fit = ilink(fit),
         upr = ilink(upr),
         lwr = ilink(lwr),
         conf = (upr-lwr)/2,
         fit = ifelse(conf<fit,fit,NA))

# Figure 5B - passage probability with respect to flow and temperature
newdat %>%
  ggplot(aes(temp_c,flow_cms)) +
  geom_raster(aes(fill = fit), interpolate=F) +
  scale_fill_gradientn(colours=c("red","blanchedalmond","blue"),
                       name=expression(bold(atop("Passage", paste("(Pr)")))),
                       limits = c(0,0.4),breaks = seq(0,0.4,0.1),na.value = "grey70") +
  geom_point(aes(temp_c,flow_cms,pch=origin),data=fish_red,alpha=0.2) +
  scale_x_continuous(limits = c(min(newdat$temp_c)-0.5,max(newdat$temp_c)+0.5),
                     expand = c(0,0),
                     breaks = seq(round(min(fish_red$temp_c)),round(max(fish_red$temp_c)),4)) +
  scale_y_continuous(limits = c(min(newdat$flow_cms)-100,max(newdat$flow_cms)+100),
                     expand = c(0,0),
                     breaks = seq(500,2500,500)) +
  geom_contour(aes(z=fit),col="black",breaks = c(0.1),lty=2) +
  theme_classic(base_size = 15) +
  guides(shape=guide_legend(title="origin",reverse = T,keywidth = 0,direction = "vertical",nrow = 1)) +
  theme(legend.title = element_text(face = "bold",size = 12),
        legend.position = "bottom",
        legend.box.spacing = unit(0.1,"cm")) +
  xlab(c("temperature (°C)")) +
  ylab(expression(paste("flow "~"(",m^3,s^-1,")")))

#ggsave("Figure5B.png",height = 5,width = 4.5)

newdat <- mutate(newdat,flow_cms = flow_cms+(10000/35.315))
prdat <-  predict.gam(g1t,newdat,se.fit = T)

newdat$fit_inc <- prdat$fit
newdat$conf_inc <- prdat$se.fit*2
newdat$lwr_inc <- newdat$fit_inc-newdat$conf_inc
newdat$upr_inc <- newdat$fit_inc+newdat$conf_inc

newdat <- newdat %>%
  mutate(fit_inc = ilink(fit_inc),
         upr_inc = ilink(upr_inc),
         lwr_inc = ilink(lwr_inc),
         conf_inc = (upr_inc-lwr_inc)/2,
         fit_inc = ifelse(conf_inc<fit_inc,fit_inc,NA))

newdat <- newdat %>%
  mutate(delta_surv = fit_inc-fit,
         flow_cms = flow_cms-(10000/35.315))

# Figure 5D - effect of 10,000 cfs flow addition
newdat %>%
  ggplot(aes(temp_c,flow_cms)) +
  geom_raster(aes(fill = delta_surv), interpolate=F) +
  scale_fill_gradientn(colours=c("red","blanchedalmond","blue"),
                       name=expression(bold(atop("passage"~Delta~"(Pr)",paste("+283"~m^3,s^-1)))),
                       limits = c(0,0.15),breaks = c(0,0.05,0.1,0.15),na.value = "grey70") +
  scale_x_continuous(limits = c(min(newdat$temp_c)-0.5,max(newdat$temp_c)+0.5),
                     expand = c(0,0),
                     breaks = seq(round(min(fish_red$temp_c)),round(max(fish_red$temp_c)),4)) +
  scale_y_continuous(limits = c(min(newdat$flow_cms)-100,max(newdat$flow_cms)+100),
                     expand = c(0,0),
                     breaks = seq(500,2500,500)) +
  geom_contour(aes(z=delta_surv),col="black",breaks = c(0.1),lty=2) +
  theme_classic(base_size = 15) +
  theme(legend.title = element_text(face = "bold",size = 12),
        legend.position = "bottom",
        legend.box.spacing = unit(0.1,"cm")) +
  xlab(c("temperature (°C)")) +
  ylab(expression(paste("flow "~"(",m^3,s^-1,")")))

#ggsave("Figure5D.png",height = 5,width = 4.5)

newdat <- expand.grid(cond = mean(fish_red$cond,na.rm=T),
                      release_river_km = median(fish_red$release_river_km),
                      temp_c = seq(min(fish_red$temp_c),max(fish_red$temp_c)+0.5,0.5),
                      flow_cms = seq(min(fish_red$flow_cms)-10,max(fish_red$flow_cms)+10,100))

prdat <-  predict.gam(g1t,newdat,se.fit = T)

newdat$fit <- prdat$fit
newdat$conf <- prdat$se.fit*2
newdat$lwr <- newdat$fit-newdat$conf
newdat$upr <- newdat$fit+newdat$conf

newdat <- newdat %>%
  mutate(fit = ilink(fit),
         upr = ilink(upr),
         lwr = ilink(lwr),
         conf = (upr-lwr)/2,
         fit = ifelse(conf<fit,fit,NA))

newdat <- mutate(newdat,flow_cms = flow_cms+(30000/35.315))
prdat <-  predict.gam(g1t,newdat,se.fit = T)

newdat$fit_inc <- prdat$fit
newdat$conf_inc <- prdat$se.fit*2
newdat$lwr_inc <- newdat$fit_inc-newdat$conf_inc
newdat$upr_inc <- newdat$fit_inc+newdat$conf_inc

newdat <- newdat %>%
  mutate(fit_inc = ilink(fit_inc),
         upr_inc = ilink(upr_inc),
         lwr_inc = ilink(lwr_inc),
         conf_inc = (upr_inc-lwr_inc)/2,
         fit_inc = ifelse(conf_inc<fit_inc,fit_inc,NA))

newdat <- newdat %>%
  mutate(delta_surv = fit_inc-fit,
         flow_cms = flow_cms-(30000/35.315))

# Figure 5F - effect of 30,000 cfs flow addition
newdat %>%
  ggplot(aes(temp_c,flow_cms)) +
  geom_raster(aes(fill = delta_surv), interpolate=F) +
  scale_fill_gradientn(colours=c("red","blanchedalmond","blue"),
                       name=expression(bold(atop("passage"~Delta~"(Pr)",paste("+849"~m^3,s^-1)))),
                       limits = c(0,0.4),breaks = seq(0,0.4,0.1),na.value = "grey70") +
  scale_x_continuous(limits = c(min(newdat$temp_c)-0.5,max(newdat$temp_c)+0.5),
                     expand = c(0,0),
                     breaks = seq(round(min(fish_red$temp_c)),round(max(fish_red$temp_c)),4)) +
  scale_y_continuous(limits = c(min(newdat$flow_cms)-100,max(newdat$flow_cms)+100),
                     expand = c(0,0),
                     breaks = seq(500,2500,500)) +
  geom_contour(aes(z=delta_surv),col="black",breaks = c(0.1),lty=2) +
  theme_classic(base_size = 15) +
  theme(legend.title = element_text(face = "bold",size = 12),
        legend.position = "bottom",
        legend.box.spacing = unit(0.1,"cm")) +
  xlab(c("temperature (°C)")) +
  ylab(expression(paste("flow "~"(",m^3,s^-1,")")))

#ggsave("Figure5F.png",height = 5,width = 4.5)
