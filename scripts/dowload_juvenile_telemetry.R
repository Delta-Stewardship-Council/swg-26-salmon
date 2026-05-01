
## This is a modified version of the 20231220_errdap_download_script.R script for our data processing purposes

## additional details: https://oceanview.pfeg.noaa.gov/CalFishTrack/

## load packages
library(rerddap)
library(ggplot2)
library(mapdata)
library(ggrepel)
library(tidyr); library(dplyr)

## Clear cache to get newest data
cache_delete_all()

#### DOWNLOAD ALL TAGGED FISH DATA ####
fish <- tabledap('FED_JSATS_taggedfish', url = "https://oceanview.pfeg.noaa.gov/erddap/")
## Format release datetime field
fish$fish_release_date <- as.POSIXct(fish$fish_release_date, "%m/%d/%Y %H:%M:%S", tz = "Etc/GMT+8")
## Create a column for fish_id prefixes, which are consistent throughout a study and can be used for quick data queries
fish$fish_id_prefix <- substr(fish$fish_id,start = 1, stop = (nchar(fish$fish_id)-4))

#### DOWNLOAD ALL RECEIVER DEPLOYMENT DATA ####
recvs <- tabledap('FED_JSATS_receivers', url = "https://oceanview.pfeg.noaa.gov/erddap/")


#### IDENTIFY RELEVANT STUDIES ####
## This will tell you unique study_id names (unifying name for all fish released in a year for a study)
unique_studies <- unique(fish[,c("study_id","email")])

## Find relevant studies
feather_river_spring_studies = c("FR_Spring_2013", "FR_Spring_2014", "FR_Spring_2015",
                                 "FR_Spring_2019", "FR_Spring_2020", "FR_Spring_2021", 
                                 "FR_Spring_2023", "FR_Spring_2024", "FR_Spring_2025",
                                 "FR_Spring_Delta_2013", "FR_Spring_Delta_2014", "FR_Spring_Delta_2015")
#feather_river_spring_studies = subset(unique_studies, startsWith(study_id, "FR_Spring"))$study_id

feather_river_spring_studies_select = c("FR_Spring_2013", "FR_Spring_2014", "FR_Spring_2015",
                                 "FR_Spring_2019", "FR_Spring_2020", "FR_Spring_2021", 
                                 "FR_Spring_2023", "FR_Spring_2024", "FR_Spring_2025")


feather_river_fall_studies = c("FR_Fall_2012", "FR_Fall_Chinook_2023", "FRH_Fall_2021", "FRW_Fall_2021")
#feather_river_studies = subset(unique_studies, startsWith(study_id, "FR"))$study_id
#feather_river_spring_studies = subset(unique_studies, startsWith(study_id, "FR_Spring"))$study_id
#feather_river_fall_studies = feather_river_studies[!(feather_river_studies %in% feather_river_spring_studies)]

nimbus_studies = c("Nimbus_Fall_2016", "Nimbus_Fall_2017", "Nimbus_Fall_2018", "Nimbus_Fall_2022",
                   "Nimbus_Fall_2023", "Nimbus_Fall_2024", "Nimbus_Fall_2025")
#nimbus_studies = subset(unique_studies, startsWith(study_id, "Nimbus"))$study_id

#yuba_studies = c("Lower_Yuba_FRH_Chinook_2021", "Lower_Yuba_FRH_Chinook_2022", "Lower_Yuba_FRH_Chinook_2023", "Lower_Yuba_FRH_Chinook_2024")


## Decide which studies to pull full detection data for
## (the rest of the script focuses on the studies in this vector)
studyids <- c(feather_river_spring_studies_select) 


#### DOWNLOAD DETECTION DATA FOR RELEVANT STUDIES ####

## Find out some general details on the detections database
db <- info('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/")

## Download data from 2 or more studyIDs (takes a minute)
datalist <- list() #Make a blank list for saving detection data in the loop below
for (i in studyids){
  constraint <-  noquote(paste0("'study_id=\"",i,"\"'"))
  datalist[[i]] <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", str2lang(constraint))
}
dat <- do.call(rbind,datalist) #append the lists into one dataframe
rm(datalist) #remove big object


#### MAKE DATA NICE ####

## format time
dat$first_time <- as.POSIXct(dat$first_time, origin = '1970-01-01', format = "%Y-%m-%d %H:%M:%OS", tz = "Etc/GMT+8")
dat$last_time <- as.POSIXct(dat$last_time, origin = '1970-01-01', format = "%Y-%m-%d %H:%M:%OS", tz = "Etc/GMT+8")

## Associate detection data to tagging data 
dat_fish <- merge(dat, fish, by = c("study_id", "fish_id"))

## Associate detection data to receiver data
dat_fish_recv <- merge(dat_fish, recvs, by = "dep_id")

## remove big objects from workspace
rm(dat); rm(dat_fish)

## make sure formatting is correct
dat_fish_recv$latitude <- as.numeric(dat_fish_recv$latitude)
dat_fish_recv$longitude <- as.numeric(dat_fish_recv$longitude)

## For convenience, pull out fish data specific to studies of interest
fish_of_interest = subset(fish, study_id %in% studyids)

## Find unique release locations
unique_release_locs = subset(fish_of_interest, !duplicated(release_location))[, c("release_location", "release_longitude", "release_latitude", "release_river_km")]

#### Find first, last, and count of detections per fish per general location ####
#detect_minmaxcount <- aggregate(data=dat, first_time~study_id+ fish_id+ receiver_general_location, FUN=min)
#detect_minmaxcount <- merge(detect_minmaxcount, aggregate(data=dat, last_time~study_id+ fish_id+ receiver_general_location, FUN=max))
#detect_minmaxcount <- merge(detect_minmaxcount, aggregate(data=dat, detection_count~study_id+ fish_id+ receiver_general_location, FUN=sum))
#detect_minmaxcount <- detect_minmaxcount[order(detect_minmaxcount$fish_id, detect_minmaxcount$first_time),]

### Fun initial analyses: for spring run feather river releases, classify fish by if they were detected out of the feather river or not
out_feather_dat = subset(dat_fish_recv, !(receiver_region %in% c("Feather_R", "Yuba R")))
fish_of_interest$survived.feather = ifelse(fish_of_interest$fish_id %in% unique(out_feather_dat$fish_id), 1, 0)
survived_by_study = fish_of_interest %>% group_by(study_id) %>% summarize(feather.survival.prob = mean(survived.feather))
ggplot(survived_by_study) + 
  geom_col(aes(y = study_id, x = feather.survival.prob)) +
  theme_classic()


#### Make a map of detection locations ####

## summarize across all studies by unique fish visits per receiver general location
detect_summary <- aggregate(list(fish_count = dat_fish_recv$fish_id), by = list(receiver_general_location = dat_fish_recv$receiver_general_location.x, latitude = dat_fish_recv$receiver_general_latitude, longitude = dat_fish_recv$receiver_general_longitude), function(x){length(unique(x))})

## for each study, summarize data by unique fish visits per receiver general location
detect_summary_study <- aggregate(list(fish_count = dat_fish_recv$fish_id), by = list(study_id = dat_fish_recv$study_id, receiver_general_location = dat_fish_recv$receiver_general_location.x, latitude = dat_fish_recv$receiver_general_latitude, longitude = dat_fish_recv$receiver_general_longitude), function(x){length(unique(x))})


## Set boundary box for map
xlim <- c(-123, -121)
ylim <- c(37.5, 40.6)
usa <- map_data("worldHires", ylim = ylim, xlim = xlim)
rivers <- map_data("rivers", ylim = ylim, xlim = xlim)
rivers <- rivers[rivers$lat < max(ylim),]

## all receivers
ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = "grey80") +
  geom_path(data = rivers, aes(x = long, y = lat, group = group), size = 1, color = "white", lineend = "round") +
  geom_point(data = detect_summary, aes(x = longitude, y = latitude), shape=23, fill="blue", color="darkred", size=1) +
  #geom_text_repel(data = subset(detect_summary, fish_count > 1000), aes(x = longitude, y = latitude, label = fish_count)) +
  theme_bw() + ylab("latitude") + xlab("longitude") +
  coord_fixed(1.3, xlim = xlim, ylim = ylim) +
  ggtitle("Location of study detections w/ count of unique fish visits")


## feather river
ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = "grey80") +
  geom_path(data = rivers, aes(x = long, y = lat, group = group), size = 1, color = "white", lineend = "round") +
  geom_point(data = subset(recvs, receiver_region %in% c("Feather_R", "Yuba R")), aes(x = longitude, y = latitude), shape=23, fill="blue", color="darkred", size=1) +
  geom_point(data = subset(recvs, !receiver_region %in% c("Feather_R", "Yuba R")), aes(x = longitude, y = latitude), shape=23, fill="black", color="black", size=1) +
  theme_bw() + ylab("latitude") + xlab("longitude") +
  coord_fixed(1.3, xlim = c(-122, -121), ylim = c(38.5, 39.7)) +
  ggtitle("Location of study detections w/ count of unique fish visits")


## facet_wrap by study
#detect_summmary_study1 = subset(detect_summary_study,   study_id %in% feather_river_spring_studies[1:4])
#detect_summmary_study2 = subset(detect_summary_study,  study_id %in% feather_river_spring_studies[5:8])
#detect_summmary_study3 = subset(detect_summary_study, study_id %in% feather_river_spring_studies[9:12])
#ggplot() +
#  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = "grey80") +
#  geom_path(data = rivers, aes(x = long, y = lat, group = group), size = 1, color = "white", lineend = "round") +
#  geom_point(data = detect_summmary_study3, aes(x = longitude, y = latitude), shape=23, fill="blue", color="darkred", size=1) +
#  #geom_text_repel(data = detect_summmary_study1, aes(x = longitude, y = latitude, label = fish_count), max.overlaps = 30) +
#  facet_wrap(~study_id, nrow = 1)+
#  theme_bw() + ylab("latitude") + xlab("longitude") +
#  coord_fixed(1.3, xlim = xlim, ylim = ylim) +
#  ggtitle("Location of study detections w/ count of unique fish visits")

## unique release locations
ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = "grey80") +
  geom_path(data = rivers, aes(x = long, y = lat, group = group), size = 1, color = "white", lineend = "round") +
  geom_point(data = unique_release_locs, aes(x = release_longitude, y = release_latitude), shape=23, fill="blue", color="darkred", size=2) +
  geom_text_repel(data = unique_release_locs, aes(x = release_longitude, y = release_latitude, label = release_location), box.padding = 1, nudge_x = 4) +
  theme_bw() + ylab("latitude") + xlab("longitude") +
  coord_fixed(1.3, xlim = xlim, ylim = ylim) +
  ggtitle("Locations of releases")


## for each study, plot release locations (for my own visualization purposes)
#for(i in 1:nrow(unique_studies)){
#  fish_study = subset(fish, study_id == unique_studies$study_id[i])
#  releases_study = subset(fish_study, !duplicated(release_location))[, c("release_location", "release_longitude", "release_latitude", "release_river_km")]
#  print(ggplot() +
#    geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = "grey80") +
#    geom_path(data = rivers, aes(x = long, y = lat, group = group), size = 1, color = "white", lineend = "round") +
#    geom_point(data = releases_study, aes(x = release_longitude, y = release_latitude), shape=23, fill="blue", color="darkred", size=2) +
#    geom_text_repel(data = releases_study, aes(x = release_longitude, y = release_latitude, label = release_location), box.padding = 1, nudge_x = 4) +
#    theme_bw() + ylab("latitude") + xlab("longitude") +
#    coord_fixed(1.3, xlim = xlim, ylim = ylim) +
#    ggtitle(paste0(unique_studies$study_id[i], "; release locations")))
#}

#### DESCRIPTION: TAGGED FISH DATA ####
## First, a quick tutorial on important fields in the "FED_JSATS_taggedfish" dataset:
#"fish_id" (Unique fish identification number, different from tag_id_hex, which can get reused by manufacturer)
#"tag_id_hex" (TagID code, in hexadecimal format)
#"study_id" (unifying name for all fish released in a year for a study)
#"fish_type" (origin of fish, describes species, population, and source. E.g., "CNFH Fall Chinook" represents Coleman National Fish Hatchery Fall-run Chinook salmon)
#"fish_release_date" (release date/time of fish in Pacific Standard Time, i.e., does not adjust for daylight savings)
#"release_location" (Release location name)
#"release_latitude" (Release location latitude, decimal degrees)
#"release_longitude" (Release location longitude, decimal degrees)
#"release_river_km" (Release River Kilometer - distance from Golden Gate by river, km)
#"email" (The PIs email address)

#### DESCRIPTION: DETECTION DATA ####
## "dep_id" a unique identifier for deployment record which can be the join column to merge in additional receiver metadata from the receiver deployment table ("recvs" in this script)
## "fish_id" (Unique fish identification number, different from tag_id_hex, which can get reused by manufacturer). Can be used as join column to merge detection data with fish metadata ("fish" in this script)
## "first_time" (First detection Date/time for that individual tag code at that individual receiver. In Pacific Standard Time, i.e., does not adjust for daylight savings.)
## "last_time" (Last detection Date/time for that individual tag code at that individual receiver. In Pacific Standard Time, i.e., does not adjust for daylight savings.)
## "time" (Detection Date/time in UTC. Note that this is only the first detection time of individual tag code at that individual receiver)
## "detection_count" (Number of detections by that individual tag code at that individual receiver covered by the first and last time columns)
## "receiver_serial_number" (The serial number of the acoustic receiver that detected the tag)
## "receiver_river_km" (Receiver River Kilometer of the individual receivers- distance from Golden Gate, km)
## "study_id" (unifying name for all fish released in a year for a study)
## "receiver_general_river_km" (Receiver River Kilometer of the group of receivers- distance from Golden Gate, km)
## "receiver_location" (Individual Receiver Location name)
## "receiver_general_location" (a unifying name describing 1 location covered by 1 or more individual receivers)
