
## Prepared by Cyril Michel on 2019-07-16; cyril.michel@noaa.gov
## Updated on 2020-06-24
## Updated on 2023-12-14 by BJA

###########################################################################
#### HOW TO PULL IN AUTONOMOUS RECEIVER FISH DETECTION DATA INTO R ########
###########################################################################

## install and load the 'rerddap' library
library(rerddap)

## It is important to delete your cache if you want the newest data. 
## If not, when you rerun the same data queries as before, the command will likely return the old cached data and not the newest data
## If data is unlikely to change or be amended, not deleting your cache will speed up some data queries
cache_delete_all()

#### FIRST LETS DOWNLOAD TAGGED FISH DATA ####
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

## download the tagged fish dataset
fish <- tabledap('FED_JSATS_taggedfish', url = "https://oceanview.pfeg.noaa.gov/erddap/")
## Format release datetime field
fish$fish_release_date <- as.POSIXct(fish$fish_release_date, "%m/%d/%Y %H:%M:%S", tz = "Etc/GMT+8")
## Create a column for fish_id prefixes, which are consistent throughout a study and can be used for quick data queries
fish$fish_id_prefix <- substr(fish$fish_id,start = 1, stop = (nchar(fish$fish_id)-4))

#### NEXT LETS DOWNLOAD RECEIVER DEPLOYMENT DATA ####
recvs <- tabledap('FED_JSATS_receivers', url = "https://oceanview.pfeg.noaa.gov/erddap/")

#### NEXT LETS DOWNLOAD DETECTION DATA ####
## Find out some general details on the detections database
db <- info('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/")
## This will tell you columns and their data types in database
vars <- db$variables
## This will tell you unique study_id names (unifying name for all fish released in a year for a study)
unique_studies <- unique(fish[,c("study_id","email")])

## This will tell you unique receiver_general_location (a unifying name describing 1 location covered by 1 or more individual receivers)
unique_genlocs <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", fields = c("receiver_general_location"), distinct = T)


#*********************************************************************************************
#********************************** IMPORTANT ************************************************
#*********************************************************************************************
# Downloading the entire detection dataset exceeds the 2GB limit per download.
# We recommend downloading only the detection data you need, or at least, downloading the detection data in batches
# See following code snippets for help
#*********************************************************************************************

## Download only data from 1 studyID, here for example, Winter-run 2018 study, with only "important fields". If all fields are desired, remove the "fields" command
dat <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", 'study_id="Winter_H_2018"')

## here are the column definitions for the detection table
##PLEASE NOTE THAT 'local_time' has been replaced with 3 new columns- 'first_time', 'last_time', and 'detection_count'

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





## Download data from 2 or more studyIDs, here for example, Red Bluff diversion dam tagged fish in 2017 and 2018, with only "important fields". If all fields are desired, remove the "fields" command
datalist <- list() #Make a blank list for saving detection data in the loop below
studyids <- c("RBDD_2017", "RBDD_2018") #Make a vector of the studyID names you are interested in
for (i in studyids){
  constraint <-  noquote(paste0("'study_id=\"",i,"\"'"))
  datalist[[i]] <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", str2lang(constraint))
}
dat <- do.call(rbind,datalist) #append the lists into one dataframe

## Download only data from 1 general receiver location. Beware: depending on the general receiver location, this could return a lot of data.
dat <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", 'receiver_general_location="GoldenGateW"')

## Download data from a combination of conditions. For example, study_id="ColemanLateFall_2018" and receiver_general_location="GoldenGateW"
dat <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", 'study_id="ColemanLateFall_2018"', 'receiver_general_location="GoldenGateW"')

## Download only data from a specific time range (in UTC time)
dat <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", 'time>=2019-01-01', 'time<=2019-01-02')

## Finally, download a summary of unique records. Say for example you want to know the unique Fish_ID codes detected somewhere in the array from a studyID. 
## This would include fish released but never detected, as they get assigned 1 row of detection data with their release location as a detection location
unique_fish <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", 'study_id="ColemanLateFall_2018"', fields = c("fish_id"), distinct = T)

## Or, number of unique fish detected at each general receiver location for a studyID. The "distinct = T" command allows us to summarize for only unique records
unique_fish_v_recvs <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", 'study_id="ColemanLateFall_2018"', fields = c("receiver_general_location","fish_id"), distinct = T)

## PLEASE NOTE: IF A DATA REQUEST ABOVE RETURNS SIMPLY "Error: ", THIS LIKELY MEANS THE DATA REQUEST CAME UP WITH ZERO RETURNS or FIELDS (column) NAMES DO NOT MATCH


#_________________________________________________________________________________________________________


#### The following code snippets can help with simple data manipulations, analyses, and visualizations once you've imported your data in R ####

## First, lets format first and last time so R reads it as a Posixct time
dat$first_time <- as.POSIXct(dat$first_time, origin = '1970-01-01', format = "%Y-%m-%d %H:%M:%OS", tz = "Etc/GMT+8")
dat$last_time <- as.POSIXct(dat$last_time, origin = '1970-01-01', format = "%Y-%m-%d %H:%M:%OS", tz = "Etc/GMT+8")

#### Associate detection data to tagging data ####
dat_fish <- merge(dat, fish, by = c("study_id", "fish_id"))

#### Associate detection data to receiver data
dat_fish_recv <- merge(dat_fish, recvs, by = "dep_id")

#### Find first, last, and count of detections per fish per general location ####
detect_minmaxcount <- aggregate(data=dat, first_time~study_id+ fish_id+ receiver_general_location, FUN=min)
detect_minmaxcount <- merge(detect_minmaxcount, aggregate(data=dat, last_time~study_id+ fish_id+ receiver_general_location, FUN=max))
detect_minmaxcount <- merge(detect_minmaxcount, aggregate(data=dat, detection_count~study_id+ fish_id+ receiver_general_location, FUN=sum))
detect_minmaxcount <- detect_minmaxcount[order(detect_minmaxcount$fish_id, detect_minmaxcount$first_time),]

#### Find number of fish released and detected at Benicia per Study ID. ####
## NOTE, if a released fish was never detected again, it will exist in the detection data. 
## These fish get assigned 1 row of detection data with their release location as a detection location
released <- nrow(tabledap('FED_JSATS_detects', url = "http://oceanview.pfeg.noaa.gov/erddap/", 'study_id="Winter_H_2019"', fields = c("fish_id"), distinct = T))
benicia <- nrow(tabledap('FED_JSATS_detects', url = "http://oceanview.pfeg.noaa.gov/erddap/", 'study_id="Winter_H_2019"', 'receiver_general_location="BeniciaW"', fields = c("fish_id"), distinct = T))
## Percent detected at Benicia:
round(benicia/released*100, 2)

#### Make a map of detection locations ####
## first format some fields as necessary
dat_fish_recv$latitude <- as.numeric(dat_fish_recv$latitude)
dat_fish_recv$longitude <- as.numeric(dat_fish_recv$longitude)
## summarize data by unique fish visits per receiver general location
detect_summary <- aggregate(list(fish_count = dat_fish_recv$fish_id), by = list(receiver_general_location = dat_fish_recv$receiver_general_location.x, latitude = dat_fish_recv$receiver_general_latitude, longitude = dat_fish_recv$receiver_general_longitude), function(x){length(unique(x))})
detect_summary$latitude <- as.numeric(detect_summary$latitude)
detect_summary$longitude <- as.numeric(detect_summary$longitude)
library(ggplot2)
library(mapdata)
library(ggrepel)

## Set boundary box for map
xlim <- c(-123, -121)
ylim <- c(37.5, 40.6)
usa <- map_data("worldHires", ylim = ylim, xlim = xlim)
rivers <- map_data("rivers", ylim = ylim, xlim = xlim)
rivers <- rivers[rivers$lat < max(ylim),]
ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = "grey80") +
  geom_path(data = rivers, aes(x = long, y = lat, group = group), size = 1, color = "white", lineend = "round") +
  geom_point(data = detect_summary, aes(x = longitude, y = latitude), shape=23, fill="blue", color="darkred", size=3) +
  geom_text_repel(data = detect_summary, aes(x = longitude, y = latitude, label = fish_count)) +
  theme_bw() + ylab("latitude") + xlab("longitude") +
  coord_fixed(1.3, xlim = xlim, ylim = ylim) +
  ggtitle("Location of study detections w/ count of unique fish visits")
