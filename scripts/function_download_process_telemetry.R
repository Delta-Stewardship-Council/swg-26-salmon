
## In its current form, this function downloads the relevant acoustic telemetry data (for the indicated tributary/season)
## and then identifies which fish 'survived' the trib (for simplicity, were last detected outside the trib)
## We can change this function later to do more complex processing for a more detailed survival model

## Outputs: a dataframe containing information about each fish released in the trib, and whether or not that fish survived ('survived.trib')

require(rerddap); require(tidyr); require(dplyr)
download_process_telemetry = function(trib, season = NULL){
  
  
  #### IDENTIFY RELEVANT STUDIES ####
  
  #unique_studies <- unique(fish[,c("study_id","email")])
  
  ## Identify relevant studies for each trib/season (add more here as we decide what tribs to analyze)
  if(trib == 'feather'){
    trib_receiver_regions = c("Feather_R", "Yuba R") #if a fish is detected outside these regions, they are considered to have survived the trib
    if(season == 'spring'){
      studyids = c("FR_Spring_2013", "FR_Spring_2014", "FR_Spring_2015",
                         "FR_Spring_2019", "FR_Spring_2020", "FR_Spring_2021", 
                         "FR_Spring_2023", "FR_Spring_2024", "FR_Spring_2025")
    }else if(season == 'fall'){
      studyids = c("FR_Fall_2012", "FR_Fall_Chinook_2023", "FRH_Fall_2021", "FRW_Fall_2021")
    }else{
      stop('specify season for feather river')
    }
  }else if(trib == 'american'){
    ## Anna--haven't double checked these yet
    studyids = c("Nimbus_Fall_2016", "Nimbus_Fall_2017", "Nimbus_Fall_2018", "Nimbus_Fall_2022",
                       "Nimbus_Fall_2023", "Nimbus_Fall_2024", "Nimbus_Fall_2025")
    trib_receiver_regions = c("Amer R")
  }else{
    stop('must identify relevant studies for trib')
  }
  
  
  #### DOWNLOAD ALL TAGGED FISH DATA ####
  fish <- tabledap('FED_JSATS_taggedfish', url = "https://oceanview.pfeg.noaa.gov/erddap/")
  fish$fish_release_date <- as.POSIXct(fish$fish_release_date, "%m/%d/%Y %H:%M:%S", tz = "Etc/GMT+8")
  fish$fish_id_prefix <- substr(fish$fish_id,start = 1, stop = (nchar(fish$fish_id)-4))
  
  #### DOWNLOAD ALL RECEIVER DEPLOYMENT DATA ####
  recvs <- tabledap('FED_JSATS_receivers', url = "https://oceanview.pfeg.noaa.gov/erddap/")
  
  #### DOWNLOAD DETECTION DATA FOR RELEVANT STUDIES ####
  datalist <- list()
  for (i in studyids){
    constraint <-  noquote(paste0("'study_id=\"",i,"\"'"))
    datalist[[i]] <- tabledap('FED_JSATS_detects', url = "https://oceanview.pfeg.noaa.gov/erddap/", str2lang(constraint))
  }
  dat <- do.call(rbind,datalist) 
  rm(datalist) #remove big object
  
  
  #### MAKE DATA NICE ####
  
  ## format time
  dat$first_time <- as.POSIXct(dat$first_time, origin = '1970-01-01', format = "%Y-%m-%d %H:%M:%OS", tz = "Etc/GMT+8")
  dat$last_time <- as.POSIXct(dat$last_time, origin = '1970-01-01', format = "%Y-%m-%d %H:%M:%OS", tz = "Etc/GMT+8")
  
  ## Associate detection data to tagging and receiver data 
  dat_fish <- merge(dat, fish, by = c("study_id", "fish_id"))
  dat_fish_recv <- merge(dat_fish, recvs, by = "dep_id")
  
  ## remove big objects from workspace
  rm(dat); rm(dat_fish)
  
  ## Fish in studies of interest
  fish_studied = subset(fish, study_id %in% studyids)
  
  ### FIND WHAT FISH SURVIVED THE TRIB ##
  # easy way: they 'survived' if they were last detected outside the trib
  # not perfect, we should do further processing here
  
  ## last detection of each fish
  detect_max <- dat_fish_recv %>% group_by(fish_id) %>% slice(which.max(last_time))
  
  ## check if last detected location was in the trib or not
  detect_max$survived.trib = ifelse(detect_max$receiver_region %in% trib_receiver_regions, 0, 1)
  fish_studied = left_join(fish_studied, detect_max[, c("fish_id", 'survived.trib')], by = 'fish_id')
  
  ## do additional filtering
  
  ## rearrange columns for convenience
  fish_studied_final = fish_studied %>% select(fish_id, study_id, survived.trib, fish_release_date, release_location, release_latitude, release_longitude, fish_length, fish_weight, everything())
  
  return(fish_studied_final)
  
  
}
