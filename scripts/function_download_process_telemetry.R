
## In its current form, this function downloads the relevant acoustic telemetry data (for the indicated tributary/season)
## and then identifies which fish 'survived' the trib (for simplicity, were last detected outside the trib)
## We can change this function later to do more complex processing for a more detailed survival model
## Inputs:
##      * trib: name of tributary (current options: 'feather', 'american')
##      * season: in cases where we have data for multiple runs, specify the season (currently, 'fall' or 'spring')
##      * save_dir: directory where outputs should be saved (created if it doesn't exist). 
##      * return_type: 'fishdat' only returns data summary for each fish, "all" returns full detection history
##      * censor_upstream: censor fish that move upstream more than X km; default = 5. (Set to a large number, e.g. 1000, to turn off censoring)
##      * speed_limit: filter observations where the fish appears to move at a speed greater than X km/day (note: very short bursts of speed are allowed, see below)
##      * waterfall_plots: T or F, create and save waterfall plots for the movement paths of each fish (for processing purposes, to spot potential false detections)

## Outputs: a dataframe containing information about each fish released in the trib, and whether or not that fish survived ('survived.trib')

require(rerddap); require(tidyr); require(dplyr); library(geosphere); library(ggplot2); library(this.path)
setwd(this.path::here()); setwd('..')
#trib = 'feather'; season = 'spring'; save_dir = 'data_processed'; waterfall_plots = F; censor_upstream = 5; speed_limit = 120; return_type = "fishdat"; waterfall_plots = F

download_process_telemetry = function(trib, season = NULL, save_dir = "data_processed",return_type = "fishdat",
                                      censor_upstream = 5, speed_limit = 120, waterfall_plots = F){
  
  ## Set up folder to save (if it doesn't exist)
  if(!dir.exists(save_dir)){ dir.create(save_dir) }
  
  
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
    ## Jessie can confirm these are correct - based on the study ids that fall within the coordinates of the Am river
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
  dat_fish_recv <- merge(dat_fish,
                         subset(recvs, select = -c(receiver_serial_number, receiver_river_km, receiver_general_river_km, receiver_location, receiver_general_location)),
                         by = "dep_id")
  
  print('done downloading and merging data')
  
  ## Reorder data by fish and time
  dat_fish_recv = dat_fish_recv[order(dat_fish_recv$fish_id, dat_fish_recv$first_time), ]
  
  ## remove big objects from workspace
  rm(dat); rm(dat_fish)
  
  ## Fish in studies of interest
  fish_studied = subset(fish, study_id %in% studyids)
  
  ## Receiver summary
  recvs_studied = recvs %>% filter(dep_id %in% unique(dat_fish_recv$dep_id)) %>% 
                            group_by(receiver_general_location) %>%
                            arrange(receiver_general_river_km, receiver_location) %>%
                            select(dep_id, receiver_general_location, receiver_location, receiver_region,
                                   receiver_general_river_km, receiver_river_km, 
                                   receiver_general_latitude, receiver_general_longitude,
                                   latitude, longitude, receiver_start, receiver_end) %>%
                           ungroup()
  
  recvs_location_summary = recvs_studied %>% filter(!duplicated(receiver_location)) %>%
                                             select(-dep_id, -receiver_start, -receiver_end)
  recvs_general_location_summary = recvs_location_summary %>% filter(!duplicated(receiver_general_location)) %>%
                                         select(-receiver_river_km, -latitude, -longitude, -receiver_location)
  
  ### FILTERING ###
  
  colnames_copy = colnames(dat_fish_recv)
  
  ## Remove observations that:
  
  # Are a day or more before release, or more than 60 days after release
  dat_fish_recv$release_tDiff = as.numeric(difftime(dat_fish_recv$first_time, dat_fish_recv$fish_release_date, units = 'days'))
  dat_fish_recv = subset(dat_fish_recv, release_tDiff <= 60 & release_tDiff >= -1)
  
  # Follow a gap in observations greater than two weeks
  filter_gaps <- dat_fish_recv %>% arrange(fish_id, first_time) %>% 
                                   group_by(fish_id) %>%
                                   mutate(obs_gap = as.numeric(difftime(first_time, lag(first_time, n = 1), units = 'days'))) %>%
                                   filter(obs_gap > 14) %>%
                                   mutate(censor_time1 = min(first_time)) %>%
                                   select(fish_id, censor_time1)
  dat_fish_recv = left_join(dat_fish_recv, filter_gaps, by = 'fish_id')
  dat_fish_recv = subset(dat_fish_recv, is.na(censor_time1) | first_time < censor_time1)
  
  ## Are biologically unrealistic (more than x distance in y time, either by lat/long or river km)
  dat_fish_recv$row_i = 1:nrow(dat_fish_recv)
  
  #initial step
  filter_speeds = dat_fish_recv %>% arrange(fish_id, first_time) %>% 
                                    group_by(fish_id) %>%
                                    mutate(obs_gap = as.numeric(difftime(first_time, lag(first_time, n = 1), units = 'days')),
                                           obs_gap2 = pmax(obs_gap, 1/24), #set minimum gap size of 1 hour (otherwise, a fish detected at receivers 100m apart 10 s apart appears to be going a speed of 864 km/day)
                                           speed_latlon = distHaversine(cbind(receiver_general_longitude, receiver_general_latitude), cbind(lag(receiver_general_longitude, n = 1), lag(receiver_general_latitude, n=1)))/(1000*obs_gap2),
                                           speed_river_km = abs(receiver_river_km - lag(receiver_river_km, n = 1))/obs_gap2,
                                           filter_point = speed_latlon > speed_limit | speed_river_km > speed_limit)
  filter_speeds_rows = filter_speeds %>% filter(filter_point == T) %>%
                                         group_by(fish_id) %>%
                                         summarize(filter_count = n(), row_filter = min(row_i)) 
  
  #remove first point flagged for each fish, then repeat process 
  filter_speeds_rows_new = filter_speeds_rows
  dat_fish_recv_filt = dat_fish_recv 
  done_filtering = F
  
  while(done_filtering == F){
    dat_fish_recv_filt = dat_fish_recv_filt %>%  filter(fish_id %in% unique(filter_speeds_rows_new$fish_id),
                                                        !(row_i %in% unique(filter_speeds_rows$row_filter) ) )
    filter_speeds_new  = dat_fish_recv_filt %>% arrange(fish_id, first_time) %>% 
                                           group_by(fish_id) %>%
                        mutate(obs_gap = as.numeric(difftime(first_time, lag(first_time, n = 1), units = 'days')),
                               obs_gap2 = pmax(obs_gap, 1/24), #set minimum gap size of 1 hour (otherwise, a fish detected at receivers 100m apart 10 s apart appears to be going a speed of 864 km/day)
                               speed_latlon = distHaversine(cbind(receiver_general_longitude, receiver_general_latitude), cbind(lag(receiver_general_longitude, n = 1), lag(receiver_general_latitude, n=1)))/(1000*obs_gap2),
                               speed_river_km = abs(receiver_river_km - lag(receiver_river_km, n = 1))/obs_gap2,
                               filter_point = speed_latlon > speed_limit | speed_river_km > speed_limit)
    
    if(sum(filter_speeds_new$filter_point, na.rm = T) > 0){
      filter_speeds_rows_new = filter_speeds_new %>% filter(filter_point == T) %>%
        group_by(fish_id) %>%
        summarize(filter_count = n(), row_filter = min(row_i))
      
      filter_speeds_rows = rbind(filter_speeds_rows, filter_speeds_rows_new)
      #print(nrow(filter_speeds_rows_new)) #fish remaining counter
      
    }else{
      done_filtering = T
    }
    
  }
  
  #remove all flagged rows from dat_fish_recv
  dat_fish_recv = dat_fish_recv %>% filter(!(row_i %in% unique(filter_speeds_rows$row_filter) ))

  
  
  ## Censor fish that swim upstream (more than X km)
  dat_fish_recv = dat_fish_recv %>%
    arrange(fish_id, first_time) %>% 
    group_by(fish_id) %>%
    mutate(receiver_general_river_km_next1 = dplyr::lead(receiver_general_river_km, n = 1),
           receiver_general_location_next1 = dplyr::lead(receiver_general_location, n = 1),
           receiver_location_next1 = dplyr::lead(receiver_location, n = 1),
           too_far_upstream = receiver_general_river_km_next1 - receiver_general_river_km > censor_upstream) %>%
    ungroup()
  
  #view censored fish/locations
  upstream_censored_fish = subset(dat_fish_recv, too_far_upstream == T)[, c("fish_id", 'first_time', 'receiver_general_location', 'receiver_general_location_next1',
                                                  'receiver_general_river_km','receiver_general_river_km_next1',
                                                  "receiver_location", "receiver_location_next1", "too_far_upstream")]
  
  #censor fish when they swim upstream
  dat_fish_recv =  dat_fish_recv %>% 
                   arrange(fish_id, first_time) %>% 
                   group_by(fish_id) %>%
                   mutate(too_far_upstream_lag1 = dplyr::lag(too_far_upstream, n = 1)) %>%
                   filter(cumsum(replace_na(too_far_upstream_lag1, F)) < 1) %>%
                   ungroup()
  
  ## remove extra columns
  dat_fish_recv = dat_fish_recv[, colnames_copy]
  
  
  print('done filtering data')
  
  
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
  
  
  ### PLOTTING ###
  
  if(waterfall_plots == T){
    for(study_i in unique(fish_studied$study_id)){
      fish_study_i = subset(fish_studied, study_id == study_i)
      detects_study_i = subset(dat_fish_recv, study_id == study_i)
      detects_study_i <- detects_study_i %>%
        arrange(fish_id, first_time) %>% 
        group_by(fish_id) %>%
        mutate(receiver_general_river_km_lag1 = dplyr::lag(receiver_general_river_km, n = 1),
               first_time_lag1 = dplyr::lag(first_time, n = 1)) %>%
        ungroup()
      detects_study_i$upstream = ifelse(detects_study_i$receiver_general_river_km > detects_study_i$receiver_general_river_km_lag1, 1, 0)
      plot_i = ggplot(detects_study_i) + 
               geom_step(aes(x = first_time, y = receiver_general_river_km, group = fish_id), alpha = .5) +
               theme_classic() + 
               labs(x = "Date", y = "Receiver Location (river km)", title = paste0("Fish Detections (", study_i, ")")) 
      ggsave(plot = plot_i, filename = paste0(save_dir, "/waterfall_plots_", study_i, ".png"), width = 6, height = 6)         
    }
  }
  
  
  
  if(return_type == 'fishdat'){
    return(fish_studied_final)
  }else if(return_type == 'all'){
    return(list(dat_fish_recv = dat_fish_recv, fishdat = fish_studied_final, recvdat = recvs_studied,
                recvs_location_summary = recvs_location_summary, recvs_general_summary = recvs_general_location_summary))
  }

  
  
}

