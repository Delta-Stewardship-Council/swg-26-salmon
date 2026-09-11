## This function takes the dat_fish_recv dataframe (downloaded and processed with the
## download_process_telemetry function) and processes it further. Detections for each fish
## are aggregated by reach, and unknown travel times are sampled.

library(dplyr)

aggregate_reaches = function(dat_fish_recv, reaches_list, release_sites_list){
  
  ## remove unneeded detections
  data = dat_fish_recv %>% filter(receiver_general_location %in% unique(unlist(reaches_list))) %>%
                           arrange(fish_id, first_time)
  
  ## data frame of receiver location names defining each reach
  reaches_df = data.frame(receiver_general_location = unlist(reaches_list),
                          rec_group_id = rep(seq_along(reaches_list), lengths(reaches_list)))
  
  data = left_join(data, reaches_df, by = 'receiver_general_location')
  
  ## summary data frame of fish
  fish_df = data[!duplicated(data$fish_id), c("fish_id", 'study_id', 'fish_type', 'fish_origin', 'fish_date_tagged',
                                              'fish_release_date','release_location') ]
  # release group
  fish_df$release_group = paste0(fish_df$study_id, "_", fish_df$release_location, "_", format(fish_df$fish_release_date, "%j"))
  
  # find which reach each fish starts in
  fish_df$first_reach = sapply(fish_df$release_location, function(val) {
    which(sapply(reaches_list, function(x) val %in% x))
  })
  
  ## For each fish, find first detection in each reach
  data_agg = data %>% group_by(fish_id, rec_group_id) %>%
                             summarize(first_detect = min(first_time)) %>%
                             ungroup() %>% complete(fish_id, rec_group_id)
  data_agg = left_join(data_agg, fish_df[, c("fish_id", "first_reach", 'fish_release_date'), by = 'first_reach'])
  data_agg = subset(data_agg, rec_group_id >= first_reach) #ignore earlier reaches for fish released in later reaches
  
  # fill in first detection time for fish not detected at release site
  data_agg = data_agg %>% mutate(first_detect_processed = if_else(is.na(first_detect) & rec_group_id == first_reach,
                                                          fish_release_date, first_detect))
  
  ## get distribution of travel times by reach and release group
  travel_times_out = travel_times_dist(data_agg, reaches_list, fish_df)
  travel_times_full = travel_times_out$travel_times_full
  travel_times_summary = travel_times_out$travel_times_summary
  
  
  ## note if specific receiver groups weren't operating in specific years
  
  
  
  
}



## for each reach and each release group, finds all the fish with known travel times and returns the distribution
travel_times_dist = function(data_agg, reaches_list, fish_df){
  
  reach_year_df = expand.grid(release_group = unique(fish_df$release_group), reach = 1:(length(reaches_list)-1), n_complete_fish = NA,
                              median_travel_time = NA, mean_travel_time = NA) 
  reach_year_df$row_id = 1:nrow(reach_year_df)
  travel_times_full = data.frame()
  
  for(i in unique(reach_year_df$row_id)){
    
    release_group_i = reach_year_df$release_group[i]
    reach_i = reach_year_df$reach[i]
    fish_df_i = subset(fish_df, release_group == release_group_i)
    
    data_agg_i = data_agg %>% filter(fish_id %in% fish_df_i$fish_id & 
                                     (rec_group_id %in% c(reach_i, reach_i + 1)) &
                                     first_reach <= reach_i) %>%
                              mutate(type = if_else(rec_group_id == reach_i, 'enter', 'leave'))
    
    if(nrow(data_agg_i) == 0){
      reach_year_df$n_complete_fish[i] = 0
    }else{
      data_agg_i_wide = data_agg_i %>% select(fish_id, type, first_detect) %>%
        pivot_wider(names_from = type, values_from = first_detect) %>%
        filter(!is.na(enter) & !is.na(leave))
      
      reach_year_df$n_complete_fish[i] = nrow(data_agg_i_wide)
      
      if(nrow(data_agg_i_wide) > 0){
        
        #calculate travel time (in hours)
        data_agg_i_wide$travel_time = as.numeric(difftime(data_agg_i_wide$leave, data_agg_i_wide$enter, units = 'hours'))
        
        #save summary stats of distribution
        reach_year_df$median_travel_time[i] = median(data_agg_i_wide$travel_time)
        reach_year_df$mean_travel_time[i] = mean(data_agg_i_wide$travel_time)
        
        #save all travel time data
        data_agg_i_wide$reach = reach_i
        data_agg_i_wide$release_group = release_group_i
        travel_times_full = rbind(travel_times_full, data_agg_i_wide)
      }
    }
    
    
  }
  
  
  return(list(travel_times_full = travel_times_full, travel_times_summary = reach_year_df))
}


