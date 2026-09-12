## This function takes the dat_fish_recv dataframe (downloaded and processed with the
## download_process_telemetry function) and processes it further. Detections for each fish
## are aggregated by reach, and unknown travel times are sampled.

library(dplyr)

aggregate_reaches = function(dat_fish_recv, reaches_list, excl_reaches_studies_list){
  
  ## remove unneeded detections
  data = dat_fish_recv %>% filter(receiver_general_location %in% unique(unlist(reaches_list))) %>%
                           arrange(fish_id, first_time)
  
  ## summary data frame of fish
  fish_df = data[!duplicated(data$fish_id), c("fish_id", 'study_id', 'fish_type', 'fish_origin', 'fish_date_tagged',
                                              'fish_release_date','release_location') ]
  # release group
  fish_df$release_group = paste0(fish_df$study_id, "_", fish_df$release_location, "_", format(fish_df$fish_release_date, "%j"))
  
  # find which reach each fish starts in
  fish_df$first_reach = sapply(fish_df$release_location, function(val) {
    which(sapply(reaches_list, function(x) val %in% x))
  })
  
  ## data frame of receiver location names defining each reach
  reaches_df = data.frame(receiver_general_location = unlist(reaches_list),
                          rec_group_id = rep(seq_along(reaches_list), lengths(reaches_list)))
  
  data = left_join(data, reaches_df, by = 'receiver_general_location')
  
  ## For each fish, find first detection in each reach
  data_agg = data %>% group_by(fish_id, rec_group_id) %>%
                             summarize(first_detect = min(first_time)) %>%
                             ungroup() %>% complete(fish_id, rec_group_id)
  data_agg = left_join(data_agg, fish_df[, c("fish_id", "first_reach", 'fish_release_date'), by = 'fish_id'])
  data_agg = subset(data_agg, rec_group_id >= first_reach) #ignore earlier reaches for fish released in later reaches
  
  # fill in first detection time for fish not detected at release site
  data_agg = data_agg %>% mutate(first_detect_processed = if_else(is.na(first_detect) & rec_group_id == first_reach,
                                                          fish_release_date, first_detect))
  
  ## get distribution of known travel times by reach and release group
  travel_times_out = travel_times_dist(data_agg, reaches_list, fish_df)
  travel_times_full = travel_times_out$travel_times_full
  travel_times_summary = travel_times_out$travel_times_summary
  
  
  ## note if specific receiver groups weren't operating in specific years
  
  
  ## associate flow data to fish with known travel times (travel_times_full)
  
  
  ## for each reach, fit lognormal distribution to travel time data, dependent on flow
  
  

  
}



## for each reach and each release group, finds all the fish with known travel times and returns the distribution
travel_times_dist = function(data_agg, reaches_list, fish_df){
  
  release_summary = fish_df[!duplicated(fish_df$release_group), c("release_group", "study_id", 'release_location', 'fish_release_date')]
  
  reach_year_df = expand.grid(release_group = unique(fish_df$release_group), reach = 1:(length(reaches_list)-1), n_complete_fish = NA,
                              median_travel_time = NA, mean_travel_time = NA) 
  reach_year_df = left_join(reach_year_df, release_summary, by = 'release_group')
  
  
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



## for each fish, fill in missing travel times
fill_travel_times = function(data_agg, fish_df, travel_times_summary, travel_times_full){

  for(fish_i in unique(fish_df$fish_id)){
    
    data_agg_i = data_agg %>% filter(fish_id == fish_i) %>%
      mutate(enter = first_detect, leave = lead(first_detect),
             travel_time = as.numeric(difftime(leave, enter, units = 'hours')),
             travel_time_cumul = as.numeric(difftime(first_detect_processed,
                                                     first(first_detect_processed), units = 'hours')),
             after_last_detect = ifelse(rec_group_id >= as.numeric(data_agg_i[which.max(data_agg_i$travel_time_cumul), "rec_group_id"]), 1, 0),
             fill_group_id = consecutive_id(travel_time, after_last_detect))
    
    for(row_i in 1:nrow(data_agg_i)){
      
      if(is.na(data_agg_i$travel_time[row_i]) & data_agg_i$after_last_detect[row_i] == F){
        
        ## identify travel times to be drawn jointly
        data_to_fill = subset(data_agg_i, fill_group_id == data_agg_i$fill_group_id[row_i])
        
        ## identify total travel time they must sum to
        total_travel_i = as.numeric(difftime(data_to_fill$leave[!is.na(data_to_fill$leave)],
                                             data_to_fill$enter[!is.na(data_to_fill$enter)], units = 'hours'))
        
        ## draw travel times with rejection sampling 
        
        
      }else if(data_agg_i$after_last_detect[row_i] == T){
        ## draw missing travel time after final detect (no need for rejection sampling)
        
        
        
      }
      
      
      
      
    }
    
    
  }
  
  
}





# for each release group, find what reaches had at least one receiver present for the period of at least X days following release
get_reach_coverage = function(recvdat, fishdat, reaches_list, minDays = 45){
  
  ## get basic fish info
  fish_df = fishdat[, c("fish_id", 'study_id', 'fish_type', 'fish_origin', 'fish_date_tagged',
                        'fish_release_date','release_location')]
  # release group
  fish_df$release_group = paste0(fish_df$study_id, "_", fish_df$release_location, "_", format(fish_df$fish_release_date, "%j"))
  
  # find which reach each fish starts in
  fish_df$first_reach = sapply(fish_df$release_location, function(val) {
    which(sapply(reaches_list, function(x) val %in% x))
  })
  
  # summary of fish releases
  release_summary = fish_df[!duplicated(fish_df$release_group), c("release_group", "study_id", 'release_location', 'fish_release_date', 'first_reach')]
  
  ## data frame of receiver location names defining each reach
  reaches_df = data.frame(receiver_general_location = unlist(reaches_list),
                          rec_group_id = rep(seq_along(reaches_list), lengths(reaches_list)))
  
  ## combinations of receiver groups and release groups
  recvs_coverage = expand.grid(rec_group_id = unique(reaches_df$rec_group_id), release_group = unique(fish_df$release_group),
                               complete_coverage = NA, partial_coverage = NA, coverage_start = NA, coverage_end = NA)

  recvs_coverage = left_join(recvs_coverage, release_summary, by = 'release_group')
  recvs_coverage = subset(recvs_coverage, first_reach <= rec_group_id)
  recvs_coverage$final_day = recvs_coverage$fish_release_date + minDays*24*60*60 #how far from release we should ideally have receiver coverage
  recvs_coverage$row_id = 1:nrow(recvs_coverage)
  
  ## for each receiver group/release group pair, see if the date range (release date -> release date + X days) is fully, 
  ##      partially, or not at all covered by at least one receiver
  for(j in 1:nrow(recvs_coverage)){
    reaches_j = subset(reaches_df, rec_group_id == recvs_coverage$rec_group_id[j])$receiver_general_location
    
    #find receiver deployments with at least some overlap with fish date range
    recvdat_j = recvdat %>% filter(receiver_general_location %in% reaches_j) %>%
      mutate(recv_start= as.POSIXct(receiver_start, format = "%m/%d/%Y %H:%M:%S"), 
             recv_end = as.POSIXct(receiver_end, format = "%m/%d/%Y %H:%M:%S"),
             before_release = ifelse(recv_start < recvs_coverage$fish_release_date[j], 1, 0), 
             after_end = ifelse(recv_end > recvs_coverage$final_day[j] , 1, 0), 
             complete_coverage = before_release & after_end) %>%
      filter(recv_start < recvs_coverage$final_day[j] & recv_end > recvs_coverage$fish_release_date[j])
    
    if(nrow(recvdat_j) == 0){
      recvs_coverage$complete_coverage[j] = F
      recvs_coverage$partial_coverage[j] = F
    }else if(sum(recvdat_j$complete_coverage) >= 1){
      recvs_coverage$complete_coverage[j] = T
      recvs_coverage$partial_coverage[j] = T
    }else{
      recvs_coverage$partial_coverage[j] = T
      
      ## see if receivers piece together to form total coverage
      recvdat_j_comb = recvdat_j %>%
                       select(dep_id, recv_start, recv_end) %>%
                       mutate(recv_start_date = as.Date(recv_start),
                              recv_end_date = as.Date(recv_end),
                              covg_range = iv(recv_start_date, recv_end_date), .keep = "unused") %>%
                       summarise(comb_range = iv_groups(covg_range)) %>%
                       mutate(range_start = iv_start(comb_range), 
                              range_end = iv_end(comb_range),
                              before_release = ifelse(range_start <= as.Date(recvs_coverage$fish_release_date[j]), 1, 0), 
                              after_end = ifelse(range_end >= as.Date(recvs_coverage$final_day[j]) , 1, 0),
                              complete_coverage = before_release & after_end)
      
      if(sum(recvdat_j_comb$complete_coverage) >= 1){
        recvs_coverage$complete_coverage[j] = T
      }else{
        recvs_coverage$complete_coverage[j] = F
        recvs_coverage$coverage_start[j] = min(recvdat_j_comb$range_start)
        recvs_coverage$coverage_end[j] = max(recvdat_j_comb$range_end)
      }
      
      
    }
    
  }
  
  recvs_coverage$coverage_start = as.Date(recvs_coverage$coverage_start)
  recvs_coverage$coverage_end = as.Date(recvs_coverage$coverage_end)
  
  return(recvs_coverage)
  
}

