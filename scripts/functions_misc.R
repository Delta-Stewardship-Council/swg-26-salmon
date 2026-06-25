## This function plots the river network based on subsequent fish observations
## For each fish leaving a general receiver location (i), this function finds the 
## next location (j) a fish is observed. The river network is then plotted, where
## the width of each line is based on the percent of fish leaving i that are next
## observed at receiver location j 
## Inputs:
##    data: full set of telemetry observations, including receiver location info
##    study_i: optional, specify a study id to filter by
##    min_count: if fewer than X fish move between i and j, disregard this connection
##    filter_perc: if fewer than X proportion of fish leaving i go to j (after filtering
##                 for min_count), the line between i and j is not plotted.
## Outputs: ggplot of river network

library(reshape2); library(igraph); library(ggplot2)

plot_river_network = function(data, study_i = NULL, filter_perc = 5, min_count = 2, min_rec_dist){
  
  if(!is.null(study_i)){
    data = subset(data, study_id == study_i)
  }
  
  ## double check data is sorted by date and fish
  data = data[order(data$fish_id, data$first_time), ]
  
  
  ## pull out release locations
  release_info = data[!duplicated(data$release_location), c("release_location",
                                                             "release_latitude",
                                                             "release_longitude",
                                                             "release_river_km") ]
  
  ## pull out receiver info
  recv_info = data[!duplicated(data$receiver_general_location), c("receiver_general_location",
                                                                  "receiver_general_latitude",
                                                                  "receiver_general_longitude",
                                                                  "receiver_general_river_km")]
  ## group receivers that are sufficiently close
  #recv_info$new_rec_name = case_match(recv_info$receiver_general_location,
  #                                    c("GoldenGateE", "GoldenGateW") ~ "GoldenGate",
  #                                    c("BeniciaE", "BeniciaW") ~ "Benicia",
  #                                    c("Chipps_1", "Chipps_2", "ChippsW", "ChippsE") ~ "Chipps",
  #                                    c("Three_Mile_North", "Three_Mile_South", "ThreeMile") ~ "ThreeMile",
  #                                    c("JerseyPoint_1", "JerseyPoint_2") ~ "JerseyPoint",
  #                                    c("Sac_Rio_Vista_1", "Sac_Rio_Vista_2") ~ "Sac_Rio_Vista",
  #                                    c("FalseRiver_1", "FalseRiver_2") ~ "FalseRiver",
  #                                    .default = recv_info$receiver_general_location)
  
  ## group receivers together that are really close
  recv.dist =dist(recv_info[, "receiver_general_river_km"])
  recv.dist=as.matrix(recv.dist, labels=TRUE)
  colnames(recv.dist) <- rownames(recv.dist) <- recv_info$receiver_general_location 
  recv.dist[recv.dist < 1] = 0
  recv.dist[recv.dist > 1] = 1
  recv.pairs = reshape2::melt(recv.dist)
  recv.pairs = subset(recv.pairs, value == 0)

  idf <- graph_from_data_frame(recv.pairs)
  recv_info = merge(recv_info, stack(components(idf)$membership), by.x = "receiver_general_location", by.y = "ind", all.x = TRUE)
  colnames(recv_info)[colnames(recv_info)=="values"] = "rec_group"
  
  ## merge back in with full dataset
  data = merge(data, recv_info[,c("receiver_general_location", "rec_group")], by = "receiver_general_location", sort = F)

  ## again double check data is sorted by date and fish
  data = data[order(data$fish_id, data$first_time), ]
  
  ## find general location of next receiver in time series
  data = data %>% group_by(fish_id) %>%
         mutate(next_rec_group = lead(rec_group, n = 1)) %>%
         ungroup()
  
  ## filter observations where the fish doesn't move between general locations
  data_sub = subset(data, next_rec_group != rec_group)
  
  ## find % of fish that go to receiver j after leaving receiver i
  data_recv = data_sub %>% group_by(rec_group, next_rec_group) %>%
                           summarise(count_next = n(), .groups = "drop_last") %>% 
                           filter(count_next >= min_count) %>%
                           mutate(perc_next = (count_next / sum(count_next)) * 100)
  
  ## merge receiver info
  data_recv_sub = subset(data_recv, perc_next > filter_perc)
  data_recv_sub = merge(data_recv_sub, recv_info[!duplicated(recv_info$rec_group),], by = 'rec_group')
  data_recv_sub = merge(data_recv_sub, recv_info[!duplicated(recv_info$rec_group),], by.x = 'next_rec_group',
                        by.y = 'rec_group', suffixes = c(".i", '.j'))

  
  
  ### Plotting
  ggplot(data_recv_sub) + 
    geom_segment(aes(x = receiver_general_longitude.i, y = receiver_general_latitude.i,
                     xend = receiver_general_longitude.j, yend = receiver_general_latitude.j,
                     linewidth = perc_next),
                 arrow = arrow(length = unit(0.2,"cm")))+
    #geom_point(aes(x = receiver_general_longitude.i, y = receiver_general_latitude.i)) +
    #geom_point(aes(x = receiver_general_longitude.j, y = receiver_general_latitude.j))+
    scale_linewidth_continuous(range = c(0.01, 1), breaks = c(0, 10, 50, 100))+
    geom_point(data = release_info, aes(x = release_longitude, y = release_latitude, color = release_location))+
    coord_fixed() + theme_classic()              

} 
