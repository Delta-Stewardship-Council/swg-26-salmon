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

library(reshape2); library(igraph); library(ggplot2); library(ggrepel)
#study_i = "FR_Spring_2013"; filter_perc = 5; min_count = 2; min_rec_dist = 1

plot_river_network = function(data, study_i = NULL, filter_perc = 5, min_count = 2, min_rec_dist = 1){
  
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
  recv.dist[recv.dist < min_rec_dist] = 0
  recv.dist[recv.dist >= min_rec_dist] = 1
  recv.pairs = reshape2::melt(recv.dist)
  recv.pairs = subset(recv.pairs, value == 0)

  idf <- graph_from_data_frame(recv.pairs)
  recv_info = merge(recv_info, stack(components(idf)$membership), by.x = "receiver_general_location", by.y = "ind", all.x = TRUE)
  recv_info = recv_info %>% arrange(receiver_general_river_km) %>% 
              mutate(rec_group = consecutive_id(values)) %>% ungroup()
  recv_info$values = NULL
  recv_info$rec_bin = cut(recv_info$receiver_general_latitude, breaks = seq(min(recv_info$receiver_general_latitude)-.1, max(recv_info$receiver_general_latitude)+.1, by = .1))

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
  
  ## filter by % 
  data_recv_sub = subset(data_recv, perc_next > filter_perc)
  
  ## assign receiver labels for plotting
  recv_info_sub = subset(recv_info, rec_group %in% unique(data_recv_sub$rec_group) | rec_group %in% unique(data_recv_sub$next_rec_group))
  recv_bins = recv_info_sub[!duplicated(recv_info_sub$rec_group), ]
  recv_bins = recv_bins %>% group_by(rec_bin) %>%
    mutate(plot_x = consecutive_id(receiver_general_longitude),
           plot_x_sc = scale(plot_x, scale = F)) %>%
    ungroup()
  recv_info_sub = merge(recv_info_sub, recv_bins[, c("rec_group", "plot_x", "plot_x_sc")], by = 'rec_group')
  
  
  # merge with receiver info
  data_recv_sub = merge(data_recv_sub, recv_info_sub[!duplicated(recv_info_sub$rec_group),], by = 'rec_group')
  data_recv_sub = merge(data_recv_sub, recv_info_sub[!duplicated(recv_info_sub$rec_group),], by.x = 'next_rec_group',
                        by.y = 'rec_group', suffixes = c(".i", '.j'))

  ## most important receivers
  recv_important = subset(data_recv_sub, count_next > 10 | rec_group == 2) #subset(data_recv_sub, rec_group %in% c(1, 2, max(data_recv_sub$rec_group))| count_next >= sort(data_recv_sub$count_next, decreasing = T)[10])
  recv_important = subset(recv_info_sub, rec_group %in% unique(recv_important$rec_group) | rec_group %in% unique(recv_important$next_rec_group))
  recv_important = subset(recv_important, !duplicated(rec_group))
  
  ### Plotting
  
  ## by lat/long
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
  
  ## simplified network map
  ggplot(data_recv_sub) + 
    geom_segment(aes(y = rec_group, x = plot_x_sc.i,
                     yend = next_rec_group, xend = plot_x_sc.j,
                     linewidth = perc_next),
                 arrow = arrow(length = unit(0.2,"cm")))+
    scale_linewidth_continuous(range = c(0.01, 1), breaks = c(0, 10, 50, 100))+
    #geom_point(data = release_info, aes(x = release_longitude, y = release_river_km, color = release_location))+
    geom_text_repel(data = recv_important, aes(y = rec_group, x = plot_x_sc, label = receiver_general_location),
                    nudge_x = ifelse(recv_important$plot_x_sc>0, 2, -2), color = 'red')+
    theme_classic() + labs(x = NULL, y = NULL, linewidth = "% of fish")

} 
