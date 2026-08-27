## This function plots the river network based on subsequent fish observations
## For each fish leaving a general receiver location (i), this function finds the 
## next location (j) a fish is observed. The river network is then plotted, where
## the width of each line is based on the percent of fish leaving i that are next
## observed at receiver location j 
## Inputs:
##    data: full set of telemetry observations, including receiver location info ('dat_fish_recv' object)
##    study_i: optional, specify a study id to filter by
##    min_count: if fewer than X fish move between i and j, disregard this connection
##    filter_perc: if fewer than X proportion of fish leaving i go to j (after filtering
##                 for min_count), the line between i and j is not plotted.
##    min_rec_dist: general receiver locations within this distance (km) will be grouped together
##    plot_only_consecutive: if T, only shows paths between consecutive receivers; if F, plots all fish paths 
##    plot_all_names: set to TRUE to label all receiver general locations on the network map; FALSE only 
##                    labels one for each group
## Outputs: ggplot of river network

library(reshape2); library(igraph); library(ggplot2); library(ggrepel); library(this.path); library(geosphere)
setwd(this.path::here()); setwd('..')
#study_i = NULL; filter_perc = 10; min_count = 2; min_rec_dist = 1; save_dir = 'figures'

plot_river_network = function(data, study_i = NULL, filter_perc = 10, min_count = 2, min_rec_dist = 1,
                              plot_only_consecutive = F, plot_all_names = T, save_dir = 'figures'){
  
  ## Set up folder to save (if it doesn't exist)
  if(!dir.exists(save_dir)){ dir.create(save_dir) }
  
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
                                                                  "receiver_general_river_km",
                                                                  'receiver_region')]
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
  #recv.dist =dist(recv_info[, "receiver_general_river_km"])
  recv.dist =distm(recv_info[, c("receiver_general_longitude", "receiver_general_latitude")], fun = distHaversine)
  recv.dist=as.matrix(recv.dist, labels=TRUE)
  recv.dist = recv.dist/1000 # m to km
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
  
  ## identify receiver groups present in most studies
  recv_studies = data %>% group_by(rec_group)  %>% 
                 summarise(unique_studies = n_distinct(study_id),
                           n_fish = n_distinct(fish_id),
                           study_ids = paste(unique(study_id), collapse = ', '))
  recv_info = left_join(recv_info, recv_studies, by = 'rec_group')
  
  ## full receiver map
  full_recv_map = ggplot(recv_info) + 
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 10)+
    coord_sf(xlim = c(-123.0, -120.7), ylim= c(37.7, 39.6), crs = 4326)+
    geom_spatial_point(aes(x = receiver_general_longitude, y = receiver_general_latitude))+
    theme_classic() + labs(x = NULL, y = NULL)
  
  full_recv_map2 = ggplot(subset(recv_info, unique_studies >= 7)) + 
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 10)+
    coord_sf(xlim = c(-123.0, -120.7), ylim= c(37.7, 39.6), crs = 4326)+
    geom_spatial_point(aes(x = receiver_general_longitude, y = receiver_general_latitude))+
    theme_classic() + labs(x = NULL, y = NULL)
  
  
  ## focus on only important receivers
  recv_studies = subset(recv_studies, unique_studies >= length(unique(data$study_id)) - 2)
  data = subset(data, rec_group %in% unique(recv_studies$rec_group))
  recv_info = subset(recv_info, rec_group %in% unique(recv_studies$rec_group))
  
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
                           mutate(perc_next = (count_next / sum(count_next)) * 100,
                                  rec_pair = paste(rec_group, next_rec_group, sep = "_"))
  
  ## find unique fish paths of *consecutive* receivers
  conscpath = data.frame()
  for(rec_i in unique(data_recv$rec_group)){
    nextgroups = subset(data_recv, rec_group == rec_i & next_rec_group < rec_group)
    next2groups = subset(data_recv, rec_group %in% unique(nextgroups$next_rec_group) & next_rec_group < rec_group)
    conscpath_i = subset(nextgroups, !(next_rec_group %in% unique(next2groups$next_rec_group)))
    conscpath = rbind(conscpath, conscpath_i)
  }
  
  ## filter by % 
  data_recv_sub = subset(data_recv, perc_next > filter_perc)
  
  ## merge consecutive path info
  data_recv_sub$consecutive = ifelse(data_recv_sub$rec_pair %in% conscpath$rec_pair, 1, 0)
  if(plot_only_consecutive == T){
    data_recv_sub = subset(data_recv_sub, consecutive == 1)
  }
  
  ## assign receiver labels for plotting
  recv_info_sub = subset(recv_info, rec_group %in% unique(data_recv_sub$rec_group) | rec_group %in% unique(data_recv_sub$next_rec_group))
  recv_bins = recv_info_sub[!duplicated(recv_info_sub$rec_group), ]
  recv_bins = recv_bins %>% group_by(rec_bin) %>%
    mutate(plot_x = consecutive_id(receiver_general_longitude),
           plot_x_sc = scale(plot_x, scale = F),
           plot_y = consecutive_id(rec_group)) %>%
    ungroup()
  recv_info_sub = merge(recv_info_sub, recv_bins[, c("rec_group", "plot_x", "plot_x_sc", "plot_y")], by = 'rec_group')
  
  
  # merge with receiver info
  data_recv_sub = merge(data_recv_sub, recv_info_sub[!duplicated(recv_info_sub$rec_group),], by = 'rec_group')
  data_recv_sub = merge(data_recv_sub, recv_info_sub[!duplicated(recv_info_sub$rec_group),], by.x = 'next_rec_group',
                        by.y = 'rec_group', suffixes = c(".i", '.j'))

  ## most important receivers
  recv_important = subset(data_recv_sub, count_next > 10 | rec_group == 2) #subset(data_recv_sub, rec_group %in% c(1, 2, max(data_recv_sub$rec_group))| count_next >= sort(data_recv_sub$count_next, decreasing = T)[10])
  recv_important = subset(recv_info_sub, rec_group %in% unique(recv_important$rec_group) | rec_group %in% unique(recv_important$next_rec_group))
  recv_important_all = recv_important %>% group_by(rec_group) %>% summarize(all_names = paste0(receiver_general_location, collapse = ', '))
  recv_important = subset(recv_important, !duplicated(rec_group))
  recv_important = left_join(recv_important, recv_important_all, by = 'rec_group')
  
  ## regions
  recv_important$receiver_general_region = ifelse(recv_important$receiver_region %in% c("Feather_R", "Yolo Bypass"), "Feather River", NA)
  recv_important$receiver_general_region = ifelse(recv_important$receiver_region %in% c("Lower Sac R"), "Lower Sac", recv_important$receiver_general_region)
  recv_important$receiver_general_region = ifelse(recv_important$receiver_region %in% c("North Delta", "West Delta"), "Delta", recv_important$receiver_general_region)
  recv_important$receiver_general_region = ifelse(recv_important$receiver_region %in% c("Carquinez Strait", "SF Bay"), "SF Bay", recv_important$receiver_general_region)
  
  ### Plotting
  
  ## by lat/long
  full_network_map = ggplot(data_recv_sub) + 
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 10)+
    coord_sf(xlim = c(-123.0, -120.7), ylim= c(37.7, 39.6), crs = 4326)+
    geom_spatial_segment(aes(x = receiver_general_longitude.i, y = receiver_general_latitude.i,
                     xend = receiver_general_longitude.j, yend = receiver_general_latitude.j),
                     linewidth = 1,
                 arrow = arrow(length = unit(0.2,"cm")))+
    geom_spatial_point(data = release_info, aes(x = release_longitude, y = release_latitude, color = release_location))+
    geom_spatial_label_repel(data = recv_important, aes(x = receiver_general_longitude, y = receiver_general_latitude,
                                                        label = receiver_general_location),
                     nudge_x = ifelse(recv_important$plot_x_sc>0 | (recv_important$plot_x_sc==0 & recv_important$rec_group%%2!=0), .5, -.5), size = 3) +
    theme_classic() + labs(x = NULL, y = NULL, linewidth = NULL, color = 'Release Location')
  
  delta_zoom = ggplot(subset(data_recv_sub, !(receiver_region.j %in% c("SF Bay", "Feather_R") ))) + 
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 10)+
    coord_sf(xlim = c(-122.3, -121.2), ylim= c(38, 38.9), crs = 4326)+
    geom_spatial_segment(aes(x = receiver_general_longitude.i, y = receiver_general_latitude.i,
                             xend = receiver_general_longitude.j, yend = receiver_general_latitude.j),
                         linewidth = 1,
                         arrow = arrow(length = unit(0.2,"cm")))+
    geom_spatial_label_repel(data = subset(recv_important, !(receiver_region %in% c("SF Bay", "Feather_R") ) ), 
                             aes(x = receiver_general_longitude, y = receiver_general_latitude,
                                 label = receiver_general_location),
                             nudge_x = ifelse(subset(recv_important, !(receiver_region %in% c("SF Bay", "Feather_R") ) )$rec_group%%2!=0, .5, -.5), size = 3) +
    theme_classic() + labs(x = NULL, y = NULL, linewidth = NULL) + theme(legend.position = 'none')

  
  ## simplified network map
  if(plot_all_names == T){
    network_plot = ggplot(data_recv_sub) + 
      geom_label_repel(data = recv_important, aes(y = rec_group, x = plot_x_sc, label = stringr::str_wrap(all_names, 20), color = receiver_general_region),
                       nudge_x = ifelse(recv_important$plot_x_sc>0 | (recv_important$plot_x_sc==0 & recv_important$rec_group%%2!=0), 2, -2), size = 3)
  }else{
    network_plot = ggplot(data_recv_sub) + 
      geom_text_repel(data = recv_important, aes(y = rec_group, x = plot_x_sc, label = receiver_general_location, color = receiver_general_region),
                       nudge_x = ifelse(recv_important$plot_x_sc>0 | (recv_important$plot_x_sc==0 & recv_important$rec_group%%2!=0), 2, -2))
  }
  
network_plot = network_plot +
    geom_segment(aes(y = rec_group, x = plot_x_sc.i,
                     yend = next_rec_group, xend = plot_x_sc.j,
                     linewidth = perc_next),
                 arrow = arrow(length = unit(0.2,"cm")))+
    scale_linewidth_binned(range = c(0.01, 1), n.breaks = 8)+
    #geom_point(data = release_info, aes(x = release_longitude, y = release_river_km, color = release_location))+
    theme_classic() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line  = element_blank(), legend.position = 'none')+
    labs(x = NULL, y = NULL, linewidth = "% of fish")

  save_plot_name1 = paste0(save_dir, "/network_maps_", gsub(" ", "", unique(data$fish_type)), "_",
                        "filter", filter_perc, "_", format(Sys.time(), "%h%d_%H%M%S"), ".png")
  save_plot_name2 = paste0(save_dir, "/network_maps_zoom_", gsub(" ", "", unique(data$fish_type)), "_",
                           "filter", filter_perc, "_", format(Sys.time(), "%h%d_%H%M%S"), ".png")
  save_plot_name3 = paste0(save_dir, "/network_plots_", gsub(" ", "", unique(data$fish_type)), "_",
                          "filter", filter_perc, "_", format(Sys.time(), "%h%d_%H%M%S"), ".png")
  ggsave(plot = full_network_map, filename = save_plot_name1, width = 8, height = 6) 
  ggsave(plot = delta_zoom, filename = save_plot_name2, width = 8, height = 6) 
  ggsave(plot = network_plot, filename = save_plot_name3, width = 8, height = 6) 
  
  return(network_plot)
} 
