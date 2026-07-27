### Map the important env and release locations

## Set up
library("ggplot2")
library("dplyr")
library("lubridate")
library("ggspatial") 
library("prettymapr")
library("sf")
library("ggrepel")
library("viridis")
library("tidyverse")
library("here")


######## Just trib stations zoom

coords_file = here("data_raw", "env_map_coords.csv")
landmarks <- read.csv(coords_file, stringsAsFactors = FALSE)

# Make datapoints into an sf
sta_to_sf <- st_as_sf(landmarks, coords= c("lon", "lat"), crs = 4326)
sta_3857_sf <- st_transform(sta_to_sf, 
  crs = 3857)
#sta_sf <- st_transform(sta_sf, st_crs(Deltaz))
#sta_wq_sf <- st_as_sf(sta[sta$WQ.type != "" ,], coords= c("Lon", "Lat"), crs = 4326)


env_map <- ggplot() + #works
  geom_sf(data= sta_to_sf, aes(color = type))+
  ggspatial::annotation_map_tile(type = "cartolight", zoom = 10)+
  geom_sf(data= sta_to_sf, aes(color = type), size =2)+
  coord_sf(xlim = c(-123.0, -120.7), ylim= c(37.7, 39.6))+
  scale_color_manual(name= "Data Type", values = c("#31688EFF", "#35B779FF"))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()

ggsave(
  here("figures/exploratory", "env_map.png"),
  plot = last_plot(),
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

################# All things together, full map 
# using env_map_coords_2.csv

coords_file = here("data_raw", "env_map_coords_2.csv")
landmarks <- read.csv(coords_file, stringsAsFactors = FALSE)

# Make datapoints into an sf
sta_to_sf <- st_as_sf(landmarks, coords= c("lon", "lat"), crs = 4326)
sta_3857_sf <- st_transform(sta_to_sf, 
  crs = 3857)


env_map_2 <- ggplot() + #works
  geom_sf(data= sta_to_sf, aes(color = type))+
  ggspatial::annotation_map_tile(type = "cartolight", zoom = 10)+
  geom_sf(data= sta_to_sf, aes(color = type), size =2)+
  coord_sf(xlim = c(-123.0, -120.5), ylim= c(37.7, 39.6))+
  scale_color_manual(name= "Data Type", values = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF"))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()

ggsave(
  here("figures/exploratory", "env_map_2.png"),
  plot = last_plot(),
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
