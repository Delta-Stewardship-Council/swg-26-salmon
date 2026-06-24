### Map the important env and release locations

## Set up
library("ggplot2")
library("dplyr")
library("lubridate")
library("leaflet")
library("ggspatial") 
library("prettymapr")
#library("ggmap")
#library("maps")
#library("mapdata")
library("deltamapr")
library("sf")
library("ggrepel")
#library("scatterpie")
library("viridis")
library("tidyverse")

coords_file = here("data_raw", "env_map_coords.csv")
landmarks <- read.csv(coords_file, stringsAsFactors = FALSE)


######### Playing around with Maps
WW_Delta <- deltamapr::WW_Delta
WW_Watershed <- deltamapr::WW_Watershed
SFB <- deltamapr::WW_Watershed[c(418, 432:433),]
#SFB$HNAME <- "SFB"
SuisunB <- deltamapr::WW_Watershed[c(308:309, 320:325, 332:334, 411:412, 414:415),]
check <- deltamapr::WW_Watershed[193,]
River <- deltamapr::WW_Watershed[c(337:338, 343,353, 367:384, 402:403, 409, 421, 425:427, 430:431, 436:437),]
Deltaw <- deltamapr::WW_Watershed[c(161:331, 438),] #c(161:163, 202, 213, 215, 237, 267, 288:289, 294:297, 311:315, 352, 438),]
#Folsom lake 434 not included right now
ggplot(WW_Watershed)+
  geom_sf(aes())+
  theme_bw()

ggplot(SFB)+
  geom_sf(fill="red")+
  theme_bw()

Deltaz <- st_crop(WW_Delta, xmin = -122.15, xmax= -121.3, 
                  ymin= 37.7, ymax= 38.7)

# Make datapoints into an sf
sta_to_sf <- st_as_sf(landmarks, coords= c("lon", "lat"), crs = 4326)
sta_3857_sf <- st_transform(sta_to_sf, 
  crs = 3857)
#sta_sf <- st_transform(sta_sf, st_crs(Deltaz))
#sta_wq_sf <- st_as_sf(sta[sta$WQ.type != "" ,], coords= c("Lon", "Lat"), crs = 4326)


sta_flow_xy <- cbind(
  sta_flow_sf,
  st_coordinates(sta_flow_sf)
)

sta_wq_xy <- cbind(
  sta_wq_sf,
  st_coordinates(sta_wq_sf)
)

sta_sed_xy <- cbind(
  sta_sed_sf,
  st_coordinates(sta_sed_sf)
)


#This gets you zoomed in delta - FLOW sta
ggplot(Deltaz)+
  geom_sf()+
  geom_sf(data= sta_flow_sf, aes(color= Flow.type), size= 2)+
  geom_label_repel(data= sta_flow_xy, aes(x= X, y= Y, label=Sta.id), box.padding = 0.3, point.padding = 0.4,
                   min.segment.length =0, size= 2, fontface= "bold")+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_color_discrete(name ="Flow Data Type")+
  theme_bw()

ggplot(WW_Watershed) +
  geom_sf(aes())+
  geom_sf(data = sta_to_sf, 
    mapping = aes(fill = type)) +
    geom_sf(data = rivers_3338_sf,
            aes(linewidth = StrOrder),
            color = 'blue', alpha = .5) +
    scale_linewidth(range = c(0.05, 0.5),
                    guide = "none") +
    geom_sf(data = pop_3338_sf, size = 0.5, color = 'grey20') +
    labs(title = "Total Population by Alaska Region",
         fill = "Total Population") +
    scale_fill_continuous(low = "khaki",
                          high =  "firebrick") +
    theme_bw()


  ggplot(data = sta_to_sf) +
    ggspatial::annotation_map_tile(type = "osm", zoom = 6, progress = 'none') + 
    # higher zoom values are more detailed 
    #geom_sf(aes(color = population),
    #        fill = NA) +
    coord_sf(xlim = c(-16900000, -16300000), ylim = c(8400000, 9000000))

ggplot() +
  ggspatial::annotation_map_tile(type = "cartolight", zoom = 16, progress = 'none')+
  geom_sf(data= sta_to_3857, aes(color = type)) +
  coord_sf(xlim = c(-13500000, -13600000), ylim = c(4500000, 4800000))

ggplot(data = sta_to_sf) +
  geom_sf()

ak_bbox <- st_bbox(pop_4326_sf) %>%
    as.numeric()

ggplot() + #works
  ggspatial::annotation_map_tile(type = "cartolight", zoom = 10)+
  geom_sf(data= sta_to_sf, aes(color = type))+
  coord_sf(xlim = c(-13500000, -13600000), ylim = c(4500000, 4800000))
  #coord_sf(xlim = c(-122.15, -121.3), ylim= c(37.7, 39.6))

  ggplot() + #works
    geom_sf(data= sta_to_sf, aes(color = type))+
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 10)+
    geom_sf(data= sta_to_sf, aes(color = type))#+
    #coord_sf(xlim = c(-13500000, -13600000), ylim = c(4500000, 4800000))

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

### All things together, full map 
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
