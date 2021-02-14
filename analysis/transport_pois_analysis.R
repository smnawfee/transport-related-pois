# Heading -----------------------------------------------------------------
# Transport-related POI analysis
#Author: Shahreen
#Date: 15/01/2021


# Clear environment-------------------------------------------------------
rm(list = ls())


#Folder setup---------------------------------------------------------
data_folder = "data/"


#call library-------------------------------------------------
library(rgdal)
library(dplyr)
library(raster)
library(sf)
library(tmap)
library(spdep)


# Load data---------------------------------------------------------------------------
#load pre-processed data for extracting transport data
london_os_refined_clusters <- readRDS(paste0(data_folder, "london_os_refined_clusters.rds"))
london_osm_refined_clusters <- readRDS(paste0(data_folder, "london_osm_refined_clusters.rds"))
#load data for sub_categories
london_os_special_pois <- readRDS(paste0(data_folder, "london_os_special_POIS.rds"))
london_osm_special_pois <- readRDS(paste0(data_folder, "london_osm_special_POIS.rds"))
#load transport-csv data
os_trans_poi <- read.csv(paste0(data_folder, "os_transport1.csv"))
osm_trans_poi <- read.csv(paste0(data_folder, "osm_transport1.csv"))

##load boundary files
london <- readOGR("data/boundary_london", "boundary_london")
london_borough <- readOGR("data/boundary_london","London_Borough_Excluding_MHW")


#Rasterize london border----------------------------------------------------------
london_border <- st_cast(st_as_sf(london), "MULTILINESTRING")
raster_template <- raster(extent(london), 
                          resolution = 1000, 
                          crs = st_crs(london)$proj4string)
london_raster <- rasterize(london, raster_template) 
london_border_raster <- rasterize(london_border, raster_template) 
london_raster_with_borders <- merge(london_raster, london_border_raster)


#Transport POIs analysis-----------------------------------------------------------------
#OSM Transport POIs
london_osm_Transport <- filter(london_osm_refined_clusters, 
                               (OSM_Cluster_No %in% c("16"))) %>% 
                        st_as_sf(coords = c("long","lat"))
osm_t <- st_transform(london_osm_Transport, proj4string(london))
osm_t <- as(osm_t, "Spatial")

#Grid-Density for OSM Transport POIs
gd_osm_t <- mask(rasterize(coordinates(osm_t), 
                           london_raster_with_borders, 
                           fun = 'count', 
                           background = 0),
                 london_raster_with_borders)
transport_osm <- gd_osm_t / cellStats(gd_osm_t, stat = "max")
transport_osm <- reclassify(transport_osm, cbind(0,NA))

#OS Transport POIs
london_os_transport <- filter(london_os_refined_clusters, 
                              (OS_Cluster_No %in% c("16"))) %>%
                       st_as_sf(coords = c("long","lat"))
os_t <- st_transform(london_os_transport, proj4string(london))
os_t <- as(os_t, "Spatial")

#GD for OS Transport POIs
gd_os_t <- mask(rasterize(coordinates(os_t), 
                          london_raster_with_borders,
                          fun = 'count', background = 0), london_raster_with_borders)
transport_os <- gd_os_t / cellStats(gd_os_t, stat = "max")
transport_os <- reclassify(transport_os, cbind(0,NA))

#find difference between os & osm
t_os_osm <- ((transport_osm) - (transport_os))/max(cellStats(transport_os, 
                                                             stat = "sum"),  
                                                   cellStats(transport_osm, 
                                                             stat="sum"))
#Transport_barplot
#convert to dataframe
t_os_osm_df <- as.data.frame(t_os_osm, 
                             xy=T, 
                             na.rm=T)
#create layer with customize bins
t_os_osm_df <- t_os_osm_df %>%
               mutate(density = cut(layer, 
               breaks = c(-0.1, -0.005, -0.001,-0.0001, 0.0001, 0.001,  0.005, 0.1)))

#Local g calculations
#convert to sf_polygons and then points
t <- rasterToPolygons(t_os_osm)
t <- st_as_sf(t)
t_po <- st_centroid(t)
#join neigbours within 0 and 3000 m 
nb_t <- dnearneigh(t_po, 0, 3000)
#create listw
nb_t <- nb2listw(nb_t)
#calculate local g
local_g_t <- localG(t$layer, nb_t)
local_g_t <- cbind(t, as.matrix(local_g_t))
#change field name
names(local_g_t)[2] <- "gstat"


#parking POI ---------------------------------
#OSM Parking POIs
london_osm_refined <- readRDS(file = "data/london_osm_refined.rds")
london_osm_parking <- filter(london_osm_refined, (Value %in% c("parking"))) %>%
                      st_as_sf(coords = c("long","lat"))
osm_p <- st_transform(london_osm_parking, proj4string(london))
osm_p <- as(osm_p, "Spatial")

#Grid-Density for OSM parking POI
gd_osm_p <- mask(rasterize(coordinates(osm_p), 
                           london_raster_with_borders, 
                           fun = 'count', 
                           background = 0),
                 london_raster_with_borders)
parking_osm <- gd_osm_p / cellStats(gd_osm_p, stat = "max")
parking_osm <- reclassify(parking_osm, cbind(0,NA))

#OS Parking POIs
london_os_refined <- readRDS(file = "data/london_os_refined.rds")
london_os_parking <- filter(london_os_refined, (classification_description %in% c("PARKING"))) %>%
                      st_as_sf(coords = c("long","lat"))
os_p <- st_transform(london_os_parking, proj4string(london))
os_p <- as(os_p, "Spatial")

#Grid-Density for OS parking POIs
gd_os_p <- mask(rasterize(coordinates(os_p), 
                          london_raster_with_borders, 
                          fun = 'count', 
                          background = 0),
                london_raster_with_borders)
parking_os <- gd_os_p / cellStats(gd_os_p, stat = "max")
parking_os <- reclassify(parking_os, cbind(0,NA))

#find difference between os & osm
p_os_osm <- ((parking_osm) - (parking_os))/max(cellStats(parking_os, 
                                                         stat = "sum"),  
                                               cellStats(parking_osm, 
                                                         stat="sum"))


#Plot the figures----------------------------------------------------------
#plot osm/os map function
osm_poi_map <- function (x) {
  tm_shape(x) +
  tm_raster(style = "fixed",
            n = 7,
            breaks = c(0,0.02, 0.05,0.1,0.2,0.4,0.7,1),
            palette = "Purples", 
            legend.hist = F,
            title = "(a) OSM") + 
  tm_layout(legend.show = T, 
            title.position = c("center", "bottom"),
            legend.outside= TRUE, 
            frame = F, 
            legend.text.size = 1,  
            legend.position = c("right", "bottom")) +
  tm_shape(london_borough) +
  tm_borders(col = "darkgrey") 
}

#function for gd difference map
gd_diff_map <- function (x) {
  tm_shape(x) +
  tm_raster(style = "fixed", 
            n= 7, 
            title = "OSM - OS (per sq. km)",
            breaks = c(-0.1, -0.005, -0.001,-0.0001,0.0001,0.001, 0.005, 0.1),
            palette = "PiYG", legend.hist = F)+ 
  tm_layout(legend.show = T,
            legend.outside = T,
            frame = F,
            legend.text.size = 1,
            legend.position = c("right", "bottom"),
            title.position = c("center", "bottom"))+
  tm_shape(london_borough)+
  tm_borders(col ="darkgrey")
}

#plot bars
barplot_of_diff <- function(df) {
  ggplot(df,  aes( x = density,
                  fill = density)) +
  geom_bar(color = "grey70")+
  xlab("density (normalised)") +
  ylab("count of categories") +
  scale_fill_brewer(palette = "PiYG", drop=F) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal()+
  theme(axis.text.x= element_blank(),legend.position = "none",
        axis.title.x = element_text( size=13),
        axis.title.y = element_text(size=13))
}

#function to plot local hotspot
local_g_plot <- function(x){
  tm_shape(x)+
  tm_fill( "gstat",
           style = "jenks",
           palette = "-RdYlBu",
           title = "Getis-Ord local GI*",
           labels = c("Hotspot(OS) 99% conf intrv",
                      "Hotspot(OS) 90% conf  intrv",
                      "Not significant",
                      "Hotspot(OSM) 90% conf intrv",
                      "Hotspot(OSM) 99% conf intrv" ))+
  tm_layout(legend.only = F,
            legend.outside= T,
            frame = F, 
            legend.text.size = 1,
            legend.title.size = 2, 
            legend.position= c("left","bottom"))+
  tm_shape(london_borough)+
  tm_borders(col ="darkgrey")
}


#Plot figures for transport poi----------------------------------------
osm_poi_map(transport_osm)
gd_diff_map(t_os_osm)
barplot_of_diff(t_os_osm_df)
local_g_plot(local_g_t)
