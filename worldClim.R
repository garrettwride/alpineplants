## Get Environmental Data
library(geodata)
library(terra)
library(rgbif)
library(sf)
library(dplyr)

rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)
rocky_vect <- vect(rocky_poly)

#Run with resolution 10
bio <- worldclim_global(var = "bio", res = 10, path = "./climate_data/climate_10")

bio_rocky <- crop(bio, rocky_vect)
bio_rocky <- mask(bio_rocky, rocky_vect)
bio_rocky

myExpl <- bio_rocky[[c(1, 3, 4, 12, 15)]]

plot(myExpl)
points(resp.xy, col="red", pch=16)

#Run with resolution 5
bio_5 <- worldclim_global(var = "bio", res = 5, path = "./climate_data/climate_5")

bio_rocky_5 <- crop(bio_5, rocky_vect)
bio_rocky_5 <- mask(bio_rocky_5, rocky_vect)
bio_rocky_5

myExpl_5 <- bio_rocky_5[[c(1, 3, 4, 12, 15)]]

plot(myExpl_5)

#Run with resolution 2.5
bio_2.5 <- worldclim_global(var = "bio", res = 2.5, path = "./climate_data/climate_2.5")

bio_rocky_2.5 <- crop(bio_2.5, rocky_vect)
bio_rocky_2.5 <- mask(bio_rocky_2.5, rocky_vect)
bio_rocky_2.5

myExpl_2.5 <- bio_rocky_2.5[[c(1, 3, 4, 12, 15)]]

plot(myExpl_2.5)