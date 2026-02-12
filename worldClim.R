## Get Environmental Data
library(geodata)
library(terra)
library(rgbif)
library(sf)
library(dplyr)

rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)

bio <- worldclim_global(var = "bio", res = 10, path = "./")
rocky_vect <- vect(rocky_poly)

bio_rocky <- crop(bio, rocky_vect)
bio_rocky <- mask(bio_rocky, rocky_vect)
bio_rocky

myExpl <- bio_rocky[[c(1, 3, 4, 12, 15)]]

plot(myExpl[[1]])
points(resp.xy, col="red", pch=16)