library(tidyverse)
library(biomod2)
library(geodata)
library(terra)
library(Rchelsa)

date <- "1981-2010"

range = c(-115, -103, 35, 54)

clim_vars = c("bio01", "bio12")

bio_clim <- getChelsa(var = "bio1",
                 extent = range,
                 date = date,
                 dataset = "chelsa-bioclim",
                 verbose = true)

bioclim_data <- unwrap(bioclim_current)

rockies_clim_data <- crop(bioclim_data, range)
print(class(rockies_clim_data))
plot(rockies_clim_data)

#using rchelsa didn't work for retreiving the data but either geodata or rchelsa comes witth the bioclim rastor dowloaded in the package
