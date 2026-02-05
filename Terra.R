library(terra)
library(geodata)

# How I got the region in the first place (download lable areas from here: https://www.naturalearthdata.com/downloads/50m-physical-vectors/50m-physical-labels/)

#labels50 <- st_read("~/Library/CloudStorage/OneDrive-BrighamYoungUniversity/rstudio/ne_50m_geography_regions_polys")
#unique(labels50$NAME)
#rockies <- labels50[grepl("Rocky Mountains", labels50$NAME, ignore.case=TRUE), ]
#class(rockies)
#plot(st_geometry(rockies), col = "darkgreen", main = "Rocky Mountains")
#st_write(rockies, "~/Desktop/rocky_mountains.shp")

# This is a map of the US and Canada rectangle. (need gadm files to run, not on git bc big)
# us <- gadm(country="USA", level=0, path=".")
# can <- gadm(country="CAN", level=0, path=".")
# map_data <- rbind(us, can)
# my_box <- ext(-128, -103, 35, 60)
# regional_map <- crop(mapData, myBox)

rocky_mnt_poly <-st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_mnt_spat <- vect(rockyMntPoly)

#plot(regional_map, col="lightgrey", main="US and Canada Map")
plot(rocky_mnt_poly, col="darkgrey", max.plot = 1) # add this argument if printing on map , add=TRUE,
