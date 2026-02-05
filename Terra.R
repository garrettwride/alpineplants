library(terra)
library(geodata)
library(sf)


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
# mapData <- rbind(us, can)
# myBox <- ext(-128, -103, 35, 60)
# regionalMap <- crop(mapData, myBox)

rockyMntPoly <-st_read("./RockyMountainsRegion/rocky_mountains.shp")
rockyMntSpat <- vect(rockyMntPoly)

#plot(regionalMap, col="lightblue", main="US and Canada Map")
plot(rockyMntPoly, col="lightgreen") # add this argument if printing on map , add=TRUE,
