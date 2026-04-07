library(terra)
library(geodata)
library(sf)

# How I got the region in the first place (download label areas from here: https://www.naturalearthdata.com/downloads/50m-physical-vectors/50m-physical-labels/)

# Replace the file path with where the downloaded file is on your personal device
labels_download_file_path = "./data/ne_50m_geography_regions_polys"

# Replace the file path with the location of the folder you created on your personal device
final_shape_folder_path = "./data/RockyMountainsRegion"


labels50 <- st_read(labels_download_file_path)
#unique(labels50$NAME)
rockies <- labels50[grepl("Rocky Mountains", labels50$NAME, ignore.case=TRUE), ]
class(rockies)

# Check that you got the right shape (should look kinda like a dinosaur/rotting guitar)
plot(st_geometry(rockies), col = "cyan", main = "Rocky Mountains")

# Save shape in folder
new_rocky_shp_file_path = paste0(final_shape_folder_path,"/rocky_mountains.shp")

if (file.exists(new_rocky_shp_file_path)) {
  message("Shapefile already exists!")
} else {
  st_write(rockies, new_rocky_shp_file_path)
}

# Check that the file can be retrieved 
if (file.access("./RockyMountainsRegion/rocky_mountains.shp", mode = 4) == 0) {
  message("The file CAN be opened.")
} else {
  message("There was an error saving the file.") 
}