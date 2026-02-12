library(rgbif)
library(sf)
library(dplyr)

hackberry <- occ_search(
  scientificName = "Celtis reticulata",
  hasCoordinate = TRUE,
  limit = 100000
)

# Extract the data frame
hack_data <- hackberry$data

occ_sf <- st_as_sf(
  occ_clean,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)

rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)  # Make sure CRS matches occurrences

occ_rocky <- st_intersection(occ_sf, rocky_poly)
occ_rocky <- occ_sf[st_within(occ_sf, rocky_poly, sparse = FALSE), ]

# occ_rocky <- occ_sf[st_within(occ_sf, rocky_poly, sparse = FALSE), ]
occ_rocky_df <- cbind(
  st_drop_geometry(occ_rocky),
  st_coordinates(occ_rocky)
)

# Rename coordinate columns correctly
colnames(occ_rocky_df)[ncol(occ_rocky_df)-1:0] <- c("decimalLongitude", "decimalLatitude")
# Presence points
resp.var <- rep(1, nrow(occ_rocky_df))

# Coordinates
resp.xy <- occ_rocky_df[, c("decimalLongitude", "decimalLatitude")]

# Name for BIOMOD
resp.name <- "Celtis_reticulata"

