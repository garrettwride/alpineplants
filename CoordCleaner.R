library(CoordinateCleaner)
library(rgbif)
library(dplyr)
library(sf)

#occ_search(taxonKey=5384754, country="US;CA", basisOfRecord="HUMAN_OBSERVATION")
#occ_search(taxonKey=5641386, country="US;CA", basisOfRecord ="HUMAN_OBSERVATION")

#herbs: Silene acaulis (G): 5384754; Primula utahensis (S): 5641386
#trees: Acer glabrum (G): 3189864; Abies lasiocarpa (S):2685313 
#shrubs: Celtis reticulata (G): 6406316; Salix petrophila (S): 5372756

#make sure that cc_inst removes fossils, and can keep all basis of records
#clean_coordinates runs: cc_val, cc_zero, cc_sea, cc_cap, cc_cen, cc_outl

#longitude (X) and latitude (Y), a vector (NA and 1), use coordinates from polygon, then check
#current length: 312
#maybe add basisOfRecord, remove duplicates, hasGeospatialIssue = FALSE

hackberry <- occ_search(
  scientificName = "Celtis reticulata",
  hasCoordinate = TRUE,
  limit = 100000
)

# Extract the data frame
hack_data <- hackberry$data

cleaned <- clean_coordinates(
  x = hack_data,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "species",
  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers", "seas",
            "zeros"),
)

#only keep good points (clean_coords() just evaluates)
#so each row has a geometry point
hack_cleaned <- hack_data[cleaned$.summary == TRUE, ]

#st_as_sf converts foreign object to a spatial object
occ_sf <- st_as_sf(
  hack_cleaned,
  #occ_clean, #can't find this function (placeholder?)
  coords = c("decimalLongitude", "decimalLatitude"),
  remove = FALSE, #do not remove coords from df
  crs = 4326
)

#reads .shp into sf object
rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)  # Make sure CRS matches occurrences

#only keep points inside polygon (only rocky mountains)
occ_rocky <- occ_sf[st_within(occ_sf, rocky_poly, sparse = FALSE), ]

#remove spatial geo col, back to normal df
occ_rocky_df <- st_drop_geometry(occ_rocky)


#BIOMOD inputs:

# Presence points
resp.var <- rep(1, nrow(occ_rocky_df))

# Coordinates
resp.xy <- occ_rocky_df[, c("decimalLongitude", "decimalLatitude")]

# Name for BIOMOD
resp.name <- "Celtis_reticulata"
