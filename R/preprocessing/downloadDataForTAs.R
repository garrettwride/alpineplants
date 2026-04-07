library(rgbif)
library(sf)
library(dplyr)
library(CoordinateCleaner)
library(purrr)
library(tidyr)
library(terra)
library(geodata)

dir.create("data", showWarnings = FALSE)
dir.create("climate", showWarnings = FALSE)

## Rocky Mountain polygon
rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)
rocky_wkt  <- st_as_text(st_union(rocky_poly))

## Species list
species_info <- tibble(
  Species = c(
    "Silene acaulis",
    "Delphinium occidentale",
    "Campanula rotundifolia"
  ),
  taxon_key = c(
    3033713,
    5410907
  )
)

clean_and_get_occurrences <- function(taxon_key, species_name) {
  
  message("Downloading: ", species_name)
  
  occ <- occ_search(
    taxonKey = taxon_key,
    hasCoordinate = TRUE,
    hasGeospatialIssue = FALSE,
    geometry = rocky_wkt,
    limit = 100000
  )
  
  df <- occ$data
  
  if(nrow(df) == 0) {
    warning("No records for ", species_name)
    return(NULL)
  }
  
  cleaned <- clean_coordinates(
    x = df,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    species = "species",
    tests = c("capitals", "centroids", "equal", "gbif",
              "institutions","outliers","seas","zeros")
  )
  
  df_cleaned <- df[cleaned$.summary == TRUE, ]
  
  occ_sf <- st_as_sf(
    df_cleaned,
    coords = c("decimalLongitude", "decimalLatitude"),
    remove = FALSE,
    crs = 4326
  )
  
  occ_rocky <- occ_sf[st_within(occ_sf, rocky_poly, sparse = FALSE), ]
  occ_df <- st_drop_geometry(occ_rocky)
  
  return(occ_df)
}

species_occurrences <- species_info %>%
  mutate(
    Occurrences = map2(taxon_key, Species, clean_and_get_occurrences)
  )

saveRDS(species_occurrences, "data/species_occurrences.rds")

message("Species download complete.")

message("Downloading WorldClim climate data...")

bio <- worldclim_global(var = "bio", res = 10, path = "climate/")

rocky_vect <- vect(rocky_poly)

bio_rocky <- crop(bio, rocky_vect)
bio_rocky <- mask(bio_rocky, rocky_vect)

# Select subset of bioclim variables
myExpl <- bio_rocky[[c(1, 3, 4, 12, 15)]]

saveRDS(myExpl, "data/myExpl.rds")

message("Downloading elevation data...")

elev <- worldclim_global(var = "elev", res = 2.5, path = "climate/")

elev_rocky <- crop(elev, rocky_vect)
elev_rocky <- mask(elev_rocky, rocky_vect)

saveRDS(elev_rocky, "data/elev_rocky.rds")

message("Climate and elevation saved.")
message("Download script complete.")