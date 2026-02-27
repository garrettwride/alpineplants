library(terra)
library(dplyr)
library(purrr)
library(rgbif)
library(sf)
library(CoordinateCleaner)
library(geodata)
## Get Rocky mountain polygon
rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)  # Make sure CRS matches occurrences
rocky_wkt <- st_as_text(st_union(rocky_poly))

## Get Environmental Data
bio <- worldclim_global(var = "bio", res = 10, path = "./")

bio_rocky <- crop(bio, rocky_poly)
bio_rocky <- mask(bio_rocky, rocky_poly)

myExpl <- bio_rocky[[c(1, 3, 4, 12, 15)]]

## Get elevation data
elev <- worldclim_global(var = "elev", res = 10, path = "./climate/")

elev_rocky <- crop(elev, rocky_poly)
elev_rocky <- mask(elev_rocky, rocky_poly)

## Get occurrence data
clean_and_get_occurrences <- function(taxon_key, species_name, species_type,
                                      rocky_poly, rocky_wkt) {
  occ <- occ_search(
    taxonKey = taxon_key, 
    hasCoordinate = TRUE,
    hasGeospatialIssue = FALSE, 
    geometry = rocky_wkt, 
    limit = 100000
  )
  
  df <- occ$data
  
  if(nrow(df) == 0) {
    stop("no record found")
  }
  
  cleaned <- clean_coordinates(
    x = df,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    species = "species",
    tests = c("capitals", "centroids", "equal", "gbif", "institutions",
              "outliers", "seas", "zeros")
  )
  
  df_cleaned <- df[cleaned$.summary == TRUE, ]
  
  occ_sf <- st_as_sf(
    df_cleaned,
    coords = c("decimalLongitude", "decimalLatitude"),
    remove = FALSE,
    crs = 4326
  )
  
  occ_rocky <- occ_sf[st_within(occ_sf, rocky_poly, sparse = FALSE), ]
  occ_rocky_df <- st_drop_geometry(occ_rocky)
  
  resp.var <- rep(1, nrow(occ_rocky_df))
  resp.xy <- occ_rocky_df[, c("decimalLongitude", "decimalLatitude")]
  resp.name <- species_name
  
  return(list(
    resp.var = resp.var,
    resp.xy = resp.xy,
    resp.name = resp.name,
    type = species_type
  ))
}

#4,297 occurrences
moss_campion_data <- clean_and_get_occurrences(5384754, "Silene acaulis", "Generalist", rocky_poly, rocky_wkt)
plot(elev_rocky, main = "Generalist Moss")
points(moss_campion_data$resp.xy, col="red", pch=16, cex=0.5)

mountain_dryad_data <- clean_and_get_occurrences(4889932, "Dryas octopetala", "Specialist", rocky_poly, rocky_wkt)
plot(elev_rocky, main = "Specialist Dryad")
points(mountain_dryad_data$resp.xy, col="red", pch=16, cex=0.5)

#6,217; 6,170
maple_data <- clean_and_get_occurrences(3189864, "Acer glabrum", "Generalist", rocky_poly, rocky_wkt)
plot(elev_rocky, main = "Generalist Maple")
points(maple_data$resp.xy, col="red", pch=16, cex=0.5)

#5,209; 5,211
fir_data <- clean_and_get_occurrences(2685313, "Abies lasiocarpa", "Specialist", rocky_poly, rocky_wkt)
plot(elev_rocky, main = "Specialist Fir")
points(fir_data$resp.xy, col="red", pch=16, cex=0.5)

#312; 323
hackberry_data <- clean_and_get_occurrences(6406316, "Celtis reticulata", "Generalist", rocky_poly, rocky_wkt)
plot(elev_rocky, main = "Generalist Hackberry")
points(hackberry_data$resp.xy, col="red", pch=16, cex=0.5)

#469; 518
willow_data <- clean_and_get_occurrences(5372756, "Salix petrophila", "Specialist", rocky_poly, rocky_wkt)
plot(elev_rocky, main = "Specialist Willow")
points(willow_data$resp.xy, col="red", pch=16, cex=0.5)

#for BIOMOD input:
#lmk if this this what garret biomod needed
species_list <- list(
  moss_campion_data,
  mountain_dryad_data,
  maple_data,
  fir_data,
  hackberry_data,
  willow_data
)

calculate_niche_metrics <- function(species_obj, climate_stack, elev_raster) {
  
  pts <- vect(species_obj$resp.xy,
              geom = c("decimalLongitude", "decimalLatitude"),
              crs = "EPSG:4326")
  
  clim_vals <- terra::extract(climate_stack, pts)[, -1]
  elev_vals <- terra::extract(elev_raster, pts)[, -1]
  
  complete <- complete.cases(clim_vals, elev_vals)
  clim_vals <- clim_vals[complete, ]
  elev_vals <- elev_vals[complete]
  
  # Elevational breadth
  elev_range <- max(elev_vals) - min(elev_vals)
  elev_sd    <- sd(elev_vals)
  
  # Multivariate climatic breadth
  clim_scaled <- scale(clim_vals)
  centroid    <- colMeans(clim_scaled)
  distances   <- sqrt(rowSums((clim_scaled - centroid)^2))
  clim_disp   <- mean(distances)
  
  # Geographic convex hull
  pts_sf <- st_as_sf(species_obj$resp.xy,
                     coords = c("decimalLongitude", "decimalLatitude"),
                     crs = 4326)
  
  hull <- st_convex_hull(st_union(pts_sf))
  hull_area <- as.numeric(st_area(st_transform(hull, 5070))) / 10^6
  
  tibble(
    Species = species_obj$resp.name,
    ElevRange = elev_range,
    ElevSD = elev_sd,
    ClimDispersion = clim_disp,
    GeoArea_km2 = hull_area
  )
}

niche_results <- map_dfr(
  species_list,
  calculate_niche_metrics,
  climate_stack = myExpl,
  elev_raster = elev_rocky
)

niche_results <- niche_results %>%
  mutate(
    Elev_z = scale(ElevSD),
    Clim_z = scale(ClimDispersion),
    Geo_z  = scale(GeoArea_km2)
  )

niche_results <- niche_results %>%
  mutate(
    NicheBreadth = Elev_z + Clim_z + Geo_z,
    Specialization = -NicheBreadth
  )