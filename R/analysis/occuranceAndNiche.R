library(terra)
library(rgbif)
library(sf)
library(dplyr)
library(CoordinateCleaner)
library(tidyr)
library(purrr)
library(geodata)
library(gridExtra)


## Get Rocky Mountain Polygon
rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)  # Make sure CRS matches occurrences
rocky_wkt <- st_as_text(st_union(rocky_poly))

## Get Environmental Data
bio <- worldclim_global(var = "bio", res = 5, path = "./")

bio_rocky <- crop(bio, rocky_poly)
bio_rocky <- mask(bio_rocky, rocky_poly)

myExpl <- bio_rocky[[c(1, 3, 4, 12, 15)]]

# Get elevation data
elev <- worldclim_global(var = "elev", res = 5, path = "./climate/")

elev_rocky <- crop(elev, rocky_poly)
elev_rocky <- mask(elev_rocky, rocky_poly)

## Get occurrence Data
species_info <- tibble(
  Species = c(
    "Silene acaulis",
    "Delphinium occidentale",
    "Dryas octopetala",
    "Acer glabrum",
    "Abies lasiocarpa",
    "Sorbus scopulina",
    "Dasiphora fruticosa",
    "Celtis reticulata",
    "Salix petrophila",
    "Carex spectabilis",
    "Carex arapahoensis",
    "Carex perglobosa",
    "Juncus drummondii",
    "Pseudoroegneria spicata",
    "Bistorta vivipara",
    "Campanula rotundifolia",
    "Artemisia tridentata",
    "Amelanchier utahensis",
    "Pinus edulis",
    "Pinus longaeva"
  ),
  taxon_key = c(
    5384754,
    3033713,
    4889932,
    3189864,
    2685313,
    3012298,
    5370380,
    6406316,
    5372756,
    2723145,
    2722910,
    2723300,
    2701717,
    2705699,
    2889299,
    5410907,
    9396703,
    3023964,
    5285796,
    5285258
  )
)

clean_and_get_occurrences <- function(taxon_key, species_name, rocky_poly, rocky_wkt) {
  
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
    resp.name = resp.name
  ))
}

species_df <- species_info %>%
  mutate(
    Data = map2(
      taxon_key,
      Species,
      ~ clean_and_get_occurrences(
        taxon_key = .x,
        species_name = .y,
        rocky_poly = rocky_poly,
        rocky_wkt = rocky_wkt
      )
    )
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
    ElevRange = elev_range,
    ElevSD = elev_sd,
    ClimDispersion = clim_disp,
    GeoArea_km2 = hull_area
  )
}

niche_metrics <- species_df %>%
  mutate(
    metrics = map(
      Data,
      calculate_niche_metrics,
      climate_stack = myExpl,
      elev_raster = elev_rocky
    )
  ) %>%
  unnest(metrics)

niche_metrics <- niche_metrics %>%
  rowwise() %>%
  mutate(
    Abundance = length(Data$resp.xy$decimalLongitude)
  ) %>%
  ungroup()

niche_metrics <- niche_metrics %>%
  mutate(
    logAbundance = log10(Abundance)
  ) %>%
  ungroup()

niche_metrics <- niche_metrics %>%
  mutate(
    Elev_z = as.numeric(scale(ElevSD)),
    Clim_z = as.numeric(scale(ClimDispersion)),
    Geo_z  = as.numeric(scale(GeoArea_km2)),
    NicheBreadth = Elev_z + Clim_z + Geo_z,
    Specialization = -NicheBreadth
  )

Vi <- niche_metrics %>%
  mutate(
    Elev_z = scale(ElevSD),
    Clim_z = scale(ClimDispersion),
    Geo_z  = scale(GeoArea_km2)
  ) %>%
  mutate(
    NicheBreadth = Elev_z + Clim_z + Geo_z,
    Specialization = -NicheBreadth
  )

correction_model <- lm(NicheBreadth ~ log10(Abundance), data = niche_results)
niche_results$NicheBreadth_Corrected <- residuals(correction_model)

niche_results <- niche_results %>%
  mutate(
    Niche_Final_z = as.numeric(scale(NicheBreadth_Corrected)),
    Specialization_Rank = -Niche_Final_z
  )

niche_results <- niche_results %>%
    select(Species, Niche_Final_z)

shapiro.test(niche_results$Niche_Final_z)
shapiro.test(niche_results$NicheBreadth)
shapiro.test(niche_results$Abundance)

cor.test(niche_results$Abundance, niche_results$NicheBreadth, method = "pearson")
cor.test(niche_results$Abundance, niche_results$Niche_Final_z, method = "pearson")


cor.test(niche_results$Abundance, niche_results$Elev_z, method = "pearson")
cor.test(niche_results$Abundance, niche_results$Clim_z, method = "pearson")
cor.test(niche_results$Abundance, niche_results$Geo_z, method = "pearson")