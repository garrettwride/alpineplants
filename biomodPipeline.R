library(terra)
library(rgbif)
library(sf)
library(dplyr)
library(CoordinateCleaner)
library(tidyr)
library(purrr)

## Get Rocky Mountain Polygon
rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)  # Make sure CRS matches occurrences
rocky_wkt <- st_as_text(st_union(rocky_poly))

## Get Environmental Data
library(geodata)

bio <- worldclim_global(var = "bio", res = 10, path = "./")

bio_rocky <- crop(bio, rocky_poly)
bio_rocky <- mask(bio_rocky, rocky_poly)

myExpl <- bio_rocky[[c(1, 3, 4, 12, 15)]]

# Get elevation data
elev <- worldclim_global(var = "elev", res = 10, path = "./climate/")

elev_rocky <- crop(elev, rocky_poly)
elev_rocky <- mask(elev_rocky, rocky_poly)

## Get occurrence Data
species_info <- tibble(
  Species = c(
    "Silene acaulis",
    "Dryas octopetala",
    "Acer glabrum",
    "Abies lasiocarpa",
    "Celtis reticulata",
    "Salix petrophila"
  ),
  taxon_key = c(
    5384754,
    4889932,
    3189864,
    2685313,
    6406316,
    5372756
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

print(species_df)
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
  mutate(
    Elev_z = as.numeric(scale(ElevSD)),
    Clim_z = as.numeric(scale(ClimDispersion)),
    Geo_z  = as.numeric(scale(GeoArea_km2)),
    NicheBreadth = Elev_z + Clim_z + Geo_z,
    Specialization = -NicheBreadth
  )

niche_results <- niche_results %>%
  mutate(
    Elev_z = scale(ElevSD),
    Clim_z = scale(ClimDispersion),
    Geo_z  = scale(GeoArea_km2)
  ) %>%
  mutate(
    NicheBreadth = Elev_z + Clim_z + Geo_z,
    Specialization = -NicheBreadth
  )

species_df <- species_df %>%
  left_join(
    niche_metrics %>% select(Species, NicheBreadth),
    by = "Species"
  )
### Biomod pipeline
library(biomod2)

## Helper functions:
# Extract model layer safely by algorithm name
get_model_layer <- function(r, algo, type = "allData_allRun") {
  
  pattern <- paste0(type, "_", algo, "$")
  lyr <- grep(pattern, names(r), value = TRUE)
  
  if (length(lyr) != 1) {
    stop(paste("Expected 1 layer for", algo, "but found", length(lyr)))
  }
  
  return(r[[lyr]])
}

# Normalize raster to sum to 1
normalize_raster <- function(r) {
  r[r < 0] <- 0
  total <- global(r, "sum", na.rm = TRUE)[1,1]
  r / total
}

# Calculate Schoener's D
schoeners_d <- function(p, q) {
  1 - 0.5 * global(abs(p - q), "sum", na.rm = TRUE)[1,1]
}

# Calculate elevation Centroid
elevation_centroid <- function(suitability, elevation) {
  suitability[suitability < 0] <- 0
  
  weighted_sum <- global(suitability * elevation, "sum", na.rm = TRUE)[1,1]
  suitability_sum <- global(suitability, "sum", na.rm = TRUE)[1,1]
  
  weighted_sum / suitability_sum
}


algorithms   <- c("GLM", "GAM", "RF") # "MAXENT"
pa_strategies <- c("random", "disk")
pa_numbers    <- c(1, 3, 5)   # × presences
n_boot        <- 5

results <- data.frame(
  Species        = character(),
  NicheBreadth    = numeric(),
  Bootstrap      = integer(),
  Algorithm      = character(),
  PAStrategy     = character(),
  PANumber       = integer(),
  MeanSuitability= numeric(),
  SuitableArea   = numeric(),
  SchoenersD     = numeric(),
  ElevationCentroid = numeric(),
  stringsAsFactors = FALSE
)

# outermost loop
for (i in seq_len(nrow(species_df))) {
  respName <- species_df$Species[i]
  species  <- species_df$Data[[i]]
  niche_breadth <- species_df$NicheBreadth[i]
  resp <- species$resp.var
  respXY <- species$resp.xy
  
  # Baseline ensemble model that uses all occurrence data
  baseline_data <- BIOMOD_FormatingData(
    dir.name = tempdir(),
    resp.var = resp,
    expl.var = myExpl,
    resp.xy  = respXY,
    resp.name = paste(respName, "baseline", sep = "_"),
    PA.nb.rep = 1,
    PA.nb.absences = length(resp),
    PA.strategy = "disk",
    PA.dist.min = 50000,
    PA.dist.max = 150000
  )
  
  baseline_model <- BIOMOD_Modeling(
    bm.format = baseline_data,
    models = algorithms,
    CV.strategy = "random",
    CV.nb.rep = 1,
    CV.perc = 0.8,
    CV.do.full.models = TRUE,
    metric.eval = "TSS"
  )
  
  baseline_ensemble_model <- BIOMOD_EnsembleModeling(
    bm.mod = baseline_model,
    models.chosen = "all",
    em.by = "all",
    em.algo = "EMwmean",
    metric.select = "TSS",
    metric.select.thresh = 0.7,
    metric.eval = "TSS"
  )
  
  baseline_ensemble_proj <- BIOMOD_EnsembleForecasting(
    bm.em = baseline_ensemble_model,
    proj.name = "BaselineEnsemble",
    new.env = myExpl
  )
  
  baseline_r <- rast(baseline_ensemble_proj@proj.out@link[1])
  baseline_ensemble <- normalize_raster(baseline_r)
  
  # Bootsrap occurance data (next loop)
  for (b in seq_len(n_boot)) {
    set.seed(b)
  
    boot_idx  <- sample(seq_len(nrow(respXY)), replace = TRUE)
    boot_xy   <- respXY[boot_idx, ]
    boot_resp <- rep(1, nrow(boot_xy))
    
    # Set up Pseudo-absences (next 2 loops)
    for (pa_strat in pa_strategies) {
      for (pa_mult in pa_numbers) {
        PA_n <- pa_mult * length(resp)
        
        myBiomodData <- BIOMOD_FormatingData(
          dir.name = tempdir(),
          resp.var = boot_resp,
          expl.var = myExpl,
          resp.xy  = boot_xy,
          resp.name = paste(
            respName, "boot", b, pa_strat, pa_mult, sep = "_"),
          PA.nb.rep = 1,
          PA.nb.absences = PA_n,
          PA.strategy = pa_strat,
          PA.dist.min = 50000,
          PA.dist.max = 150000
        )
        
        myBiomodModelOut <- BIOMOD_Modeling(
          bm.format = myBiomodData,
          models = algorithms,
          CV.strategy = "random",
          CV.nb.rep = 1,
          CV.perc = 0.8,
          CV.do.full.models = TRUE,
          metric.eval = c("TSS")
        )
        
        myBiomodProj <- BIOMOD_Projection(
          bm.mod = myBiomodModelOut,
          proj.name = paste0("Current_boot", b),
          new.env = myExpl, models.chosen = "all",
          metric.binary = "TSS",
          metric.filter = "all",
          build.clamping.mask = TRUE 
        )
        
        # Extract the SpatRaster object
        r <- rast(myBiomodProj@proj.out@link[1])
        r_bin <- rast(myBiomodProj@proj.out@link[2])
        # Extract the layers we want
        glm_raster <- get_model_layer(r, "GLM")
        gam_raster <- get_model_layer(r, "GAM")
        rf_raster  <- get_model_layer(r, "RF")
  
        glm_bin <- get_model_layer(r_bin, "GLM")
        gam_bin <- get_model_layer(r_bin, "GAM")
        rf_bin  <- get_model_layer(r_bin, "RF")
        
        # Calculate Metrics
        mean_glm <- global(glm_raster, "mean", na.rm = TRUE)[1,1]
        area_glm <- global(glm_bin, "sum", na.rm = TRUE)[1,1] * prod(res(glm_bin))
        mean_gam <- global(gam_raster, "mean", na.rm = TRUE)[1,1]
        area_gam <- global(gam_bin, "sum", na.rm = TRUE)[1,1] * prod(res(gam_bin))
        mean_rf <- global(rf_raster, "mean", na.rm = TRUE)[1,1]
        area_rf <- global(rf_bin, "sum", na.rm = TRUE)[1,1] * prod(res(rf_bin))
        
        # Normalize rasters for Schoener's D
        glm_norm <- normalize_raster(glm_raster)
        gam_norm <- normalize_raster(gam_raster)
        rf_norm  <- normalize_raster(rf_raster)
        D_glm <- schoeners_d(glm_norm, baseline_ensemble)
        D_gam <- schoeners_d(gam_norm, baseline_ensemble)
        D_rf  <- schoeners_d(rf_norm, baseline_ensemble)
        
        # Calculate Elevation Centroid
        elev_glm <- elevation_centroid(glm_raster, elev_rocky)
        elev_gam <- elevation_centroid(gam_raster, elev_rocky)
        elev_rf  <- elevation_centroid(rf_raster, elev_rocky)
        
        results <- rbind(results, data.frame(
          Species = respName,
          NicheBreadth = niche_breadth,
          Bootstrap = b,
          Algorithm = "GLM",
          PAStrategy = pa_strat,
          PANumber = pa_mult,
          MeanSuitability = mean_glm,
          SuitableArea = area_glm,
          SchoenersD = D_glm,
          ElevationCentroid = elev_glm
        ))
        results <- rbind(results, data.frame(
          Species = respName,
          NicheBreadth = niche_breadth,
          Bootstrap = b,
          Algorithm = "GAM",
          PAStrategy = pa_strat,
          PANumber = pa_mult,
          MeanSuitability = mean_gam,
          SuitableArea = area_gam,
          SchoenersD = D_gam,
          ElevationCentroid = elev_gam
        ))
        results <- rbind(results, data.frame(
          Species = respName,
          NicheBreadth = niche_breadth,
          Bootstrap = b,
          Algorithm = "RF",
          PAStrategy = pa_strat,
          PANumber = pa_mult,
          MeanSuitability = mean_rf,
          SuitableArea = area_rf,
          SchoenersD = D_rf,
          ElevationCentroid = elev_glm
        ))
      }
    }
  }
}

write.csv(results, "./testFiles/SDM_results_for_LMM.csv", row.names = FALSE)
