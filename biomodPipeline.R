library(terra)
library(rgbif)
library(sf)
library(dplyr)

## Get Occurrence Data
hackberry <- occ_search(
  scientificName = "Celtis reticulata",
  hasCoordinate = TRUE,
  limit = 100000
)

# Extract the data frame
hack_data <- hackberry$data

occ_sf <- st_as_sf(
  hack_data,
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
resp.name <- "celtis_reticulata"

## Get Environmental Data
library(geodata)

bio <- worldclim_global(var = "bio", res = 10, path = "./climate_data/climate_10")

bio_rocky <- crop(bio, rocky_vect)
bio_rocky <- mask(bio_rocky, rocky_vect)
bio_rocky

myExpl <- bio_rocky[[c(1, 3, 4, 12, 15)]]

# Get elevation data
elev <- worldclim_global(var = "elev", res = 10, path = "./climate/")

elev_rocky <- crop(elev, rocky_vect)
elev_rocky <- mask(elev_rocky, rocky_vect)


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
n_boot        <- 10

results <- data.frame(
  Species        = character(),
  SpeciesType    = character(),
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
# for (species in species_list)
respName <- resp.name
resp <- resp.var
respXY <- resp.xy
SpeciesType <- "generalist" # Generalist/Specialist

# Baseline ensemble model that uses all occurrence data

baseline_ensemble <- normalize_raster(baseline_r)


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

  boot_occ <- occ_rocky_df[
    sample(seq_len(nrow(occ_rocky_df)), replace = TRUE),
  ]
  
  boot_resp <- rep(1, nrow(boot_occ))
  boot_xy   <- boot_occ[, c("decimalLongitude", "decimalLatitude")]
  
  # Set up Pseudo-absences (next 2 loops)
  for (pa_strat in pa_strategies) {
    for (pa_mult in pa_numbers) {
      PA_n <- pa_mult * length(resp)
      
      myBiomodData <- BIOMOD_FormatingData(
        dir.name = tempdir(),
        resp.var = boot_resp,
        expl.var = myExpl,
        resp.xy  = respXY,
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
        SpeciesType = SpeciesType,
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
        SpeciesType = SpeciesType,
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
        SpeciesType = SpeciesType,
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

write.csv(results, "./testFiles/SDM_results_for_LMM.csv", row.names = FALSE)
