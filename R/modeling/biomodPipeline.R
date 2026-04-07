library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)

## Get Rocky Mountain Polygon
rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)  # Make sure CRS matches occurrences
rocky_wkt <- st_as_text(st_union(rocky_poly))
species_df <- readRDS("./data/species_occurrences.rds")
myExpl     <- readRDS("./data/myExpl.rds")
elev_rocky <- readRDS("./data/elev_rocky.rds")

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
    metrics = map(Occurrences, function(df) {
      species_obj <- list(
        resp.xy = df[, c("decimalLongitude","decimalLatitude")]
      )
      
      calculate_niche_metrics(
        species_obj,
        climate_stack = myExpl,
        elev_raster = elev_rocky
      )
    })
  ) %>%
  unnest(metrics) %>%
  mutate(
    Elev_z = as.numeric(scale(ElevSD)),
    Clim_z = as.numeric(scale(ClimDispersion)),
    Geo_z  = as.numeric(scale(GeoArea_km2)),
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


algorithms   <- c("GLM", "GAM", "RF", "MAXENT")
pa_strategies <- c("random", "disk")
pa_numbers    <- c(1, 3, 5)   # × presences
n_boot        <- 8

results <- data.frame(
  Species = character(),
  NicheBreadth = numeric(),
  Occurrences = numeric(),
  LogOccurrences = numeric(),
  Bootstrap = integer(),
  Algorithm = character(),
  PAStrategy = character(),
  PANumber = integer(),
  MeanSuitability = numeric(),
  SuitableArea = numeric(),
  SchoenersD = numeric(),
  ElevationCentroid = numeric(),
  TSS = numeric(),
  stringsAsFactors = FALSE
)

# outermost loop
for (i in seq_len(nrow(species_df))) {
  
  respName <- species_df$Species[i]
  niche_breadth <- species_df$NicheBreadth[i]
  
  occ_df <- species_df$Occurrences[[i]]
  
  respXY <- occ_df[, c("decimalLongitude", "decimalLatitude")]
  resp   <- rep(1, nrow(respXY))
  
  num_occ <- nrow(respXY)
  log_occ <- log10(num_occ)
  
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
        PA_n <- pa_mult * length(boot_resp)
        
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
        # Extract the SpatRaster objects
        r <- rast(myBiomodProj@proj.out@link[1])
        r_bin <- rast(myBiomodProj@proj.out@link[2])
        
        # Loop over algorithms
        for (algo in algorithms) {
          
          # Extract raster layers
          algo_raster <- get_model_layer(r, algo)
          algo_bin    <- get_model_layer(r_bin, algo)
          
          # ---- Metrics ----
          mean_suit <- global(algo_raster, "mean", na.rm = TRUE)[1,1]
          
          suitable_area <- global(algo_bin, "sum", na.rm = TRUE)[1,1] *
            prod(res(algo_bin))
          
          # Schoener's D
          algo_norm <- normalize_raster(algo_raster)
          D_val <- schoeners_d(algo_norm, baseline_ensemble)
          
          # Elevation centroid
          elev_cent <- elevation_centroid(algo_raster, elev_rocky)
          
          # Model Performance (TSS)
          evals_df <- as.data.frame(get_evaluations(myBiomodModelOut))
          
          tss_row <- evals_df[
            evals_df$algo == algo &
              evals_df$PA == "allData" &
              evals_df$run == "allRun" &
              evals_df$metric.eval == "TSS",
          ]
          
          if(nrow(tss_row) == 0){
            tss_val <- NA
          } else {
            tss_val <- tss_row$cutoff
          }
          
          # Save results
          results <- rbind(
            results,
            data.frame(
              Species = respName,
              NicheBreadth = niche_breadth,
              Occurrences = num_occ,
              LogOccurrences = log_occ,
              Bootstrap = b,
              Algorithm = algo,
              PAStrategy = pa_strat,
              PANumber = pa_mult,
              MeanSuitability = mean_suit,
              SuitableArea = suitable_area,
              SchoenersD = D_val,
              ElevationCentroid = elev_cent,
              TSS = tss_val
            )
          )
        }
      }
    }
  }
}
print(results)
write.csv(results, "./data/SDM_results_for_LMM.csv", row.names = FALSE)
