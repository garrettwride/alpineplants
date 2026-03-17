library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(biomod2)
install.packages("R.utils")
install.packages("randomForest")

#expl = environment data
species_df <- readRDS("data/species_occurrences.rds")
myExpl     <- readRDS("data/myExpl.rds")

#choose one species, pull out lon and lat occurrences, create response var
sp <- species_df %>% 
  filter(Species == "Artemisia tridentata")
respXY <- sp$Occurrences[[1]][, c("decimalLongitude", "decimalLatitude")]
resp   <- rep(1, nrow(respXY))

#preps for modeling
myBiomodData <- BIOMOD_FormatingData(
  resp.var = resp,
  expl.var = myExpl,
  resp.xy  = respXY,
  resp.name = as.character(sp$Species[1]),
  #biomod creates fake absences
  PA.nb.rep = 1,
  PA.nb.absences = length(resp),
  #disk = sampled in a ring around the presences, 50-150km
  PA.strategy = "disk",
  PA.dist.min = 50000,
  PA.dist.max = 150000#,
  #filter.raster = TRUE #if automatic filtering
)

algorithms <- c("GLM", "GAM", "RF", "MAXENT")

#fit models using data
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  models = algorithms,
  CV.strategy = "random",
  CV.nb.rep = 1,
  #80% training, 20% testing
  CV.perc = 0.8,
  metric.eval = c("TSS")
)
#*artemisia tridentata maxent failed, no executable file

#apply trained models to environmental space and produce map of habitat suitability
myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = "current",
  new.env = myExpl,
  models.chosen = "all"
)

#read projection file into raster object
r <- rast(myBiomodProj@proj.out@link[1])

# Extract model layer safely by algorithm name
#find raster layer for each algorithm
get_model_layer <- function(r, algo) {
  pattern <- paste0("allData_allRun_", algo)
  lyr <- grep(pattern, names(r), value = TRUE)
  r[[lyr[1]]]
}

#print(names(r))

par(mfrow = c(2,2))

for (algo in algorithms) {
  plot(get_model_layer(r, algo), main = algo)
}
