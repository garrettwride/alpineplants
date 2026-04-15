library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(biomod2)
library(geodata)

species_df <- readRDS("data/species_occurrences.rds")
myExpl     <- readRDS("data/myExpl.rds")

us_states <- gadm(country="USA", level=1, path=tempdir())
can_provinces <- gadm(country="CAN", level=1, path=tempdir())
admin_lines <- rbind(us_states, can_provinces)
admin_lines <- project(admin_lines, crs(myExpl))

#1
png("Campanula_Projection_Big1.png", width = 10, height = 15, units = "in", res = 300)

par(mar = c(10, 8, 12, 8), bg = "white")

sp <- species_df %>% filter(Species == "Campanula rotundifolia")
respXY <- sp$Occurrences[[1]][, c("decimalLongitude", "decimalLatitude")]
resp   <- rep(1, nrow(respXY))

myBiomodData <- BIOMOD_FormatingData(
  resp.var = resp,
  expl.var = myExpl,
  resp.xy  = respXY,
  resp.name = "Campanula.rotundifolia",
  PA.nb.rep = 1,
  PA.nb.absences = length(resp),
  PA.strategy = "random",
  PA.dist.min = 50000,
  PA.dist.max = 150000
)

myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  models = "RF",
  CV.strategy = "random",
  CV.nb.rep = 1,
  CV.perc = 0.8,
  metric.eval = c("TSS")
)

myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = "current",
  new.env = myExpl,
  models.chosen = "all"
)

# Plotting with your requested style
r <- rast(myBiomodProj@proj.out@link[1])
get_model_layer <- function(r, algo) {
  pattern <- paste0("allData_allRun_", algo)
  lyr <- grep(pattern, names(r), value = TRUE)
  r[[lyr[1]]]
}

plot(get_model_layer(r, "RF"), 
     main = "", 
     xlab = "Longitude", ylab = "Latitude", 
     cex.lab = 1.5,
     plg = list(title = "Suitability", title.cex = 1.2, title.font = 2))

plot(admin_lines, add = TRUE, border = "gray40", lwd = 0.5)

mtext("Campanula rotundifolia\nRF projection", side = 3, line = 3, cex = 2, font = 2)

#2
png("Delphinium_Projection_Big1.png", width = 10, height = 15, units = "in", res = 300)

par(mar = c(10, 8, 12, 8), bg = "white")

sp <- species_df %>% filter(Species == "Delphinium occidentale")
respXY <- sp$Occurrences[[1]][, c("decimalLongitude", "decimalLatitude")]
resp   <- rep(1, nrow(respXY))

myBiomodData <- BIOMOD_FormatingData(
  resp.var = resp,
  expl.var = myExpl,
  resp.xy  = respXY,
  resp.name = "Delphinium occidentale",
  PA.nb.rep = 1,
  PA.nb.absences = length(resp),
  PA.strategy = "random",
  PA.dist.min = 50000,
  PA.dist.max = 150000
)

myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  models = "RF",
  CV.strategy = "random",
  CV.nb.rep = 1,
  CV.perc = 0.8,
  metric.eval = c("TSS")
)

myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = "current",
  new.env = myExpl,
  models.chosen = "all"
)

# Plotting with your requested style
r <- rast(myBiomodProj@proj.out@link[1])
get_model_layer <- function(r, algo) {
  pattern <- paste0("allData_allRun_", algo)
  lyr <- grep(pattern, names(r), value = TRUE)
  r[[lyr[1]]]
}

plot(get_model_layer(r, "RF"), 
     main = "", 
     xlab = "Longitude", ylab = "Latitude", 
     cex.lab = 1.5,
     plg = list(title = "Suitability", title.cex = 1.2, title.font = 2))

plot(admin_lines, add = TRUE, border = "gray40", lwd = 0.5)

mtext("Delphinium occidentale\nRF projection", side = 3, line = 3, cex = 2, font = 2)
dev.off()

#3
png("Campanula_Projection_Big2.png", width = 10, height = 15, units = "in", res = 300)

par(mar = c(10, 8, 12, 8), bg = "white")

sp <- species_df %>% filter(Species == "Campanula rotundifolia")
respXY <- sp$Occurrences[[1]][, c("decimalLongitude", "decimalLatitude")]
resp   <- rep(1, nrow(respXY))

myBiomodData <- BIOMOD_FormatingData(
  resp.var = resp,
  expl.var = myExpl,
  resp.xy  = respXY,
  resp.name = "Campanula.rotundifolia",
  PA.nb.rep = 1,
  PA.nb.absences = length(resp),
  PA.strategy = "random",
  PA.dist.min = 50000,
  PA.dist.max = 150000
)

myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  models = "GLM",
  CV.strategy = "random",
  CV.nb.rep = 1,
  CV.perc = 0.8,
  metric.eval = c("TSS")
)

myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = "current",
  new.env = myExpl,
  models.chosen = "all"
)

# Plotting with your requested style
r <- rast(myBiomodProj@proj.out@link[1])
get_model_layer <- function(r, algo) {
  pattern <- paste0("allData_allRun_", algo)
  lyr <- grep(pattern, names(r), value = TRUE)
  r[[lyr[1]]]
}

plot(get_model_layer(r, "GLM"), 
     main = "", 
     xlab = "Longitude", ylab = "Latitude", 
     cex.lab = 1.5,
     plg = list(title = "Suitability", title.cex = 1.2, title.font = 2))

plot(admin_lines, add = TRUE, border = "gray40", lwd = 0.5)

mtext("Campanula rotundifolia\nGLM projection", side = 3, line = 3, cex = 2, font = 2)

#4
png("Delphinium_Projection_Big2.png", width = 10, height = 15, units = "in", res = 300)

par(mar = c(10, 8, 12, 8), bg = "white")

sp <- species_df %>% filter(Species == "Delphinium occidentale")
respXY <- sp$Occurrences[[1]][, c("decimalLongitude", "decimalLatitude")]
resp   <- rep(1, nrow(respXY))

myBiomodData <- BIOMOD_FormatingData(
  resp.var = resp,
  expl.var = myExpl,
  resp.xy  = respXY,
  resp.name = "Delphinium occidentale",
  PA.nb.rep = 1,
  PA.nb.absences = length(resp),
  PA.strategy = "random",
  PA.dist.min = 50000,
  PA.dist.max = 150000
)

myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData,
  models = "GLM",
  CV.strategy = "random",
  CV.nb.rep = 1,
  CV.perc = 0.8,
  metric.eval = c("TSS")
)

myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = "current",
  new.env = myExpl,
  models.chosen = "all"
)

# Plotting with your requested style
r <- rast(myBiomodProj@proj.out@link[1])
get_model_layer <- function(r, algo) {
  pattern <- paste0("allData_allRun_", algo)
  lyr <- grep(pattern, names(r), value = TRUE)
  r[[lyr[1]]]
}

plot(get_model_layer(r, "GLM"), 
     main = "", 
     xlab = "Longitude", ylab = "Latitude", 
     cex.lab = 1.5,
     plg = list(title = "Suitability", title.cex = 1.2, title.font = 2))

plot(admin_lines, add = TRUE, border = "gray40", lwd = 0.5)

mtext("Delphinium occidentale\nGLM projection", side = 3, line = 3, cex = 2, font = 2)
dev.off()