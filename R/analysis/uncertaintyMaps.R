library(terra)
library(dplyr)

# Choose species
# species_name <- "Delphinium occidentale"
species_name <- "Campanula rotundifolia"  
species_dir <- file.path("projections", species_name)

# Get all suitability rasters
files <- list.files(
  species_dir,
  pattern = "_suitability\\.tif$",
  full.names = TRUE
)
# Check
length(files)

# Stack rasters
rast_stack <- rast(files)

# Compute pixel-wise SD
sd_map <- app(rast_stack, sd, na.rm = TRUE)

writeRaster(
  sd_map,
  filename = paste0(species_dir, "/", species_name, "_SD_total.tif"),
  overwrite = TRUE
)

## Methodological Uncertainty
# Example: fix bootstrap 1
subset_files <- files[grepl("boot_1_", files)]
length(files)

rast_stack <- rast(subset_files)

sd_method <- app(rast_stack, sd, na.rm = TRUE)

writeRaster(
  sd_method,
  filename = paste0(species_dir, "/", species_name, "_SD_method_uncertainty.tif"),
  overwrite = TRUE
)
## Stochastic Uncertainty (Bootstraps)
# Example: fix algorithm + PA settings
subset_files <- files[grepl("random_PA1_GLM", files)]
length(files)

rast_stack <- rast(subset_files)

sd_stochastic <- app(rast_stack, sd, na.rm = TRUE)

writeRaster(
  sd_stochastic,
  filename = paste0(species_dir, "/", species_name, "_SD_stochastic_uncertainty.tif"),
  overwrite = TRUE
)

# Load rasters
total_specialist <- rast("./projections/Delphinium occidentale_SD_total.tif")
total_generalist <- rast("./projections/Campanula rotundifolia_SD_total.tif")

method_specialist <- rast("./projections/Delphinium occidentale_SD_method_uncertainty.tif")
method_generalist <- rast("./projections/Campanula rotundifolia_SD_method_uncertainty.tif")

stochastic_specialist <- rast("./projections/Delphinium occidentale_SD_stochastic_uncertainty.tif")
stochastic_generalist <- rast("./projections/Campanula rotundifolia_SD_stochastic_uncertainty.tif")

# Load administrative boundaries
us_states <- gadm(country = "USA", level = 1, path = tempdir())
can_provinces <- gadm(country = "CAN", level = 1, path = tempdir())

admin_lines <- rbind(us_states, can_provinces)

# match CRS to rasters
admin_lines <- project(admin_lines, crs(total_generalist))

# Shared color scale for comparison
all_vals <- c(
  values(total_specialist),
  values(total_generalist),
  values(method_specialist),
  values(method_generalist),
  values(stochastic_specialist),
  values(stochastic_generalist)
)

zlim <- range(all_vals, na.rm = TRUE)

uncertainty_palette <- colorRampPalette(
  c("bisque", "burlywood1", "sienna3", "saddlebrown", "grey10"))(100)

# Plot function
plot_map <- function(r, title) {
  plot(
    r,
    col = uncertainty_palette,
    zlim = zlim,
    main = title,
    axes = TRUE,
    legend = FALSE
  )
  
  plot(admin_lines, add = TRUE, border = "gray40", lwd = 0.5)
}

# Create figure (FIXED LAYOUT)
png("Uncertainty_Maps.png", width = 16, height = 10, units = "in", res = 300)

layout(matrix(c(1,2,3,7,
                4,5,6,7), nrow = 2, byrow = TRUE),
       widths = c(1,1,1,0.3))

par(mar = c(3, 3, 3, 2))

plot_map <- function(r, title, show_x = FALSE, show_y = FALSE) {
  plot(
    r,
    col = uncertainty_palette,
    zlim = zlim,
    main = title,
    axes = TRUE,
    legend = FALSE
  )
  
  plot(admin_lines, add = TRUE, border = "gray40", lwd = 0.5)
  
  if (show_x) mtext("Longitude", side = 1, line = 2, cex = 0.9)
  if (show_y) mtext("Latitude", side = 2, line = 2, cex = 0.9)
}

# Row 1 (Generalist)
plot_map(total_generalist, "A) Generalist – Total", show_y = TRUE)
plot_map(method_generalist, "B) Generalist – Method")
plot_map(stochastic_generalist, "C) Generalist – Stochastic")

# Row 2 (Specialist)
plot_map(total_specialist, "D) Specialist – Total", show_x = TRUE, show_y = TRUE)
plot_map(method_specialist, "E) Specialist – Method", show_x = TRUE)
plot_map(stochastic_specialist, "F) Specialist – Stochastic", show_x = TRUE)

dev.off()

library(fields)

png("Uncertainty_Legend.png", width = 3, height = 6, units = "in", res = 300)

par(mar = c(2, 2, 3, 2))  # balanced margins

fields::image.plot(
  legend.only = TRUE,
  zlim = zlim,
  col = uncertainty_palette,
  
  smallplot = c(0.48, 0.52, 0.1, 0.9),

  legend.lab = "Uncertainty\n(SD of predicted habitat suitability)",
  legend.line = 3,
  
  axis.args = list(
    at = seq(zlim[1], zlim[2], length.out = 5),
    labels = round(seq(zlim[1], zlim[2], length.out = 5), 0),
    las = 1,
    cex.axis = 0.8
  )
)

dev.off()
