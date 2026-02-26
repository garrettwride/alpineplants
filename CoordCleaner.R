library(CoordinateCleaner)
library(rgbif)
library(dplyr)
library(sf)
library(tidyr)

#herbs: Silene acaulis (G): 5384754; Primula utahensis (S): 5641386
#trees: Acer glabrum (G): 3189864; Abies lasiocarpa (S):2685313 
#shrubs: Celtis reticulata (G): 6406316; Salix petrophila (S): 5372756

#check diff between using occ_search and occ_download
#maybe add basisOfRecord, remove duplicates

#reads .shp into sf object
rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)  # Make sure CRS matches occurrences
rocky_wkt <- st_as_text(st_union(rocky_poly))

clean_and_get_occurrences <- function(taxon_key, species_name, rocky_poly, rocky_wkt) {

  occ <- occ_search(
    taxonKey = taxon_key, 
    #scientificName = species_name, not needed
    hasCoordinate = TRUE,
    hasGeospatialIssue = FALSE, 
    geometry = rocky_wkt, 
    limit = 100000
    #basisOfRecord = "HUMAN_OBSERVATION", "MACHINE_OBSERVATION",
  )
  
  # Extract the data frame
  df <- occ$data
  
  if(nrow(df) == 0) {
    stop("no record found")
  }
  
  cleaned <- clean_coordinates(
    x = df,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    species = "species",
    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers", "seas",
              "zeros"),
  )
  
  #only keep good points (clean_coords() just evaluates)
  #so each row has a geometry point
  df_cleaned <- df[cleaned$.summary == TRUE, ]
  
  #if remove dupes, do here
  
  #st_as_sf converts foreign object to a spatial object
  occ_sf <- st_as_sf(
    df_cleaned,
    coords = c("decimalLongitude", "decimalLatitude"),
    remove = FALSE, #do not remove coords from df
    crs = 4326 #same for each?
  )
  
  #only keep points inside polygon (only rocky mountains)
  #don't need bc of gbif geometry parameter,
  #but there may be differences bc gbif uses its own spatial engine
  occ_rocky <- occ_sf[st_within(occ_sf, rocky_poly, sparse = FALSE), ]
  
  #remove spatial geo col, back to normal df
  occ_rocky_df <- st_drop_geometry(occ_rocky)
  
  #BIOMOD inputs:
  
  # Presence points
  resp.var <- rep(1, nrow(occ_rocky_df))
  
  # Coordinates
  resp.xy <- occ_rocky_df[, c("decimalLongitude", "decimalLatitude")]
  
  # Name for BIOMOD
  resp.name <- species_name
  
  return(list(
    resp.var = resp.var,
    resp.xy = resp.xy,
    resp.name = resp.name
  ))

} 

#4,304 occurrences
moss_campion_data <- clean_and_get_occurrences(5384754, "Silene acaulis", rocky_poly, rocky_wkt)
moss_campion_data$resp.var
moss_campion_data$resp.xy
moss_campion_data$resp.name

#1,608
tall_larkspur <- clean_and_get_occurrences(3033713, "Delphinium occidentale", rocky_poly, rocky_wkt)
tall_larkspur$resp.xy

#10
shootingstar_data <- clean_and_get_occurrences(5641386, "Primula utahensis", rocky_poly, rocky_wkt)
shootingstar_data$resp.xy

#6,170
maple_data <- clean_and_get_occurrences(3189864, "Acer glabrum", rocky_poly, rocky_wkt)
maple_data$resp.xy

#5,211
fir_data <- clean_and_get_occurrences(2685313, "Abies lasiocarpa", rocky_poly, rocky_wkt)
fir_data$resp.xy

#1,470
green_ash <- clean_and_get_occurrences(3012298, "Sorbus scopulina", rocky_poly, rocky_wkt)
green_ash$resp.xy

#9,059
cinque_foil <- clean_and_get_occurrences(5370380, "Dasiphora fruticosa", rocky_poly, rocky_wkt)
cinque_foil$resp.xy

#323
hackberry_data <- clean_and_get_occurrences(6406316, "Celtis reticulata", rocky_poly, rocky_wkt)
hackberry_data$resp.xy

#518
willow_data <- clean_and_get_occurrences(5372756, "Salix petrophila", rocky_poly, rocky_wkt)
willow_data$resp.xy

#238
showy_sedge <- clean_and_get_occurrences(2723145, "Carex spectabilis", rocky_poly, rocky_wkt)
showy_sedge$resp.xy

#110
rocky_sedge <- clean_and_get_occurrences(2722910, "Carex arapahoensis", rocky_poly, rocky_wkt)
rocky_sedge$resp.xy

#147
baldy_sedge <- clean_and_get_occurrences(2723300, "Carex perglobosa", rocky_poly, rocky_wkt)
baldy_sedge$resp.xy

#for BIOMOD input:
species_list <- list(
  moss_campion_data,
  shootingstar_data,
  maple_data,
  fir_data,
  hackberry_data,
  willow_data#,
  #start with 6 species, add more when using supercomputer(?)
  #showy_sedge,
  #rocky_sedge,
  #baldy_sedge
)

long_data <- species_list %>%
  map_df(function(species) {
    data.frame(
      X_WGS84 = species$resp.xy$decimalLongitude,
      Y_WGS84 = species$resp.xy$decimalLatitude,
      Species = make.names(species$resp.name),
      Presence = 1
    )
  })

DataSpecies <- long_data %>%
  distinct() %>%
  pivot_wider(
    names_from = Species,
    values_from = Presence,
    values_fill = 0
  )

head(DataSpecies)


