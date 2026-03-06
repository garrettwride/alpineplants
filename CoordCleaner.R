library(CoordinateCleaner)
library(rgbif)
library(dplyr)
library(sf)
library(tidyr)

#check diff between using occ_search and occ_download
#maybe add basisOfRecord, remove duplicates

#reads .shp into sf object
rocky_poly <- st_read("./RockyMountainsRegion/rocky_mountains.shp")
rocky_poly <- st_transform(rocky_poly, 4326)  # Make sure CRS matches occurrences
rocky_wkt <- st_as_text(st_union(rocky_poly))

clean_and_get_occurrences <- function(taxon_key, species_name, rocky_poly, rocky_wkt) {

  occ <- occ_search(
    taxonKey = taxon_key, 
    hasCoordinate = TRUE,
    hasGeospatialIssue = FALSE, 
    geometry = rocky_wkt, 
    limit = 100000
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

#graminoids:
showy_sedge <- clean_and_get_occurrences(2723145, "Carex spectabilis", rocky_poly, rocky_wkt)
rm_sedge <- clean_and_get_occurrences(2722910, "Carex arapahoensis", rocky_poly, rocky_wkt)
baldy_sedge <- clean_and_get_occurrences(2723300, "Carex perglobosa", rocky_poly, rocky_wkt)
drummond_rush <- clean_and_get_occurrences(2701717, "Juncus drummondii", rocky_poly, rocky_wkt)
beardless_wheatgrass <- clean_and_get_occurrences(2705699, "Pseudoroegneria spicata", rocky_poly, rocky_wkt)
tufted_hairgrass <- clean_and_get_occurrences(8059811, "Deschampsia cespitosa", rocky_poly, rocky_wkt)
baltic_rush <- clean_and_get_occurrences(2702048, "Juncus balticus", rocky_poly, rocky_wkt)
spiked_rush <- clean_and_get_occurrences(2700795, "Luzula spicata", rocky_poly, rocky_wkt)
single_spike_sedge <- clean_and_get_occurrences(2721925, "Carex scirpoidea", rocky_poly, rocky_wkt)
tufted_hair_grass <- clean_and_get_occurrences(8059811, "Deschampsia cespitosa", rocky_poly, rocky_wkt)

#herbs:
moss_campion_data <- clean_and_get_occurrences(5384754, "Silene acaulis", rocky_poly, rocky_wkt)
moss_campion_data$resp.var
moss_campion_data$resp.xy
moss_campion_data$resp.name

tall_larkspur <- clean_and_get_occurrences(3033713, "Delphinium occidentale", rocky_poly, rocky_wkt)
shootingstar_data <- clean_and_get_occurrences(5641386, "Primula utahensis", rocky_poly, rocky_wkt)
bistort <- clean_and_get_occurrences(2889299, "Bistorta vivipara", rocky_poly, rocky_wkt)
bluebell <- clean_and_get_occurrences(5410907, "Campanula rotundifolia", rocky_poly, rocky_wkt)
mint <- clean_and_get_occurrences(2927192, "Mentha arvensis", rocky_poly, rocky_wkt)
woolly_pussytoes <- clean_and_get_occurrences(5385629, "Antennaria lanata", rocky_poly, rocky_wkt)
butter_and_eggs <- clean_and_get_occurrences(5415020, "Linaria vulgaris", rocky_poly, rocky_wkt)
monkshood <- clean_and_get_occurrences(3033668, "Aconitum columbianum", rocky_poly, rocky_wkt)
alpine_avens <- clean_and_get_occurrences(5369874, "Geum rossii", rocky_poly, rocky_wkt)

#shrubs:
cinque_foil <- clean_and_get_occurrences(5370380, "Dasiphora fruticosa", rocky_poly, rocky_wkt)
hackberry_data <- clean_and_get_occurrences(6406316, "Celtis reticulata", rocky_poly, rocky_wkt)
willow_data <- clean_and_get_occurrences(5372756, "Salix petrophila", rocky_poly, rocky_wkt)
sagebrush <- clean_and_get_occurrences(9396703, "Artemisia tridentata", rocky_poly, rocky_wkt)
serviceberry <- clean_and_get_occurrences(3023964, "Amelanchier utahensis", rocky_poly, rocky_wkt)
bitterbrush <- clean_and_get_occurrences(5370313, "Purshia tridentata", rocky_poly, rocky_wkt)
green_normon_tea <- clean_and_get_occurrences(2653272, "Ephedra viridis", rocky_poly, rocky_wkt)
oregon_boxleaf <- clean_and_get_occurrences(3169107, "Paxistima myrsinites", rocky_poly, rocky_wkt)
rabbitbrush <- clean_and_get_occurrences(3148624, "Chrysothamnus", rocky_poly, rocky_wkt)
gambel_oak <- clean_and_get_occurrences(2878268, "Quercus gambelii", rocky_poly, rocky_wkt)

#trees:
maple_data <- clean_and_get_occurrences(3189864, "Acer glabrum", rocky_poly, rocky_wkt)
fir_data <- clean_and_get_occurrences(2685313, "Abies lasiocarpa", rocky_poly, rocky_wkt)
green_ash <- clean_and_get_occurrences(3012298, "Sorbus scopulina", rocky_poly, rocky_wkt)
ash_leaf_maple <- clean_and_get_occurrences(3189866, "Acer negundo", rocky_poly, rocky_wkt)
pinon <- clean_and_get_occurrences(5285796, "Pinus edulis", rocky_poly, rocky_wkt)
bristlecone <- clean_and_get_occurrences(5285258, "Pinus longaeva", rocky_poly, rocky_wkt)
bull_pine <- clean_and_get_occurrences(5285053, "Pinus ponderosa", rocky_poly, rocky_wkt)
colombian_spruce <- clean_and_get_occurrences(5284917, "Picea engelmannii", rocky_poly, rocky_wkt)
bolander_pine <- clean_and_get_occurrences(5285750, "Pinus contorta", rocky_poly, rocky_wkt)
juniper <- clean_and_get_occurrences(2684709, "Juniperus communis", rocky_poly, rocky_wkt)

#for BIOMOD input:
species_list <- list(
  moss_campion_data,
  shootingstar_data,
  maple_data,
  fir_data,
  hackberry_data,
  willow_data,
  showy_sedge,
  rocky_sedge,
  baldy_sedge,
  tall_larkspur,
  green_ash,
  cinque_foil#,
  
  #all the rest of the species, 10 each:
  # drummond_rush, 
  # beardless_wheatgrass, 
  # tufted_hairgrass, 
  # baltic_rush, 
  # spiked_rush, 
  # single_spike_sedge, 
  # tufted_hair_grass,
  # bistort,
  # bluebell, 
  # mint, 
  # woolly_pussytoes, 
  # butter_and_eggs, 
  # monkshood, 
  # alpine_avens,
  # sagebrush,
  # serviceberry,
  # bitterbrush, 
  # green_normon_tea,
  # oregon_boxleaf,
  # rabbitbrush,
  # gambel_oak,
  # ash_leaf_maple,
  # pinon, 
  # bristlecone, 
  # bull_pine,
  # colombian_spruce, 
  # bolander_pine,
  # juniper
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


