install.packages("CoordinateCleaner")
install.packages("rgbif")
install.packages("dplyr")
install.packages("usethis")
usethis::edit_r_environ()

library(CoordinateCleaner)
library(rgbif)
library(dplyr)

#herbs:
#occ_search is for small, datasets; for real project use occ_download
#hasCoordinate=TRUE
#use coordCleaner functions: cc_val() valid coords, cc_equ() removes 0s, cc_inst() remove invalid or university/garden coords
#maybe use cc_cen to remove centroids of countries
#moss campion: occ_search(taxonKey=5384754, country="US;CA", basisOfRecord="HUMAN_OBSERVATION")
#occ_search(taxonKey=5641386, country="US;CA", basisOfRecord ="HUMAN_OBSERVATION")

#req gbif data
res <- occ_download(
  pred("taxonKey", 5384754),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=GBIF_USER
)

occ_download_get(res, path = ".")
data = occ_download_import(res)

#clean data here
#use coordCleaner functions: cc_val() valid coords, cc_equ() removes 0s, cc_inst() remove invalid or university/garden coords
#cc_sea removes coords from ocean
#maybe use cc_cen to remove centroids of countries
cleaned_data = cc_val(data, lon = "decimalLongitude", lat = "decimalLatitude")
cleaned_data = cc_zero(cleaned_data, lon = "decimalLongitude", lat = "decimalLatitude")
cleaned_data = cc_inst(cleaned_data, lon = "decimalLongitude", lat = "decimalLatitude", buffer = 100) #buffer is radius in m around institution
cleaned_data = cc_sea(cleaned_data, lon = "decimalLongitude", lat = "decimalLatitude")


