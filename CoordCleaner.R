install.packages("CoordinateCleaner")
install.packages("rgbif")
install.packages("sf")

library(CoordinateCleaner)
library(rgbif)

occ_search(taxonKey=5384754, country="US;CA", basisOfRecord="HUMAN_OBSERVATION")
