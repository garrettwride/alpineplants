pkgs <- c(
  "biomod2", "CoordinateCleaner", "dismo", "dplyr", "elevatr", 
  "emmeans", "geodata", "ggplot2", "gridExtra", "lme4", 
  "purrr", "r2glmm", "reshape2", "rgbif", "rJava", 
  "R.utils", "sf", "terra", "tidyr", "partR2"
)

if (!requireNamespace("biomod2", quietly = TRUE)) install.packages("biomod2")
if (!requireNamespace("CoordinateCleaner", quietly = TRUE)) install.packages("CoordinateCleaner")
if (!requireNamespace("dismo", quietly = TRUE)) install.packages("dismo")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("elevatr", quietly = TRUE)) install.packages("elevatr")
if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
if (!requireNamespace("geodata", quietly = TRUE)) install.packages("geodata")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")
if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
if (!requireNamespace("r2glmm", quietly = TRUE)) install.packages("r2glmm")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("rgbif", quietly = TRUE)) install.packages("rgbif")
if (!requireNamespace("rJava", quietly = TRUE)) install.packages("rJava")
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("partR2", quietly = TRUE)) install.packages("partR2")
