Install Dependencies
- Run `./R/preprocessing/DownloadDependencies.R`

Obtain Rocky Mountain Region Raster Shape
- Go to https://www.naturalearthdata.com/downloads/50m-physical-vectors/50m-physical-labels/
  - Download “Download label areas (1.55 MB) version 5.0.0”
  - Open the zip file on your computer
    - Make sure it is called `ne_50m_geography_regions_polys`
  - Save the unzipped file in `./data` folder
- Run `./R/preprocessing/RockyMountainShp.R` script
  - This will plot the shape so you can verify it looks correct, save the shape if not already saved, and make sure that the saved file can be opened
  - The shape will be saved in `./data/RockyMountainsRegion/rocky_mountains.shp`

Get and Clean Species Occurrence Data and Environmental Data
- Run `./R/preprocessing/downloadDataForTAs.R`
  - This is a modified version of `./R/preprocessing/downloadData.R` for TAs to download 2 species rather than 28 in order to reduce runtime for the reproducibility challenge
  - This also downloads WorldClim environmental variables as well as elevation data in the format of rasters
  - This script uses the Rocky Mountain polygon from the previous step
  - These files should be in the `./data` directory and are used as the input data for the BIOMOD pipeline:
    - `species_occurrences.rds`
    - `myExpl.rds`
    - `elev_rocky.rds`

Run BioMOD Pipeline
- Run `./R/modeling/biomodPipelineForTAs.R`
  - This is a modified version of `./R/modeling/biomodPipeline.R`. The script for TAs trains fewer models to reduce runtime for the reproducibility challenge
  - This script will produce a file called `./data/SDM_results_for_TAs.csv`. This file demonstrates that the pipeline works, but should not be used in later steps
  - Use `./data/SDM_results_for_LMM_3.csv`, which are the actual results used in our analysis

Run Analysis and Save results
- Run `./R/modeling/LMMandEMTrends.R`
- Results will be located in:
  - `./output/EMTrendsResults`
  - `./output/LMMResults`
  - `./output/R2BetaResults`

Create Analysis Tables
- Run `./R/analysis/PlottingResults.R`
- Charts can be looked at within the RStudio window

Generate Figure 2
- Run `./R/analysis/Generate_Fig_2.R`
- View the figures in `./output/Fig2`

Generate Figure 3
- Run `./R/analysis/Generate_Fig_3.R`
- Figure will be saved in `./output/Fig3`