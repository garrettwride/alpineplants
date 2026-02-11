library(biomod2)
library(terra)

algorithms   <- c("GLM", "GAM", "RF") # "MAXENT"
pa_strategies <- c("random", "disk")
pa_numbers    <- c(1, 5, 10)   # × presences
n_boot        <- 10

# outermost loop
# for (species in species_list) 
#     myRespName <- species$name
      # myResp <- species$response
      # myRespXY <- species$coordinates
      # SpeciesType <- species$type # Generalist/Specialist

# Bootsrap occurance data (next loop)
# for (b in range(n_boot))
  # set.seed(b)
  # pres_idx <- which(myResp == 1)
  # boot_pres <- sample(pres_idx, replace = TRUE)
  # 
  # boot_resp <- rep(0, length(myResp))
  # boot_resp[boot_pres] <- 1

# Set up Pseudo-absences (next 2 loops)
for (pa_strat in pa_strategies) {
  for (pa_mult in pa_numbers) {
    PA_n <- pa_mult * length(pres_idx)
    
    myBiomodData <- BIOMOD_FormatingData(
      resp.var = boot_resp,
      expl.var = myExpl,
      resp.xy  = myRespXY,
      resp.name = paste(
        myRespName, "boot", b, pa_strat, pa_mult, sep = "_"),
      PA.nb.rep = 1,
      PA.nb.absences = PA_n,
      PA.strategy = pa_strat,
      PA.dist.min = 10,
      PA.dist.max = 20
    )
    
    myBiomodModelOut <- BIOMOD_Modeling(
      bm.format = myBiomodData,
      models = algorithms,
      CV.strategy = "random",
      CV.nb.rep = 1,
      CV.perc = 0.8,
      CV.do.full.models = TRUE,
      metric.eval = c("TSS", "AUC")
    )
    ## Save the rasters/analyze them
    myBiomodProj <- BIOMOD_Projection(
      bm.mod = myBiomodModelOut,
      proj.name = paste0("Current_boot", b),
      new.env = myExpl,
      models.chosen = "all",
      metric.binary = "TSS",
      metric.filter = "all",
      build.clamping.mask = TRUE
    )
  }  
}

#     BIOMOD_FormatingData (PA.nb.rep = 1)
#     BIOMOD_Modeling (CV optional)
#     BIOMOD_Projection

# Load species occurrences (6 species available)
data('DataSpecies')
head(DataSpecies)

# Select the name of the studied species
myRespName <- 'GuloGulo'

# Get corresponding presence/absence data
myResp <- as.numeric(DataSpecies[, myRespName])

# Get corresponding XY coordinates
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]

# Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
data('bioclim_current')
myExpl <- rast(bioclim_current)
print(myExpl)

# # Transform true absences into potential pseudo-absences
myResp.PA <- ifelse(myResp == 1, 1, NA)

# # Format data with pseudo-absences : random method
myBiomodData.r <- BIOMOD_FormatingData(dir.name = 'testFiles/',
                                       resp.var = myResp.PA,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName,
                                       PA.nb.rep = 1,
                                       PA.nb.absences = 100,
                                       PA.strategy = 'disk',
                                       PA.dist.min = 10,
                                       PA.dist.max = 20)

myBiomodData.r
plot(myBiomodData.r)


# Model single models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData.r,
                                    models = c('GLM'),
                                    CV.strategy = 'random',
                                    CV.perc = 0.8,
                                    CV.do.full.models = FALSE,
                                    OPT.strategy = 'default',
                                    metric.eval = 'TSS')
# seed.val = 123)
# nb.cpu = 8)
myBiomodModelOut

# Get evaluation scores & variables importance
get_evaluations(myBiomodModelOut)
get_variables_importance(myBiomodModelOut)

# Project single models
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = myExpl,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
myBiomodProj
plot(myBiomodProj)
get_predictions(myBiomodProj)

# 
# # Model ensemble models
# myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
#                                       models.chosen = 'all',
#                                       em.by = 'all',
#                                       em.algo = c('EMmedian', 'EMmean', 'EMwmean',
#                                                   'EMca', 'EMci', 'EMcv'),
#                                       metric.select = c('TSS'),
#                                       metric.select.thresh = c(0.7),
#                                       metric.eval = c('TSS', 'AUCroc', 'BOYCE'),
#                                       var.import = 3,
#                                       EMci.alpha = 0.05,
#                                       EMwmean.decay = 'proportional')
# myBiomodEM
# 
# # Get evaluation scores & variables importance
# get_evaluations(myBiomodEM)
# get_variables_importance(myBiomodEM)
# 
# # Represent evaluation scores & variables importance
# bm_PlotEvalMean(bm.out = myBiomodEM, dataset = 'calibration', group.by = 'full.name')
# bm_PlotEvalBoxplot(bm.out = myBiomodEM, dataset = 'calibration', group.by = c('full.name', 'full.name'))
# bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'full.name', 'full.name'))
# bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.run'))
# bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'merged.by.run'))
# 
# # Represent response curves
# bm_PlotResponseCurves(bm.out = myBiomodEM, 
#                       models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
#                       fixed.var = 'median')
# bm_PlotResponseCurves(bm.out = myBiomodEM, 
#                       models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
#                       fixed.var = 'min')
# bm_PlotResponseCurves(bm.out = myBiomodEM, 
#                       models.chosen = get_built_models(myBiomodEM)[7],
#                       fixed.var = 'median',
#                       do.bivariate = TRUE)