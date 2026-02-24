library(lme4)
library(emmeans)
library(ggplot2)
library(dplyr)

biomod_data <- read.csv("./testFiles/testLme4Data.csv")

# for SPECIALISTS VS GENERALISTS as a whole group: 

# Response Variable = Schoener's D
# Random Effects = Species, Bootstrap
# Fixed Effects = SpeciesType, Algorithm, PAStrategy, PANumber
schoener_resp_var <- lmer(SchoenersD ~ 
        SpeciesType + Algorithm + PAStrategy + PANumber 
        + (1 | Species) + (1 | Bootstrap),
        data = biomod_data,
        REML = FALSE
)

# Response Variable = Mean Suitability
meansuit_resp_var <- lmer(MeanSuitability ~ 
        SpeciesType + Algorithm + PAStrategy + PANumber 
        + (1 | Species) + (1 | Bootstrap),
        data = biomod_data,
        REML = FALSE
)

# Response Variable = Suitable Area
suitarea_resp_var <- lmer(SuitableArea ~ 
        SpeciesType + Algorithm  + PAStrategy + PANumber 
        + (1 | Species) + (1 | Bootstrap),
        data = biomod_data,
        REML = FALSE
)

# Response Variable =  Elevation Centroid
elevcentroid_resp_var <- lmer(ElevationCentroid ~ 
        SpeciesType + Algorithm + PAStrategy + PANumber 
        + (1 | Species) + (1 | Bootstrap),
        data = biomod_data,
        REML = FALSE
)

# for SPECIES VS SPECIES
# # Random Effects = Bootstrap
# Fixed Effects = SpeciesType, Algorithm, PAStrategy, PANumber, *Species* (because we care about specific species)

summary(schoener_resp_var)
# summary(meansuit_resp_var)
# summary(suitarea_resp_var)
# summary(elevcentroid_resp_var)

# Predicted Response variable based on 2nd argument (so SpeciesType) 
emm_schoener = emmeans(schoener_resp_var, ~SpeciesType)
emm_meansuit = emmeans(meansuit_resp_var, ~SpeciesType)
emm_suitarea = emmeans(suitarea_resp_var, ~SpeciesType)
emm_elevcentroid = emmeans(elevcentroid_resp_var, ~SpeciesType)


schoener_df = as.data.frame(emm_schoener) %>%
  mutate(ResponseVariable = "Schoener's D")

meansuit_df = as.data.frame(emm_meansuit) %>%
  mutate(ResponseVariable = "Mean Suitablilty")

suitarea_df = as.data.frame(emm_suitarea) %>%
  mutate(ResponseVariable = "Suitable Area")

elevcentroid_df = as.data.frame(emm_elevcentroid) %>%
  mutate(ResponseVariable = "Elevation Centroid")

combined_df = bind_rows(schoener_df, meansuit_df, suitarea_df, elevcentroid_df)
# combined_df

ggplot(combined_df, aes(x = SpeciesType, y = emmean)) +
  geom_bar(stat = "identity",color = "black") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.5) +
  facet_wrap(~ResponseVariable, scales = "free") +
  labs(y = "Response Variable Predictions", x = "Species Type")

