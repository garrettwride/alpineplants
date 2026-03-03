library(lme4)
library(r2glmm)
library(emmeans)
library(ggplot2)
library(dplyr)
library(gridExtra)

biomod_data <- read.csv("./testFiles/testLme4Data.csv")

# Normalize MeanSuitability, SuitableArea, SchoenersD, ElevationCentroid
# (Generalists and Specialists combined)
mean_mean_suitability <- mean(biomod_data$MeanSuitability)
sd_mean_suitability <- sd(biomod_data$MeanSuitability)
biomod_data$normalizedMeanSuitability <- scale(biomod_data$MeanSuitability)

mean_suitable_area <- mean(biomod_data$SuitableArea)
sd_suitable_area <- sd(biomod_data$SuitableArea)
biomod_data$normalizedSuitableArea <- scale(biomod_data$SuitableArea)

mean_schoenersD <- mean(biomod_data$SchoenersD)
sd_schoenersD <- sd(biomod_data$SchoenersD)
biomod_data$normalizedSchoenersD <- scale(biomod_data$SchoenersD)

mean_elevation_centroid <- mean(biomod_data$ElevationCentroid)
sd_elevation_centroid <- sd(biomod_data$ElevationCentroid)
biomod_data$normalizedElevationCentroid <- scale(biomod_data$ElevationCentroid)

# for SPECIALISTS VS GENERALISTS as a whole: 
# Random Effects = Species, Bootstrap
# Fixed Effects = SpeciesType, Algorithm, PAStrategy, PANumber

# Response Variable = Schoener's D
schoener_resp_var <- lmer(normalizedSchoenersD ~ 
        SpeciesType + Algorithm + PAStrategy + PANumber +
        SpeciesType:Algorithm + 
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# Response Variable = Mean Suitability
meansuit_resp_var <- lmer(normalizedMeanSuitability ~ 
        SpeciesType + Algorithm + PAStrategy + PANumber +
        SpeciesType:Algorithm + 
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# Response Variable = Suitable Area
suitarea_resp_var <- lmer(normalizedSuitableArea ~ 
        SpeciesType + Algorithm  + PAStrategy + PANumber +
        SpeciesType:Algorithm + 
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# Response Variable =  Elevation Centroid
elevcentroid_resp_var <- lmer(normalizedElevationCentroid ~ 
        SpeciesType + Algorithm + PAStrategy + PANumber +
        SpeciesType:Algorithm + 
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# summary(schoener_resp_var)
# summary(meansuit_resp_var)
# summary(suitarea_resp_var)
# summary(elevcentroid_resp_var)

# CONSIDER PULLING APART THE GROUPS FOR GENERALIST AND SPECIALISTS
#   would have to filter BEFORE putting into lmer() ?

# MARGINAL calculation of PERCENT VARIANCE (only fixed variables)
schoener_var_r2 = r2beta(schoener_resp_var, 
                              partial = TRUE, 
                              method = "nsj" )

meansuit_var_r2 = r2beta(meansuit_resp_var, 
                              partial = TRUE, 
                              method = "nsj" )

suitarea_var_r2 = r2beta(suitarea_resp_var, 
                         partial = TRUE, 
                         method = "nsj" )

elevcentroid_var_r2 = r2beta(elevcentroid_resp_var, 
                         partial = TRUE, 
                         method = "nsj" )

# plotting variances 
schoener_plot <- plot(schoener_var_r2) + ggtitle("Marginal Schoener R²")
meansuit_plot <- plot(meansuit_var_r2) + ggtitle("Marginal Mean Predicted Suitability R²")
suitarea_plot <- plot(suitarea_var_r2) + ggtitle("Marginal Suitable Area R²")
elevcentroid_plot <- plot(elevcentroid_var_r2) + ggtitle("Marginal Elevartion Centroid R²")

grid.arrange(schoener_plot, meansuit_plot, suitarea_plot, elevcentroid_plot, ncol = 2, nrow = 2)

# CONDITIONAL calculation of PERCENT VARIANCE (fixed AND random variables)

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

