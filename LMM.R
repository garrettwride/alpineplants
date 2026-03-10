library(lme4)
library(r2glmm)
library(emmeans)
library(ggplot2)
library(dplyr)
library(gridExtra)

biomod_data <- read.csv("./testFiles/SDM_results_for_LMM.csv")


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
# Fixed Effects = NicheBreadth, Algorithm, PAStrategy, PANumber

# Response Variable = Schoener's D
schoener_resp_var <- lmer(normalizedSchoenersD ~ 
        NicheBreadth + Algorithm + PAStrategy + PANumber +
        NicheBreadth:Algorithm + 
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )
#schoener_resp_var

# Response Variable = Mean Suitability
meansuit_resp_var <- lmer(normalizedMeanSuitability ~ 
        NicheBreadth + Algorithm + PAStrategy + PANumber +
        NicheBreadth:Algorithm + 
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# Response Variable = Suitable Area
suitarea_resp_var <- lmer(normalizedSuitableArea ~ 
        NicheBreadth + Algorithm  + PAStrategy + PANumber +
        NicheBreadth:Algorithm + 
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# Response Variable =  Elevation Centroid
elevcentroid_resp_var <- lmer(normalizedElevationCentroid ~ 
        NicheBreadth + Algorithm + PAStrategy + PANumber +
        NicheBreadth:Algorithm + 
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

summary(schoener_resp_var)
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
schoener_plot
meansuit_plot <- plot(meansuit_var_r2) + ggtitle("Marginal Mean Predicted Suitability R²")
suitarea_plot <- plot(suitarea_var_r2) + ggtitle("Marginal Suitable Area R²")
elevcentroid_plot <- plot(elevcentroid_var_r2) + ggtitle("Marginal Elevartion Centroid R²")

grid.arrange(schoener_plot, meansuit_plot, suitarea_plot, elevcentroid_plot, ncol = 2, nrow = 2)

# CONDITIONAL calculation of PERCENT VARIANCE (fixed AND random variables)

# How each ResponceVariable changes when NicheBreadth changes, within each Bootstrap
emt_algorithm_schoener = emtrends(schoener_resp_var, ~Bootstrap, var = "NicheBreadth")
emt_algorithm_meansuit = emtrends(meansuit_resp_var, ~Bootstrap, var = "NicheBreadth")
emt_algorithm_suitarea = emtrends(suitarea_resp_var, ~Bootstrap, var = "NicheBreadth")
emt_algorithm_elevcentroid = emtrends(elevcentroid_resp_var, ~Bootstrap, var = "NicheBreadth")

# How each ResponceVariable changes when NicheBreadth changes, within each Algorithm
emt_algorithm_schoener = emtrends(schoener_resp_var, ~Algorithm, var = "NicheBreadth")
emt_algorithm_meansuit = emtrends(meansuit_resp_var, ~Algorithm, var = "NicheBreadth")
emt_algorithm_suitarea = emtrends(suitarea_resp_var, ~Algorithm, var = "NicheBreadth")
emt_algorithm_elevcentroid = emtrends(elevcentroid_resp_var, ~Algorithm, var = "NicheBreadth")

# How each ResponceVariable changes when NicheBreadth changes, within each PAStrategy
emt_paStrategy_schoener = emtrends(schoener_resp_var, ~PAStrategy, var = "NicheBreadth")
emt_paStrategy_meansuit = emtrends(meansuit_resp_var, ~PAStrategy, var = "NicheBreadth")
emt_paStrategy_suitarea = emtrends(suitarea_resp_var, ~PAStrategy, var = "NicheBreadth")
emt_paStrategy_elevcentroid = emtrends(elevcentroid_resp_var, ~PAStrategy, var = "NicheBreadth")

# How each ResponceVariable changes when NicheBreadth changes, within each PANumber
emt_paNumber_schoener = emtrends(schoener_resp_var, ~PANumber, var = "NicheBreadth")
emt_paNumber_meansuit = emtrends(meansuit_resp_var, ~PANumber, var = "NicheBreadth")
emt_paNumber_suitarea = emtrends(suitarea_resp_var, ~PANumber, var = "NicheBreadth")
emt_paNumber_elevcentroid = emtrends(elevcentroid_resp_var, ~PANumber, var = "NicheBreadth")


schoener_df = as.data.frame(emt_algorithm_schoener) %>%
  mutate(ResponseVariable = "Schoener's D")

meansuit_df = as.data.frame(emt_algorithm_meansuit) %>%
  mutate(ResponseVariable = "Mean Suitablilty")

suitarea_df = as.data.frame(emt_algorithm_suitarea) %>%
  mutate(ResponseVariable = "Suitable Area")

elevcentroid_df = as.data.frame(emt_algorithm_elevcentroid) %>%
  mutate(ResponseVariable = "Elevation Centroid")

combined_df = bind_rows(schoener_df, meansuit_df, suitarea_df, elevcentroid_df)
combined_df

ggplot(combined_df, aes(x = NicheBreadth.trend, y = Algorithm, color = ResponseVariable)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "NicheBreadth Effect (Slope)",
    y = "Algorithm",
    title = "Effect of NicheBreadth on Response Variables by Algorithm")


