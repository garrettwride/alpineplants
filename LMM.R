library(lme4)
library(r2glmm)
library(emmeans)
library(ggplot2)
library(dplyr)
library(gridExtra)

Create_R_Graph <- function(r2_object, title) {
  df <- data.frame(r2_object)
  df <- df[!(df$Effect == "Model" | df$Effect == "NicheBreadth:Algorithm"),]
  
  plot <- ggplot(data = df, aes(x = Effect, y = Rsq)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL)) +
    ggtitle(title) +
    ylab("Variance Explained") +
    xlab("Predictor") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)
}


biomod_data <- read.csv("./testFiles/SDM_results_for_LMM_2.csv")

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
schoener_resp_var_lmer <- lmer(normalizedSchoenersD ~ 
        NicheBreadth + Algorithm + PAStrategy + PANumber +
        NicheBreadth:Algorithm + NicheBreadth:Occurrences +
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# Response Variable = Mean Suitability
meansuit_resp_var_lmer <- lmer(normalizedMeanSuitability ~ 
        NicheBreadth + Algorithm + PAStrategy + PANumber +
        NicheBreadth:Algorithm + NicheBreadth:Occurrences +
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# Response Variable = Suitable Area
suitarea_resp_var_lmer <- lmer(normalizedSuitableArea ~ 
        NicheBreadth + Algorithm  + PAStrategy + PANumber +
        NicheBreadth:Algorithm +  NicheBreadth:Occurrences +
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# Response Variable =  Elevation Centroid
elevcentroid_resp_var_lmer <- lmer(normalizedElevationCentroid ~ 
        NicheBreadth + Algorithm + PAStrategy + PANumber + 
        NicheBreadth:Algorithm + NicheBreadth:Occurrences +
        (1 | Species) + (1 | Species:Bootstrap),
        data = biomod_data,
        REML = FALSE )

# summary(schoener_resp_var)
# summary(meansuit_resp_var)
# summary(suitarea_resp_var)
# summary(elevcentroid_resp_var)

# MARGINAL calculation of PERCENT VARIANCE (only fixed variables)
schoener_var_r2 = r2beta(schoener_resp_var_lmer, 
                              partial = TRUE, 
                              method = "nsj" )

meansuit_var_r2 = r2beta(meansuit_resp_var_lmer, 
                              partial = TRUE, 
                              method = "nsj" )

suitarea_var_r2 = r2beta(suitarea_resp_var_lmer, 
                         partial = TRUE, 
                         method = "nsj" )

elevcentroid_var_r2 = r2beta(elevcentroid_resp_var_lmer, 
                         partial = TRUE, 
                         method = "nsj" )

# plotting variances 
schoener_plot <- Create_R_Graph(schoener_var_r2, "Marginal R² for Schoener's D")
meansuit_plot <- Create_R_Graph(meansuit_var_r2, "Marginal R² for Mean Suitablity")
suitarea_plot <- Create_R_Graph(suitarea_var_r2, "Marginal R² for Suitable Area")
elevcentroid_plot <- Create_R_Graph(elevcentroid_var_r2, "Marginal R² for Elevation Centroid")

grid.arrange(schoener_plot, meansuit_plot, suitarea_plot, elevcentroid_plot, ncol = 2, nrow = 2)

# CONDITIONAL calculation of PERCENT VARIANCE (fixed AND random variables)
schoener_resp_var_lmer

# How each ResponceVariable changes when NicheBreadth changes, within each Bootstrap
#emt_bootstrap_schoener = emtrends(schoener_resp_var_lmer, ~Species:Bootstrap, var = "NicheBreadth")
#emt_bootstrap_meansuit = emtrends(meansuit_resp_var_lmer, ~Species:Bootstrap, var = "NicheBreadth")
#emt_bootstrap_suitarea = emtrends(suitarea_resp_var_lmer, ~Species:Bootstrap, var = "NicheBreadth")
#emt_bootstrap_elevcentroid = emtrends(elevcentroid_resp_var_lmer, ~Species:Bootstrap, var = "NicheBreadth")

# How each ResponceVariable changes when NicheBreadth changes, within each Algorithm
emt_algorithm_schoener = emtrends(schoener_resp_var_lmer, ~Algorithm, var = "NicheBreadth")
emt_algorithm_meansuit = emtrends(meansuit_resp_var_lmer, ~Algorithm, var = "NicheBreadth")
emt_algorithm_suitarea = emtrends(suitarea_resp_var_lmer, ~Algorithm, var = "NicheBreadth")
emt_algorithm_elevcentroid = emtrends(elevcentroid_resp_var_lmer, ~Algorithm, var = "NicheBreadth")

# How each ResponceVariable changes when NicheBreadth changes, within each PAStrategy
emt_paStrategy_schoener = emtrends(schoener_resp_var_lmer, ~PAStrategy, var = "NicheBreadth")
emt_paStrategy_meansuit = emtrends(meansuit_resp_var_lmer, ~PAStrategy, var = "NicheBreadth")
emt_paStrategy_suitarea = emtrends(suitarea_resp_var_lmer, ~PAStrategy, var = "NicheBreadth")
emt_paStrategy_elevcentroid = emtrends(elevcentroid_resp_var_lmer, ~PAStrategy, var = "NicheBreadth")

# How each ResponceVariable changes when NicheBreadth changes, within each PANumber
emt_paNumber_schoener = emtrends(schoener_resp_var_lmer, ~PANumber, var = "NicheBreadth")
emt_paNumber_meansuit = emtrends(meansuit_resp_var_lmer, ~PANumber, var = "NicheBreadth")
emt_paNumber_suitarea = emtrends(suitarea_resp_var_lmer, ~PANumber, var = "NicheBreadth")
emt_paNumber_elevcentroid = emtrends(elevcentroid_resp_var_lmer, ~PANumber, var = "NicheBreadth")

# uses asymptotic confidence intervals (normal distribution)
schoener_df = as.data.frame(emt_algorithm_schoener) %>%
  mutate(ResponseVariable = "Schoener's D")
schoener_df

meansuit_df = as.data.frame(emt_algorithm_meansuit) %>%
  mutate(ResponseVariable = "Mean Suitablilty")

suitarea_df = as.data.frame(emt_algorithm_suitarea) %>%
  mutate(ResponseVariable = "Suitable Area")

elevcentroid_df = as.data.frame(emt_algorithm_elevcentroid) %>%
  mutate(ResponseVariable = "Elevation Centroid")

combined_bootstrap_df = bind_rows(schoener_df, meansuit_df, suitarea_df, elevcentroid_df)

write.csv(combined_bootstrap_df, file = "LMMResults/combined_bootstrap_df.csv", row.names = FALSE)



ggplot() +
  geom_point(
    data = biomod_data,
    aes(NicheBreadth, SchoenersD, color = Algorithm),
    alpha = 0.15
  ) +
  geom_point(
    data = combined_bootstrap_df,
    aes(NicheBreadth.trend, Algorithm, color = ResponseVariable),
    size = 4
  ) +
  geom_errorbarh(data = combined_bootstrap_df, aes(xmin = asymp.LCL, xmax = asymp.UCL), width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "NicheBreadth Effect (Slope)", y = "Algorithm", title = "Effect of NicheBreadth on Response Variables by Algorithm")


