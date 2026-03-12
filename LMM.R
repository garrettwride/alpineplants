library(lme4)
library(r2glmm)
library(emmeans)
library(ggplot2)
library(dplyr)
library(gridExtra)

Create_R_Graph <- function(r2_object, title) {
  df <- data.frame(r2_object)
  df <- df[!(df$Effect == "Model" | df$Effect == "NicheBreadth:Algorithm" | df$Effect == "NicheBreadth:PAStrategy" | df$Effect == "NicheBreadth:LogOccurrences" | df$Effect == "NicheBreadth:PANumber"),]
  
  plot <- ggplot(data = df, aes(x = Effect, y = Rsq)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL)) +
    ggtitle(title) +
    ylab("Variance Explained") +
    xlab("Predictor") +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 20, hjust = 0.8))
  
  return(plot)
}


Run_LMM <- function(biomod_data, response_variable, output_file_path) {
  # for SPECIALISTS VS GENERALISTS as a whole (different Niche Scores): 
  # Random Effects = Species, Bootstrap
  # Fixed Effects = NicheBreadth, Algorithm, PAStrategy, PANumber
  lmer_output <- lmer(as.formula(paste0(response_variable, " ~ 
                    NicheBreadth + Algorithm + PAStrategy + PANumber+ LogOccurrences +
                    NicheBreadth:Algorithm + 
                    NicheBreadth:PAStrategy + 
                    NicheBreadth:PANumber + 
                    NicheBreadth:LogOccurrences +
                    (1 | Species) + (1 | Species:Bootstrap)")),
                    data = biomod_data,
                    REML = FALSE )
  
  return(lmer_output)
}


Compute_and_Graph_EMTrends <- function(fixed_var, schoener_lmer, meansuit_lmer, suitarea_lmer,elevcentroid_lmer, output_file_path) {
  emt_schoener = emtrends(schoener_lmer, as.formula(paste0("~", fixed_var)), var = "NicheBreadth")
  emt_meansuit = emtrends(meansuit_lmer, as.formula(paste0("~", fixed_var)), var = "NicheBreadth")
  emt_suitarea = emtrends(suitarea_lmer, as.formula(paste0("~", fixed_var)), var = "NicheBreadth")
  emt_elevcentroid = emtrends(elevcentroid_lmer, as.formula(paste0("~", fixed_var)), var = "NicheBreadth")
  
  schoener_df = as.data.frame(emt_schoener) %>%
    mutate(ResponseVariable = "Schoener's D")
  
  meansuit_df = as.data.frame(emt_meansuit) %>%
    mutate(ResponseVariable = "Mean Suitablilty")
  
  suitarea_df = as.data.frame(emt_suitarea) %>%
    mutate(ResponseVariable = "Suitable Area")
  
  elevcentroid_df = as.data.frame(emt_elevcentroid) %>%
    mutate(ResponseVariable = "Elevation Centroid")
  
  combined_df = bind_rows(schoener_df, meansuit_df, suitarea_df, elevcentroid_df)
  
  write.csv(combined_df, file = output_file_path, row.names = FALSE)
  
  return(combined_df)
}


biomod_data <- read.csv("./testFiles/SDM_results_for_LMM_2.csv")

# makes PAStrategy and PANumber categorical variables 
biomod_data$PAStrategy <- as.factor(biomod_data$PAStrategy)
biomod_data$PANumber <- as.factor(biomod_data$PANumber)


#run LMM 
schoener_resp_var_lmer <- Run_LMM(biomod_data, "SchoenersD") 
meansuit_resp_var_lmer <- Run_LMM(biomod_data, "MeanSuitability")
suitarea_resp_var_lmer <- Run_LMM(biomod_data, "SuitableArea")
elevcentroid_resp_var_lmer <- Run_LMM(biomod_data, "ElevationCentroid")

summary(schoener_resp_var_lmer)
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

# CONDITIONAL calculation of PERCENT VARIANCE (fixed AND random variables) <- this is no longer correct(change)


algorithm_emtrends_results <- Compute_and_Graph_EMTrends("Algorithm", schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, "LMMResults/algorithm_emtrends_results.csv")
paStrategy_emtrends_results <- Compute_and_Graph_EMTrends("PAStrategy", schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, "LMMResults/paStrategy_emtrends_results.csv")
paNumber_emtrends_results <- Compute_and_Graph_EMTrends("PANumber", schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, "LMMResults/paNumber_emtrends_results.csv")


algorithm_emtrends_results



# ggplot(algorithm_emtrends_results, aes(x = NicheBreadth.trend, y = Algorithm, color = ResponseVariable)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL), width = 0.2) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   labs(
#     x = "NicheBreadth Effect (Slope)",
#     y = "Algorithm",
#     title = "Effect of NicheBreadth on Response Variables by Algorithm")
# 

