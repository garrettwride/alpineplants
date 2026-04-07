library(ggplot2)
library(readr)
library(gridExtra)
library(dplyr)
library(reshape2)


Plot_Residual_Normality_of_Response_Variables <- function(schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, tss_resp_var_lmer) {
  # Plot Normality of the Residuals (we want straight line)
  qqnorm(resid(schoener_resp_var_lmer), main = "Schoener's D Normal Q-Q Plot")
  qqline(resid(schoener_resp_var_lmer), col = "red")
  
  qqnorm(resid(meansuit_resp_var_lmer), main = "Mean Suitablilty Normal Q-Q Plot")
  qqline(resid(meansuit_resp_var_lmer), col = "red")
  
  qqnorm(resid(suitarea_resp_var_lmer), main = "Suitable Area Normal Q-Q Plot")
  qqline(resid(suitarea_resp_var_lmer), col = "red")
  
  qqnorm(resid(elevcentroid_resp_var_lmer), main = "Elevation Centroid Normal Q-Q Plot")
  qqline(resid(elevcentroid_resp_var_lmer), col = "red")
  
  qqnorm(resid(tss_resp_var_lmer), main = "TSS Normal Q-Q Plot")
  qqline(resid(tss_resp_var_lmer), col = "red")
  
}

Create_Correlation_Heat_Map <- function(lmer_results, fixed_var) {
  summary(lmer_results)
  correlation_matrix <- cov2cor(vcov(lmer_results))
  correlation_matrix <- as.matrix(correlation_matrix)
  correlation_long <- melt(correlation_matrix)
  names(correlation_long) <- c("Coef1", "Coef2", "Correlation")
  correlation_long <- correlation_long[!grepl("^NicheBreadth:", correlation_long$Coef1), ]
  correlation_long <- correlation_long[!grepl("^NicheBreadth:", correlation_long$Coef2), ]
  correlation_long <- correlation_long[!(correlation_long$Coef1 == "(Intercept)"),]
  correlation_long <- correlation_long[!(correlation_long$Coef2 == "(Intercept)"),]
  
  lmerHeatMap <- ggplot(correlation_long, aes(x = Coef1, y = Coef2, fill = Correlation)) +
    geom_tile(color = "lightgrey", linewidth = 0.1) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0,
                         limits = c(-1, 1), 
                         name = "Correlation") +
    labs(title = paste("Fixed Variable Correlations with", fixed_var)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          legend.title = element_text(size = 8),
          axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          axis.title = element_blank()) +
    coord_fixed()
  
  #print(lmerHeatMap)
  return(lmerHeatMap)
}

Create_Variance_Graph_WITHOUT_Interactions <- function(r2_object, title) {
  df <- data.frame(r2_object)
  df <- df[!(df$Effect == "Model"),] 
  df <- df[!grepl("^NicheBreadth:", df$Effect), ] 
  
  plot <- ggplot(data = df, aes(x = Effect, y = Rsq)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL)) +
    coord_cartesian(ylim = c(0, 0.34)) +
    labs(x = "Predictor", y = "Variance Explained", title = title) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 20, hjust = 0.8))
  
  return(plot)
}

Create_Variance_Graph_JUST_Interactions <- function(r2_object, title) {
  df <- data.frame(r2_object)
  df <- df[!(df$Effect == "Model" | df$Effect == "Algorithm" | df$Effect == "LogOccurrences" | df$Effect == "NicheBreadth" | df$Effect == "PANumber" | df$Effect == "PAStrategy"),]
  
  plot <- ggplot(data = df, aes(x = Effect, y = Rsq)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL)) +
    coord_cartesian(ylim = c(0, 0.34)) +
    labs(x = "Predictor", y = "Variance Explained", title = title) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 35, hjust = 1))
  
  return(plot)
}

Plot_R2Beta_Results <- function(schoener_var_r2, meansuit_var_r2, suitarea_var_r2, elevcentroid_var_r2, tss_var_r2) {
  # plotting R2Beta variances WITHOUT interactions
  schoener_plot_without_interactions <- Create_Variance_Graph_WITHOUT_Interactions(schoener_var_r2, "Marginal R² for Schoener's D")
  meansuit_plot_without_interactions <- Create_Variance_Graph_WITHOUT_Interactions(meansuit_var_r2, "Marginal R² for Mean Suitablity")
  suitarea_plot_without_interactions <- Create_Variance_Graph_WITHOUT_Interactions(suitarea_var_r2, "Marginal R² for Suitable Area")
  elevcentroid_plot_without_interactions <- Create_Variance_Graph_WITHOUT_Interactions(elevcentroid_var_r2, "Marginal R² for Elevation Centroid")
  tss_plot_without_interactions <- Create_Variance_Graph_WITHOUT_Interactions(tss_var_r2, "Marginal R² for Model Performance (TSS)")
  
  grid.arrange(schoener_plot_without_interactions, 
               meansuit_plot_without_interactions, 
               suitarea_plot_without_interactions, 
               elevcentroid_plot_without_interactions, 
               tss_plot_without_interactions,
               ncol = 3, nrow = 2)
  
  # plotting R2Beta variances WITH interactions
  schoener_plot_with_interactions <- Create_Variance_Graph_JUST_Interactions(schoener_var_r2, "Marginal R² for Schoener's D")
  meansuit_plot_with_interactions <- Create_Variance_Graph_JUST_Interactions(meansuit_var_r2, "Marginal R² for Mean Suitablity")
  suitarea_plot_with_interactions <- Create_Variance_Graph_JUST_Interactions(suitarea_var_r2, "Marginal R² for Suitable Area")
  elevcentroid_plot_with_interactions <- Create_Variance_Graph_JUST_Interactions(elevcentroid_var_r2, "Marginal R² for Elevation Centroid")
  tss_plot_with_interactions <- Create_Variance_Graph_JUST_Interactions(tss_var_r2, "Marginal R² for Model Performance (TSS)")
  
  grid.arrange(schoener_plot_with_interactions, 
               meansuit_plot_with_interactions, 
               suitarea_plot_with_interactions, 
               elevcentroid_plot_with_interactions, 
               tss_plot_with_interactions,
               ncol = 3, nrow = 2)
}

Plot_EMTrends_Results <- function(emtrends_results, fixed_var) {
  plotEMT <- ggplot(emtrends_results, aes(x = NicheBreadth.trend, y = .data[[fixed_var]])) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL), width = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~ ResponseVariable, scales = "free_x") +
    labs(x = "NicheBreadth Effect (Slope)", y = fixed_var,
         title = paste0("Effect of NicheBreadth on Response Variables by ", fixed_var))
  
  plotEMT
} 


# open LMM results
schoener_resp_var_lmer <- readRDS("output/LMMResults/SchoenersD_lmm_results.rds")
meansuit_resp_var_lmer <- readRDS("output/LMMResults/MeanSuitability_lmm_results.rds")
suitarea_resp_var_lmer <- readRDS("output/LMMResults/SuitableArea_lmm_results.rds")
elevcentroid_resp_var_lmer <- readRDS("output/LMMResults/ElevationCentroid_lmm_results.rds")
tss_resp_var_lmer <- readRDS("output/LMMResults/TSS_lmm_results.rds")

# Creates a plot for each Response Variable (will be 5)
Plot_Residual_Normality_of_Response_Variables(schoener_resp_var_lmer, 
                                              meansuit_resp_var_lmer, 
                                              suitarea_resp_var_lmer,
                                              elevcentroid_resp_var_lmer,
                                              tss_resp_var_lmer)

# # Correlation Matrices (helped to show what fixed variables were redundant)
# schoener_correlation_heat_map <- Create_Correlation_Heat_Map(schoener_resp_var_lmer, "Schoener's D")
# meansuit_correlation_heat_map <- Create_Correlation_Heat_Map(meansuit_resp_var_lmer, "Mean Suitablity")
# suitarea_correlation_heat_map <- Create_Correlation_Heat_Map(suitarea_resp_var_lmer, "Suitable Area")
# elevcentroid_resp_correlation_heat_map <- Create_Correlation_Heat_Map(elevcentroid_resp_var_lmer, "Elevation Centroid")
# tss_correlation_heat_map <- Create_Correlation_Heat_Map(tss_resp_var_lmer, "Model Performance (TSS)")

# grid.arrange(schoener_correlation_heat_map,
#              meansuit_correlation_heat_map,
#              suitarea_correlation_heat_map,
#              elevcentroid_resp_correlation_heat_map,
#              tss_correlation_heat_map,
#              ncol = 3, nrow = 2)


# opening R2Beta Results
schoener_var_r2 <- readRDS("output/R2BetaResults/schoener_var_r2.rds")
meansuit_var_r2 <- readRDS("output/R2BetaResults/meansuit_var_r2.rds")
suitarea_var_r2 <- readRDS("output/R2BetaResults/suitarea_var_r2.rds")
elevcentroid_var_r2 <- readRDS("output/R2BetaResults/elevcentroid_var_r2.rds")
tss_var_r2 <- readRDS("output/R2BetaResults/tss_var_r2.rds")

# make graphs (one with and one without interactions)
Plot_R2Beta_Results(schoener_var_r2, meansuit_var_r2, suitarea_var_r2, elevcentroid_var_r2, tss_var_r2)


# Open Estimated Marginal Trends Results
algorithm_emtrends_results <- read_tsv("output/EMTrendsResults/algorithm_emtrends_results.tsv")
paStrategy_emtrends_results <- read_tsv("output/EMTrendsResults/paStrategy_emtrends_results.tsv")
paNumber_emtrends_results <- read_tsv("output/EMTrendsResults/paNumber_emtrends_results.tsv")

algorithm_emtrends_results
# Plot Estimated Marginal Trends 
Plot_EMTrends_Results(algorithm_emtrends_results, "Algorithm") 
Plot_EMTrends_Results(paStrategy_emtrends_results, "PAStrategy") 
Plot_EMTrends_Results(paNumber_emtrends_results, "PANumber") 

