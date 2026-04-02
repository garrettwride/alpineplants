library(lme4)
library(r2glmm)
library(emmeans)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(partR2)

Run_LMM <- function(biomod_data, response_variable, output_file_path) {
  # for SPECIALISTS VS GENERALISTS as a whole (different Niche Scores): 
  # Random Effects = Species
  # Fixed Effects = NicheBreadth, Algorithm, PAStrategy, PANumber
  lmer_output <- lmer(as.formula(paste0(response_variable, " ~ 
                    NicheBreadth + Algorithm + PAStrategy + PANumber +
                    NicheBreadth:Algorithm + 
                    NicheBreadth:PAStrategy + 
                    NicheBreadth:PANumber + 
                    (1 | Species)")),
                      data = biomod_data,
                      REML = FALSE )
  
  # short hand for "NicheBreadth + Algorithm + NicheBreadth:Algorithm" = "NicheBreadth*Algorithm"
  
  # save results as txt file
  summary <- summary(lmer_output, correlation = TRUE)
  capture.output(summary, file = output_file_path)
  return(lmer_output)
}

Create_Correlation_Heat_Map <- function(lmer_results, fixed_var) {
  summary(lmer_results)
  correlation_matrix <- cov2cor(vcov(lmer_results))
  correlation_matrix <- as.matrix(correlation_matrix)
  correlation_long <- reshape2::melt(correlation_matrix)
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

Run_PartR2 <- function(lmer_results) {
  part_outcome <- partR2(
    lmer_results,
    partvars = c("NicheBreadth", "Algorithm", "PAStrategy", "PANumber"),
    R2_type = "conditional",                    # conditional or marginal??
    data = biomod_data,
    nboot = 50
  )
  
  return(part_outcome)
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

Compute_EMTrends <- function(fixed_var, schoener_lmer, meansuit_lmer, suitarea_lmer, elevcentroid_lmer, tss_lmer, output_file_path) {
  emt_schoener = emtrends(schoener_lmer, as.formula(paste0("~", fixed_var)), var = "NicheBreadth")
  emt_meansuit = emtrends(meansuit_lmer, as.formula(paste0("~", fixed_var)), var = "NicheBreadth")
  emt_suitarea = emtrends(suitarea_lmer, as.formula(paste0("~", fixed_var)), var = "NicheBreadth")
  emt_elevcentroid = emtrends(elevcentroid_lmer, as.formula(paste0("~", fixed_var)), var = "NicheBreadth")
  emt_tss = emtrends(tss_lmer, as.formula(paste0("~", fixed_var)), var = "NicheBreadth")
  
  
  schoener_df = as.data.frame(emt_schoener) %>%
    mutate(ResponseVariable = "Schoener's D")
  
  meansuit_df = as.data.frame(emt_meansuit) %>%
    mutate(ResponseVariable = "Mean Suitablilty")
  
  suitarea_df = as.data.frame(emt_suitarea) %>%
    mutate(ResponseVariable = "Suitable Area")
  
  elevcentroid_df = as.data.frame(emt_elevcentroid) %>%
    mutate(ResponseVariable = "Elevation Centroid")
  
  tss_df = as.data.frame(emt_tss) %>%
    mutate(ResponseVariable = "Model Performance (TSS)")
  
  
  combined_df = bind_rows(schoener_df, meansuit_df, suitarea_df, elevcentroid_df, tss_df)
  
  # save results
  write.csv(combined_df, file = output_file_path, row.names = FALSE)
  
  return(combined_df)
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


# csv of Biomod output 
# Change the file path to where the file is saved on your computer
biomod_results_file_path <- "./testFiles/SDM_results_for_LMM_3.csv"
biomod_data <- read.csv(biomod_results_file_path)

# makes PAStrategy and PANumber categorical variables 
biomod_data$PAStrategy <- as.factor(biomod_data$PAStrategy)
biomod_data$PANumber <- as.factor(biomod_data$PANumber)
biomod_data$Algorithm <- as.factor(biomod_data$Algorithm)
biomod_data$Species <- as.factor(biomod_data$Species)
biomod_data$Bootstrap <- as.factor(biomod_data$Bootstrap)

#run LMM 
schoener_resp_var_lmer <- Run_LMM(biomod_data, "SchoenersD", "LMMResults/SchoenersD_lmm_results.txt")
meansuit_resp_var_lmer <- Run_LMM(biomod_data, "MeanSuitability", "LMMResults/MeanSuitability_lmm_results.txt")
suitarea_resp_var_lmer <- Run_LMM(biomod_data, "SuitableArea", "LMMResults/SuitableArea_lmm_results.txt")
elevcentroid_resp_var_lmer <- Run_LMM(biomod_data, "ElevationCentroid", "LMMResults/ElevationCentroid_lmm_results.txt")
tss_resp_var_lmer <- Run_LMM(biomod_data, "TSS", "LMMResults/TSS_lmm_results.txt")
# vcov(schoener_resp_var_lmer)
# vcov(meansuit_resp_var_lmer)

# Correlation Matrices 
schoener_correlation_heat_map <- Create_Correlation_Heat_Map(schoener_resp_var_lmer, "Schoener's D")
meansuit_correlation_heat_map <- Create_Correlation_Heat_Map(meansuit_resp_var_lmer, "Mean Suitablity")
suitarea_correlation_heat_map <- Create_Correlation_Heat_Map(suitarea_resp_var_lmer, "Suitable Area")
elevcentroid_resp_correlation_heat_map <- Create_Correlation_Heat_Map(elevcentroid_resp_var_lmer, "Elevation Centroid")
tss_correlation_heat_map <- Create_Correlation_Heat_Map(tss_resp_var_lmer, "Model Performance (TSS)")

grid.arrange(schoener_correlation_heat_map, 
             meansuit_correlation_heat_map, 
             suitarea_correlation_heat_map, 
             elevcentroid_resp_correlation_heat_map, 
             tss_correlation_heat_map,
             ncol = 3, nrow = 2)


# calculation of PERCENT VARIANCE 
schoener_var_r2 <- r2beta(schoener_resp_var_lmer, partial = TRUE, method = "nsj")
meansuit_var_r2 <- r2beta(meansuit_resp_var_lmer, partial = TRUE, method = "nsj")
suitarea_var_r2 <- r2beta(suitarea_resp_var_lmer, partial = TRUE, method = "nsj")
elevcentroid_var_r2 <- r2beta(elevcentroid_resp_var_lmer, partial = TRUE, method = "nsj")
tss_var_r2 = r2beta(tss_resp_var_lmer, partial = TRUE, method = "nsj")

# plotting variances WITHOUT interactions
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

# plotting variances WITH interactions
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


# marginal
# partr2_schoener <- Run_PartR2(schoener_resp_var_lmer)
# #partr2_meansuit <- Run_PartR2(meansuit_resp_var_lmer)
# #partr2_suitarea <- Run_PartR2(suitarea_resp_var_lmer)
# #partr2_elevcentroid <- Run_PartR2(elevcentroid_resp_var_lmer)
# #partr2_tss <- Run_PartR2(tss_resp_var_lmer)
# 
# # We care about the "Inclusive R2 (SC^2 * R2):" section
# summary(partr2_schoener)
# sr2 <- partr2_schoener$IR2
# sr2

# 
# ggplot(sr2, aes(x = reorder(Predictor, IR2), y = IR2)) +
#   geom_col(fill = "darkgreen") +
#   geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
#   coord_flip() +
#   labs(
#     x = "Predictor",
#     y = "Inclusive R²",
#     title = "Total contribution of predictors (including shared variance)"
#   ) +
#   theme_minimal()
# # partr2_schoener$boot_warnings
# 
# 
# 







## LOOK AT THIS AND FIX

# Estimated Marginal Trends of linear trends 
algorithm_emtrends_results <- Compute_EMTrends("Algorithm", schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, tss_resp_var_lmer, "LMMResults/algorithm_emtrends_results.csv")
paStrategy_emtrends_results <- Compute_EMTrends("PAStrategy", schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, tss_resp_var_lmer, "LMMResults/paStrategy_emtrends_results.csv")
paNumber_emtrends_results <- Compute_EMTrends("PANumber", schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, tss_resp_var_lmer, "LMMResults/paNumber_emtrends_results.csv")


# Plot Results 
Plot_EMTrends_Results(algorithm_emtrends_results, "Algorithm") 
Plot_EMTrends_Results(paStrategy_emtrends_results, "PAStrategy") 
Plot_EMTrends_Results(paNumber_emtrends_results, "PANumber") 



