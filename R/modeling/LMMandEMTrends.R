library(lme4)
library(r2glmm)
library(emmeans)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(partR2)
library(readr)

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
  #summary <- summary(lmer_output, correlation = TRUE)
  saveRDS(lmer_output, file = output_file_path)
  return(lmer_output)
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


# csv of Biomod output 
# Change the file path to where the file is saved on your computer
biomod_results_file_path <- "./data/SDM_results_for_LMM_3.csv"
biomod_data <- read.csv(biomod_results_file_path)

# makes PAStrategy and PANumber categorical variables 
biomod_data$PAStrategy <- as.factor(biomod_data$PAStrategy)
biomod_data$PANumber <- as.factor(biomod_data$PANumber)
biomod_data$Algorithm <- as.factor(biomod_data$Algorithm)
biomod_data$Species <- as.factor(biomod_data$Species)
biomod_data$Bootstrap <- as.factor(biomod_data$Bootstrap)


#run LMM 
if (!dir.exists("./output/LMMResults")) {
  dir.create("./output/LMMResults")
}

schoener_resp_var_lmer <- Run_LMM(biomod_data, "SchoenersD", "output/LMMResults/SchoenersD_lmm_results.rds")
meansuit_resp_var_lmer <- Run_LMM(biomod_data, "MeanSuitability", "output/LMMResults/MeanSuitability_lmm_results.rds")
suitarea_resp_var_lmer <- Run_LMM(biomod_data, "SuitableArea", "output/LMMResults/SuitableArea_lmm_results.rds")
elevcentroid_resp_var_lmer <- Run_LMM(biomod_data, "ElevationCentroid", "output/LMMResults/ElevationCentroid_lmm_results.rds")
tss_resp_var_lmer <- Run_LMM(biomod_data, "TSS", "output/LMMResults/TSS_lmm_results.rds")
# vcov(schoener_resp_var_lmer)


# calculation of PERCENT VARIANCE and save results
if (!dir.exists("./output/R2BetaResults")) {
  dir.create("./output/R2BetaResults")
}

schoener_var_r2 <- r2beta(schoener_resp_var_lmer, partial = TRUE, method = "nsj")
saveRDS(schoener_var_r2, file = "output/R2BetaResults/schoener_var_r2.rds")

meansuit_var_r2 <- r2beta(meansuit_resp_var_lmer, partial = TRUE, method = "nsj")
saveRDS(meansuit_var_r2, file = "output/R2BetaResults/meansuit_var_r2.rds")

suitarea_var_r2 <- r2beta(suitarea_resp_var_lmer, partial = TRUE, method = "nsj")
saveRDS(suitarea_var_r2, file = "output/R2BetaResults/suitarea_var_r2.rds")

elevcentroid_var_r2 <- r2beta(elevcentroid_resp_var_lmer, partial = TRUE, method = "nsj")
saveRDS(elevcentroid_var_r2, file = "output/R2BetaResults/elevcentroid_var_r2.rds")

tss_var_r2 = r2beta(tss_resp_var_lmer, partial = TRUE, method = "nsj")
saveRDS(tss_var_r2, file = "output/R2BetaResults/tss_var_r2.rds")


# # We were going to use this but not anymore
# partr2_schoener <- Run_PartR2(schoener_resp_var_lmer)
# #partr2_meansuit <- Run_PartR2(meansuit_resp_var_lmer)
# #partr2_suitarea <- Run_PartR2(suitarea_resp_var_lmer)
# #partr2_elevcentroid <- Run_PartR2(elevcentroid_resp_var_lmer)
# #partr2_tss <- Run_PartR2(tss_resp_var_lmer)
# 
# # We care about the "Inclusive R2 (SC^2 * R2):" section
# summary(partr2_schoener)
# sr2 <- partr2_schoener$IR2


# Estimated Marginal Trends of linear trends 
if (!dir.exists("./output/EMTrendsResults")) {
  dir.create("./output/EMTrendsResults")
}

algorithm_emtrends_results <- Compute_EMTrends("Algorithm", schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, tss_resp_var_lmer, "output/LMMResults/algorithm_emtrends_results.csv")
write_tsv(algorithm_emtrends_results, file = "output/EMTrendsResults/algorithm_emtrends_results.tsv")

paStrategy_emtrends_results <- Compute_EMTrends("PAStrategy", schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, tss_resp_var_lmer, "output/LMMResults/paStrategy_emtrends_results.csv")
write_tsv(paStrategy_emtrends_results, file = "output/EMTrendsResults/paStrategy_emtrends_results.tsv")

paNumber_emtrends_results <- Compute_EMTrends("PANumber", schoener_resp_var_lmer, meansuit_resp_var_lmer, suitarea_resp_var_lmer, elevcentroid_resp_var_lmer, tss_resp_var_lmer, "output/LMMResults/paNumber_emtrends_results.csv")
write_tsv(paNumber_emtrends_results, file = "output/EMTrendsResults/paNumber_emtrends_results.tsv")


# Estimated Marginal Means
emm_schoener = emmeans(schoener_resp_var_lmer, ~NicheBreadth)
emm_meansuit = emmeans(meansuit_resp_var_lmer, ~NicheBreadth)
emm_suitarea = emmeans(suitarea_resp_var_lmer, ~NicheBreadth)
emm_elevcentroid = emmeans(elevcentroid_resp_var_lmer, ~NicheBreadth)
emm_tss = emmeans(tss_resp_var_lmer, ~NicheBreadth)

schoener_df = as.data.frame(emm_schoener) %>%
  mutate(ResponseVariable = "Schoener's D")

meansuit_df = as.data.frame(emm_meansuit) %>%
  mutate(ResponseVariable = "Mean Suitablilty")

suitarea_df = as.data.frame(emm_suitarea) %>%
  mutate(ResponseVariable = "Suitable Area")

elevcentroid_df = as.data.frame(emm_elevcentroid) %>%
  mutate(ResponseVariable = "Elevation Centroid")

tss_df = as.data.frame(emm_tss) %>%
  mutate(ResponseVariable = "TSS")

combined_emm_df = bind_rows(schoener_df, meansuit_df, suitarea_df, elevcentroid_df, tss_df)
write_tsv(combined_emm_df, file = "output/EMMeansCombinedResults.tsv")

# combined_emm_df

# ggplot(combined_emm_df, aes(x = NicheBreadth, y = emmean)) +
#   geom_bar(stat = "identity",color = "black") +
#   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.5) +
#   facet_wrap(~ResponseVariable, scales = "free") +
#   labs(y = "Response Variable Predictions", x = "Species Type")

