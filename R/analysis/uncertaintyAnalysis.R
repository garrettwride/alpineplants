library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

biomod_data <- read.csv("./data/SDM_results_for_LMM.csv")

# Ensure factors
biomod_data <- biomod_data %>%
  mutate(
    Species = as.factor(Species),
    Algorithm = as.factor(Algorithm),
    PAStrategy = as.factor(PAStrategy),
    PANumber = as.factor(PANumber),
    Bootstrap = as.factor(Bootstrap)
  )

nb_lookup <- biomod_data %>%
  select(Species, NicheBreadth) %>%
  distinct()

stochastic_uncertainty <- biomod_data %>%
  group_by(Species, Algorithm, PAStrategy, PANumber) %>%
  summarize(
    sd_SchoenersD = sd(SchoenersD),
    sd_MeanSuitability = sd(MeanSuitability),
    sd_SuitableArea = sd(SuitableArea),
    sd_ElevationCentroid = sd(ElevationCentroid),
    sd_TSS = sd(TSS),
    .groups = "drop"
  ) %>%
  left_join(nb_lookup, by = "Species")

biomod_summary <- biomod_data %>%
  group_by(Species, Algorithm, PAStrategy, PANumber) %>%
  summarize(
    SchoenersD = mean(SchoenersD),
    MeanSuitability = mean(MeanSuitability),
    SuitableArea = mean(SuitableArea),
    ElevationCentroid = mean(ElevationCentroid),
    TSS = mean(TSS),
    .groups = "drop"
  )

method_uncertainty <- biomod_summary %>%
  group_by(Species) %>%
  summarize(
    sd_SchoenersD = sd(SchoenersD),
    sd_MeanSuitability = sd(MeanSuitability),
    sd_SuitableArea = sd(SuitableArea),
    sd_ElevationCentroid = sd(ElevationCentroid),
    sd_TSS = sd(TSS),
    .groups = "drop"
  ) %>%
  left_join(nb_lookup, by = "Species")

run_lm_models <- function(df) {
  responses <- c(
    "sd_SchoenersD",
    "sd_MeanSuitability",
    "sd_SuitableArea",
    "sd_ElevationCentroid",
    "sd_TSS"
  )
  
  models <- map(responses, ~ lm(as.formula(paste(.x, "~ NicheBreadth")), data = df))
  names(models) <- responses
  
  return(models)
}

lm_results <- run_lm_models(method_uncertainty)
lm_results_stochastic <- run_lm_models(stochastic_uncertainty)

# Print summaries
map(lm_results, summary)
map(lm_results_stochastic, summary)

plot_uncertainty <- function(df, title) {
  df %>%
    pivot_longer(
      cols = starts_with("sd_"),
      names_to = "Response",
      values_to = "Uncertainty"
    ) %>%
    ggplot(aes(x = NicheBreadth, y = Uncertainty)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~ Response, scales = "free_y") +
    labs(
      x = "Niche Breadth",
      y = "Uncertainty (SD)",
      title = title
    ) +
    theme_minimal()
}

plot_uncertainty(method_uncertainty, "Methodological Uncertainty")
plot_uncertainty(stochastic_uncertainty, "Stochastic Uncertainty")

uncertainty_by_alg <- biomod_summary %>%
  group_by(Species, Algorithm) %>%
  summarize(across(
    c(SchoenersD, MeanSuitability, SuitableArea, ElevationCentroid, TSS),
    sd,
    .names = "sd_{.col}"
  ), .groups = "drop") %>%
  left_join(nb_lookup, by = "Species")

uncertainty_by_pa <- biomod_summary %>%
  group_by(Species, PAStrategy) %>%
  summarize(across(
    c(SchoenersD, MeanSuitability, SuitableArea, ElevationCentroid, TSS),
    sd,
    .names = "sd_{.col}"
  ), .groups = "drop") %>%
  left_join(nb_lookup, by = "Species")

uncertainty_by_pan <- biomod_summary %>%
  group_by(Species, PANumber) %>%
  summarize(across(
    c(SchoenersD, MeanSuitability, SuitableArea, ElevationCentroid, TSS),
    sd,
    .names = "sd_{.col}"
  ), .groups = "drop") %>%
  left_join(nb_lookup, by = "Species")

lmer(sd_SchoenersD ~ NicheBreadth * Algorithm + (1 | Species), data = uncertainty_by_alg)
lmer(sd_SchoenersD ~ NicheBreadth * PAStrategy + (1 | Species), data = uncertainty_by_pa)
lmer(sd_SchoenersD ~ NicheBreadth * PANumber + (1 | Species), data = uncertainty_by_pan)

lmer(sd_MeanSuitability ~ NicheBreadth * Algorithm + (1 | Species), data = uncertainty_by_alg)
lmer(sd_MeanSuitability ~ NicheBreadth * PAStrategy + (1 | Species), data = uncertainty_by_pa)
lmer(sd_MeanSuitability ~ NicheBreadth * PANumber + (1 | Species), data = uncertainty_by_pan)

lmer(sd_SuitableArea ~ NicheBreadth * Algorithm + (1 | Species), data = uncertainty_by_alg)
lmer(sd_SuitableArea ~ NicheBreadth * PAStrategy + (1 | Species), data = uncertainty_by_pa)
lmer(sd_SuitableArea ~ NicheBreadth * PANumber + (1 | Species), data = uncertainty_by_pan)

lmer(sd_ElevationCentroid ~ NicheBreadth * Algorithm + (1 | Species), data = uncertainty_by_alg)
lmer(sd_ElevationCentroid ~ NicheBreadth * PAStrategy + (1 | Species), data = uncertainty_by_pa)
lmer(sd_ElevationCentroid ~ NicheBreadth * PANumber + (1 | Species), data = uncertainty_by_pan)

lmer(sd_TSS ~ NicheBreadth * Algorithm + (1 | Species), data = uncertainty_by_alg)
lmer(sd_TSS ~ NicheBreadth * PAStrategy + (1 | Species), data = uncertainty_by_pa)
lmer(sd_TSS ~ NicheBreadth * PANumber + (1 | Species), data = uncertainty_by_pan)
