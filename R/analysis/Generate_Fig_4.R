library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(gridExtra)

biomod_data <- read.csv("./data/SDM_results_for_LMM_3.csv")

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

plot_schoenersD <- function(df, title) {
  ggplot(df, aes(x = NicheBreadth, y = sd_SchoenersD)) +
    geom_point(color = "black", alpha = 0.8, size = 2.0) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "blue") +
    labs(
      x = "Niche Breadth",
      y = "Uncertainty (SD of Schoener’s D overlap)",
      title = title
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold",size = 14)
    )
}

p_method <- plot_schoenersD(
  method_uncertainty,
  "A. Methodological Uncertainty"
)

p_stochastic <- plot_schoenersD(
  stochastic_uncertainty,
  "B. Stochastic Uncertainty"
) +
  labs(y = NULL)

grid.arrange(p_method, p_stochastic, ncol = 2)
