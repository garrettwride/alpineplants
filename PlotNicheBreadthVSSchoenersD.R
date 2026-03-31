library(ggplot2)
library(gridExtra)

# Change the file path to where BioMOD output csv is saved on your computer
biomod_results_file_path <- "./testFiles/SDM_results_for_LMM_3.csv"
biomod_data <- read.csv(biomod_results_file_path)

# makes PAStrategy and PANumber categorical variables 
biomod_data$PAStrategy <- as.factor(biomod_data$PAStrategy)
biomod_data$PANumber <- as.factor(biomod_data$PANumber)
biomod_data$Algorithm <- as.factor(biomod_data$Algorithm)
biomod_data$Species <- as.factor(biomod_data$Species)
biomod_data$Bootstrap <- as.factor(biomod_data$Bootstrap)


# Plot NicheBreadth and Schoener's D based on Algorithm
algorithm_plot <- ggplot(biomod_data, aes(x = NicheBreadth, y = SchoenersD, color = Algorithm)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(size = 2, method = "lm", se = TRUE, linetype = "solid") +
  labs(
    x = "Niche Breadth",
    y = "Schoener's D",
    title = "Schoener's D vs Niche Breadth by Algorithm",
    color = "Algorithm"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.position = "right"
  )


# Plot NicheBreadth and Schoener's D based on PAStrategy
paStrategy_plot <- ggplot(biomod_data, aes(x = NicheBreadth, y = SchoenersD, color = PAStrategy)) +
  geom_point(size = 1, alpha = 0.8) +   
  geom_smooth(size = 2, method = "lm", se = TRUE, linetype = "solid") +
  labs(
    x = "Niche Breadth",
    y = "Schoener's D",
    title = "Schoener's D vs Niche Breadth by PAStrategy",
    color = "PAStrategy"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.position = "right"
  )


# Plot NicheBreadth and Schoener's D based on PANumber
paNumber_plot <- ggplot(biomod_data, aes(x = NicheBreadth, y = SchoenersD, color = PANumber)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(size = 2, method = "lm", se = TRUE, linetype = "solid") +
  labs(
    x = "Niche Breadth",
    y = "Schoener's D",
    title = "Schoener's D vs Niche Breadth by PANumber",
    color = "PANumber"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.position = "right"
  )


grid.arrange(algorithm_plot, 
             paStrategy_plot, 
             paNumber_plot,
             ncol = 3, nrow = 1)
