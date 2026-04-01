library(ggplot2)
library(gridExtra)


Make_Graph <- function(y_var, color_var) {
  my_plot <- ggplot(biomod_data, aes_string(x = "NicheBreadth", y = y_var,color = color_var)) +
    geom_point(size = 1, alpha = 0.1) +
    geom_smooth(size = 2, method = "lm", se = TRUE, linetype = "solid") +
    labs(
      x = "Niche Breadth",
      y = y_var,
      title = paste0(y_var, " vs Niche Breadth by ",color_var),
      color = color_var
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.position = "right"
    )
  
  return(my_plot)
}


# Change the file path to where BioMOD output csv is saved on your computer
biomod_results_file_path <- "./testFiles/SDM_results_for_LMM_3.csv"
biomod_data <- read.csv(biomod_results_file_path)

# makes PAStrategy and PANumber categorical variables 
biomod_data$PAStrategy <- as.factor(biomod_data$PAStrategy)
biomod_data$PANumber <- as.factor(biomod_data$PANumber)
biomod_data$Algorithm <- as.factor(biomod_data$Algorithm)
biomod_data$Species <- as.factor(biomod_data$Species)
biomod_data$Bootstrap <- as.factor(biomod_data$Bootstrap)


# Plot NicheBreadth and Schoener's D 
schoeners_algorithm_plot <- Make_Graph("SchoenersD", "Algorithm")
schoeners_paStrategy_plot <- Make_Graph("SchoenersD", "PAStrategy")
schoeners_paNumber_plot <- Make_Graph("SchoenersD", "PANumber")

grid.arrange(schoeners_algorithm_plot, 
             schoeners_paStrategy_plot, 
             schoeners_paNumber_plot,
             ncol = 3, nrow = 1)

# Plot NicheBreadth and MeanSuitability
meanSuit_algorithm_plot <- Make_Graph("MeanSuitability", "Algorithm")
meanSuit_paStrategy_plot <- Make_Graph("MeanSuitability", "PAStrategy")
meanSuit_paNumber_plot <- Make_Graph("MeanSuitability", "PANumber")

grid.arrange(meanSuit_algorithm_plot, 
             meanSuit_paStrategy_plot, 
             meanSuit_paNumber_plot,
             ncol = 3, nrow = 1)

# Plot NicheBreadth and SuitableArea
suitArea_algorithm_plot <- Make_Graph("SuitableArea", "Algorithm")
suitArea_paStrategy_plot <- Make_Graph("SuitableArea", "PAStrategy")
suitArea_paNumber_plot <- Make_Graph("SuitableArea", "PANumber")

grid.arrange(suitArea_algorithm_plot, 
             suitArea_paStrategy_plot, 
             suitArea_paNumber_plot,
             ncol = 3, nrow = 1)

# Plot NicheBreadth and ElevationCentroid
elevCentroid_algorithm_plot <- Make_Graph("ElevationCentroid", "Algorithm")
elevCentroid_paStrategy_plot <- Make_Graph("ElevationCentroid", "PAStrategy")
elevCentroid_paNumber_plot <- Make_Graph("ElevationCentroid", "PANumber")

grid.arrange(elevCentroid_algorithm_plot, 
             elevCentroid_paStrategy_plot, 
             elevCentroid_paNumber_plot,
             ncol = 3, nrow = 1)

# Plot NicheBreadth and TSS
tss_algorithm_plot <- Make_Graph("TSS", "Algorithm")
tss_paStrategy_plot <- Make_Graph("TSS", "PAStrategy")
tss_paNumber_plot <- Make_Graph("TSS", "PANumber")

grid.arrange(tss_algorithm_plot, 
             tss_paStrategy_plot, 
             tss_paNumber_plot,
             ncol = 3, nrow = 1)

