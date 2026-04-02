library(ggplot2)
library(gridExtra)


Make_Graph <- function(y_var, color_var, color_palette) {
  my_plot <- ggplot(biomod_data, aes_string(x = "NicheBreadth", y = y_var,color = color_var)) +
    geom_point(size = 1, alpha = 0.1) +
    geom_smooth(size = 2, method = "lm", se = TRUE, linetype = "solid") +
    labs(
      x = "Niche Breadth",
      y = y_var,
      title = paste0(y_var, " vs Niche Breadth by ",color_var),
      color = color_var
    ) +
    scale_color_manual(values = color_palette) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
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


algorithm_palette <- c("#F0E442", "#E69F00", "#56B4E9", "#000000")
paStrategy_palette <- c("#009E73", "#CC79A7")
paNumber_palette <- c("#D55E00", "#0072B2", "#999999")


# Plot NicheBreadth and Schoener's D 
schoeners_algorithm_plot <- Make_Graph("SchoenersD", "Algorithm", algorithm_palette)
schoeners_paStrategy_plot <- Make_Graph("SchoenersD", "PAStrategy", paStrategy_palette)
schoeners_paNumber_plot <- Make_Graph("SchoenersD", "PANumber", paNumber_palette)

grid.arrange(schoeners_algorithm_plot, 
             schoeners_paStrategy_plot, 
             schoeners_paNumber_plot,
             ncol = 3, nrow = 1)

# Plot NicheBreadth and MeanSuitability
meanSuit_algorithm_plot <- Make_Graph("MeanSuitability", "Algorithm", algorithm_palette)
meanSuit_paStrategy_plot <- Make_Graph("MeanSuitability", "PAStrategy", paStrategy_palette)
meanSuit_paNumber_plot <- Make_Graph("MeanSuitability", "PANumber", paNumber_palette)

grid.arrange(meanSuit_algorithm_plot, 
             meanSuit_paStrategy_plot, 
             meanSuit_paNumber_plot,
             ncol = 3, nrow = 1)

# Plot NicheBreadth and SuitableArea
suitArea_algorithm_plot <- Make_Graph("SuitableArea", "Algorithm", algorithm_palette)
suitArea_paStrategy_plot <- Make_Graph("SuitableArea", "PAStrategy", paStrategy_palette)
suitArea_paNumber_plot <- Make_Graph("SuitableArea", "PANumber", paNumber_palette)

grid.arrange(suitArea_algorithm_plot, 
             suitArea_paStrategy_plot, 
             suitArea_paNumber_plot,
             ncol = 3, nrow = 1)

# Plot NicheBreadth and ElevationCentroid
elevCentroid_algorithm_plot <- Make_Graph("ElevationCentroid", "Algorithm", algorithm_palette)
elevCentroid_paStrategy_plot <- Make_Graph("ElevationCentroid", "PAStrategy", paStrategy_palette)
elevCentroid_paNumber_plot <- Make_Graph("ElevationCentroid", "PANumber", paNumber_palette)

grid.arrange(elevCentroid_algorithm_plot, 
             elevCentroid_paStrategy_plot, 
             elevCentroid_paNumber_plot,
             ncol = 3, nrow = 1)

# Plot NicheBreadth and TSS
tss_algorithm_plot <- Make_Graph("TSS", "Algorithm", algorithm_palette)
tss_paStrategy_plot <- Make_Graph("TSS", "PAStrategy", paStrategy_palette)
tss_paNumber_plot <- Make_Graph("TSS", "PANumber", paNumber_palette)

grid.arrange(tss_algorithm_plot, 
             tss_paStrategy_plot, 
             tss_paNumber_plot,
             ncol = 3, nrow = 1)

